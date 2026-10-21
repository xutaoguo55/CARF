"""Fixed-scope factorial measurement-diagnostic feasibility pilot.

Binomial thinning is conditional on observed counts, not a new experiment.
Dm/Dc cells are disjoint; historical feature selection nevertheless used all D.
Held-out responses and the historical upstream all-cell HVG scope stay fixed.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix, vstack
from scipy.stats import spearmanr

from carf_core import residualize, topk_weights
from io_utils import HERE, ROOT, column, gene_symbols, input_path, sha256, write_json

CONFIG = HERE / "config/measurement_diagnostic_pilot_v1.json"
OUT = HERE / "results/measurement_diagnostic_pilot_v1"
RUNNER = HERE / "run_measurement_model.R"


def check_hash(path, expected):
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"Hash mismatch: {path}")
    return actual


def read_sparse_rows(path, rows):
    """Read retained raw counts in bounded slices; preserve requested row order."""
    rows = np.asarray(rows, dtype=int)
    if len(np.unique(rows)) != len(rows):
        raise ValueError("Repeated cell rows")
    order = np.argsort(rows)
    ordered = rows[order]
    blocks = []
    with h5py.File(path) as f:
        g = f["layers/counts"]
        shape = tuple(int(x) for x in g.attrs["shape"])
        ptr = g["indptr"][:]
        for start in range(0, len(rows), 64):
            wanted = ordered[start:start+64]
            first, last = int(wanted[0]), int(wanted[-1]) + 1
            a, b = int(ptr[first]), int(ptr[last])
            block = csr_matrix((g["data"][a:b], g["indices"][a:b], ptr[first:last+1]-a),
                               shape=(last-first, shape[1]))
            blocks.append(block[wanted-first])
    result = vstack(blocks, format="csr")[np.argsort(order)]
    if not np.isfinite(result.data).all() or (result.data < 0).any() or not np.equal(result.data, np.floor(result.data)).all():
        raise ValueError("Noninteger or invalid source UMI")
    result.data = result.data.astype(np.int64)
    return result


def stream_rng(dataset, q, seed, stream):
    payload = f"measurement_v1|{dataset}|{q}|{seed}|{stream}".encode()
    return np.random.default_rng(int.from_bytes(hashlib.sha256(payload).digest()[:16], "little"))


def thin_counts(counts, q, rng):
    if not 0 < q <= 1:
        raise ValueError("Invalid thinning probability")
    out = counts.copy().tocsr()
    if q < 1:
        out.data = rng.binomial(out.data.astype(np.int64), q)
        out.eliminate_zeros()
    return out


def covariate_mean(counts):
    totals = np.asarray(counts.sum(axis=1)).ravel()
    if (totals <= 0).any():
        raise ValueError("Zero-library Dc cell; protocol forbids silent selection")
    normalized = counts.astype(float).multiply(10000.0/totals[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    return np.asarray(normalized.mean(axis=0)).ravel(), totals


def ranking_scores(distance, cov):
    if not np.isfinite(distance).all() or (distance < 0).any():
        raise ValueError("Distances must be finite and nonnegative")
    logscore = np.log10(np.maximum(distance, 1e-10))
    return {"raw_distance": distance, "EDR": np.vstack([residualize(x, cov) for x in logscore])}, logscore


def overlap(weights, baseline_weights):
    return np.minimum(weights, baseline_weights).sum(axis=1)/weights.sum(axis=1)


def evaluate_scores(scores, cov, logscore, B, C, null_response, baseline_weights, targets, k):
    weights = np.vstack([topk_weights(s, k) for s in scores])
    utility = (B*C) @ weights.T / k
    same = np.diag(utility)
    other = (utility.sum(axis=1)-same)/(len(targets)-1)
    null_utility = null_response @ weights.T/k
    ov = overlap(weights, baseline_weights)
    rows = []
    for i, target in enumerate(targets):
        rho = float(spearmanr(scores[i], cov).statistic) if np.ptp(scores[i]) > 0 and np.ptp(cov)>0 else np.nan
        sst = float(np.square(logscore[i]-logscore[i].mean()).sum())
        residual = residualize(logscore[i], cov)
        r2 = 1-float(np.square(residual).sum())/sst if sst>1e-20 else np.nan
        rows.append(dict(target=target,top50_overlap=ov[i],score_expression_rho=rho,
                         logscore_expression_r2=r2,total_BC=same[i],other_target_BC=other[i],
                         identity_gain=same[i]-other[i],distance_floor_count=int((logscore[i]<=-10).sum()),
                         null_response_mean=float(null_utility[:,i].mean()),
                         null_response_q95=float(np.quantile(null_utility[:,i],.95)),
                         null_response_abs_mean=float(np.abs(null_utility[:,i]).mean())))
    return pd.DataFrame(rows), weights, utility, null_utility


def prepare(dataset, cfg):
    print(f"Preparing verified fixed inputs: {dataset}", flush=True)
    old = HERE / "results/model_rankings_sctenifoldknk_pilot" / dataset
    prep_path, model_path = old/"preparation_manifest.json", old/"model_manifest.json"
    prep, model = json.loads(prep_path.read_text()), json.loads(model_path.read_text())
    check_hash(prep_path, model["input_manifest_sha256"])
    for name, expected in prep["input_sha256"].items():
        check_hash(old/name, expected)
    check_hash(old/"scTenifoldKnk_heat_distances.csv",model["output_sha256"]["distances"])
    source = input_path(dataset)
    check_hash(source, prep["source_sha256"])
    responses = HERE / "results/library_calibration_v2" / f"{dataset}_disjoint_library_pools_seed{cfg['base_response_seed']}"
    check_hash(responses/"manifest.json", prep["response_manifest_sha256"])
    receipt = json.loads((responses/"manifest.json").read_text())
    for name in ["responses.npz","cell_splits.csv.gz"]:
        check_hash(responses/name,receipt["outputs"][name])
    split = pd.read_csv(responses/"cell_splits.csv.gz")
    d = split.loc[split.arm=="covariates_D",["cell_index","gemgroup"]].drop_duplicates()
    dm = pd.read_csv(old/"training_cells.csv")
    dc = d.loc[~d.cell_index.isin(dm.cell_index)].copy()
    if len(dm)!=400 or not len(dc) or not set(dm.cell_index).issubset(d.cell_index) or set(dm.cell_index)&set(dc.cell_index):
        raise ValueError("Invalid Dm/Dc partition")
    forbidden = set(split.loc[split.arm!="covariates_D","cell_index"])
    if (set(dm.cell_index)|set(dc.cell_index)) & forbidden:
        raise ValueError("D cells overlap A/B/C")
    genes = (old/"genes.txt").read_text().splitlines()
    targets = (old/"targets.txt").read_text().splitlines()
    all_targets = set((old/"all_real_targets.txt").read_text().splitlines())
    z = np.load(responses/"responses.npz",allow_pickle=False)
    all_genes = list(z["genes"].astype(str))
    tasks = list(z["tasks"].astype(str))
    with h5py.File(source) as f:
        if list(gene_symbols(f).astype(str)) != all_genes:
            raise ValueError("Frozen response genes changed")
        labels = column(f["obs"],"perturbation")
        if not np.all(labels[np.r_[dm.cell_index,dc.cell_index]]=="control"):
            raise ValueError("Noncontrol cells in D")
    selected = np.array([all_genes.index(g) for g in genes])
    model_counts = mmread(old/"control_counts.mtx").tocsr().astype(np.int64)
    current_dm = read_sparse_rows(source,dm.cell_index.to_numpy(int))[:,selected].T.tocsr()
    if model_counts.shape!=current_dm.shape or (model_counts!=current_dm).nnz:
        raise ValueError("Archived Dm counts do not match source cells/genes")
    dc_counts = read_sparse_rows(source,dc.cell_index.to_numpy(int))
    common = np.array([g not in all_targets for g in genes])
    common_genes = np.asarray(genes)[common]
    response_cols = np.array([all_genes.index(g) for g in common_genes])
    response_rows = np.array([tasks.index(t) for t in targets])
    null_rows = np.array([i for i,t in enumerate(tasks) if t.startswith("control_vs_control_")])
    data = dict(dataset=dataset,old=old,genes=genes,targets=targets,common=common,
        common_genes=common_genes,model_counts=model_counts,dc_counts=dc_counts,selected=selected,
        A=z["A"][response_rows][:,response_cols],B=z["B"][response_rows][:,response_cols],
        C=z["C"][response_rows][:,response_cols],
        null=z["B"][null_rows][:,response_cols]*z["C"][null_rows][:,response_cols])
    folder=OUT/dataset
    folder.mkdir(parents=True,exist_ok=True)
    dm.assign(role="model_Dm").to_csv(folder/"model_cells.csv",index=False)
    dc.assign(role="covariate_Dc").to_csv(folder/"covariate_cells.csv",index=False)
    pd.DataFrame({"gene_symbol":common_genes}).to_csv(folder/"fixed_candidate_genes.csv",index=False)
    write_json(folder/"source_manifest.json",dict(dataset=dataset,model_cells=len(dm),covariate_cells=len(dc),
        common_genes=len(common_genes),targets=len(targets),source_sha256=prep["source_sha256"],
        original_preparation_sha256=sha256(prep_path),original_model_sha256=sha256(model_path),
        response_manifest_sha256=sha256(responses/"manifest.json"),responses_sha256=sha256(responses/"responses.npz"),
        Dm_Dc_cell_disjoint=True,D_disjoint_ABC=True,
        additional_limitation="Historical model-feature and target selection used all D, including this Dc. Cell-disjointness does not remove that fixed-feature dependence."))
    return data


def native_fit(data,q,thin_seed,model_seed,cfg):
    name=f"q{q:g}_thin{thin_seed}_model{model_seed}"
    folder=OUT/data["dataset"]/"fits"/name
    folder.mkdir(parents=True,exist_ok=True)
    counts=thin_counts(data["model_counts"],q,stream_rng(data["dataset"],q,thin_seed,"Dm"))
    meta=folder/"model_manifest.json"
    if meta.exists():
        record=json.loads(meta.read_text())
        check_hash(folder/"input_manifest.json",record["input_manifest_sha256"])
        check_hash(folder/"distances.csv",record["distances_sha256"])
        previous=json.loads((folder/"input_manifest.json").read_text())
        if previous["protocol_sha256"]!=sha256(CONFIG) or previous["runner_sha256"]!=sha256(RUNNER):
            raise ValueError("Completed fit from changed code/config")
        for filename,expected in previous["input_sha256"].items():
            check_hash(folder/filename,expected)
    else:
        mmwrite(folder/"control_counts.mtx",counts)
        for filename in ["genes.txt","targets.txt"]:
            shutil.copyfile(data["old"]/filename,folder/filename)
        write_json(folder/"input_manifest.json",dict(dataset=data["dataset"],q=q,thinning_seed=thin_seed,
            model_seed=model_seed,protocol_sha256=sha256(CONFIG),runner_sha256=sha256(RUNNER),
            input_umi_mass=int(counts.sum()),input_shape=list(counts.shape),
            input_sha256={p:sha256(folder/p) for p in ["control_counts.mtx","genes.txt","targets.txt"]}))
        print(f"Fitting {data['dataset']} {name} sequentially",flush=True)
        with (folder/"fit.log").open("w") as log:
            result=subprocess.run(["Rscript",str(RUNNER),str(folder)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Native fit failed: {folder/'fit.log'}")
    frame=pd.read_csv(folder/"distances.csv",index_col=0)
    if set(frame.index)!=set(data["targets"]) or set(frame.columns)!=set(data["genes"]):
        raise ValueError("Native result identity mismatch")
    return frame.loc[data["targets"],data["genes"]].to_numpy(float)


def dataset_run(dataset,cfg):
    data=prepare(dataset,cfg)
    rows=[]
    null_rows=[]
    c0_all,totals=covariate_mean(data["dc_counts"])
    c0=c0_all[data["selected"]][data["common"]]
    pd.DataFrame({"gene_symbol":data["common_genes"],"C0_mean":c0}).to_csv(OUT/dataset/"baseline_covariates.csv",index=False)
    s0_full=native_fit(data,1.0,0,cfg["primary_model_seed"],cfg)
    archived=pd.read_csv(data["old"]/"scTenifoldKnk_heat_distances.csv",index_col=0).loc[data["targets"],data["genes"]].to_numpy(float)
    equal=np.allclose(s0_full,archived,rtol=1e-10,atol=1e-12)
    write_json(OUT/dataset/"unthinned_refit_check.json",dict(matches_archived_within_tolerance=bool(equal),rtol=1e-10,atol=1e-12,
        max_absolute_difference=float(np.max(np.abs(s0_full-archived)))))
    if not equal:
        raise ValueError(f"{dataset} q1 native refit does not reproduce archived scores")
    s0=s0_full[:,data["common"]]
    base_methods,_=ranking_scores(s0,c0)
    base_weights={m:np.vstack([topk_weights(x,cfg["primary_k"]) for x in s]) for m,s in base_methods.items()}

    def evaluate(s,c,q,ts,ms,factor,kind):
        methods,logs=ranking_scores(s,c)
        for method,scores in methods.items():
            result,w,utility,null=evaluate_scores(scores,c,logs,data["B"],data["C"],data["null"],base_weights[method],data["targets"],cfg["primary_k"])
            context=dict(dataset=dataset,q=q,thinning_seed=ts,model_seed=ms,factor=factor,kind=kind,method=method)
            for key,value in context.items(): result[key]=value
            rows.append(result)
            for ni in range(null.shape[0]):
                null_rows.append(dict(**context,null_partition=ni,mean_over_target_rankings=float(null[ni].mean()),
                    q95_over_target_rankings=float(np.quantile(null[ni],.95))))
            saved=OUT/dataset/"evaluations"/f"{kind}_q{q:g}_thin{ts}_model{ms}_{factor}_{method}.npz"
            saved.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(saved,targets=np.asarray(data["targets"]),genes=data["common_genes"],scores=scores,
                covariate=c,topk_weights=w,utility_query_source=utility,null_utility=null)

    evaluate(s0,c0,1.0,0,42,"S0_C0","baseline")
    # Target-invariant expression and empirical A anchor: check identity algebra.
    anchors=[]
    for name,scores in [("control_expression",np.tile(c0,(len(data["targets"]),1))),
                        ("empirical_A_repeatability",np.abs(data["A"]))]:
        bw=np.vstack([topk_weights(x,cfg["primary_k"]) for x in scores])
        result,_,_,_=evaluate_scores(scores,c0,np.log10(np.maximum(np.abs(scores),1e-10)),data["B"],data["C"],data["null"],bw,data["targets"],cfg["primary_k"])
        result["anchor"]=name
        anchors.append(result)
        if name=="control_expression" and result.identity_gain.abs().max()>1e-12:
            raise AssertionError("Target-invariant expression has nonzero identity gain")
    pd.concat(anchors).to_csv(OUT/dataset/"anchors.csv",index=False)
    for ms in cfg["unthinned_model_seed_controls"]:
        ss=native_fit(data,1.0,0,ms,cfg)[:,data["common"]]
        evaluate(ss,c0,1.0,0,ms,"Sq_C0","algorithm_seed_control")
    cov_checks=[]
    for q in cfg["thinning_probabilities"]:
        if q==1: continue
        for ts in cfg["thinning_seeds"]:
            ct=thin_counts(data["dc_counts"],q,stream_rng(dataset,q,ts,"Dc"))
            cq_all,ctot=covariate_mean(ct)
            cq=cq_all[data["selected"]][data["common"]]
            cov_checks.append(dict(q=q,seed=ts,cells=len(ctot),zero_library_cells=int((ctot==0).sum()),
                UMI_mass=int(ctot.sum()),median_library=float(np.median(ctot)),
                Cq_C0_spearman=float(spearmanr(cq,c0).statistic),Cq_C0_mean_absolute_difference=float(np.abs(cq-c0).mean())))
            sq=native_fit(data,q,ts,42,cfg)[:,data["common"]]
            for factor,s,c in [("S0_C0",s0,c0),("S0_Cq",s0,cq),("Sq_C0",sq,c0),("Sq_Cq",sq,cq)]:
                evaluate(s,c,q,ts,42,factor,"factorial")
            print(f"Evaluated {dataset} q={q} thinning_seed={ts}: all four arms",flush=True)
    pd.DataFrame(cov_checks).to_csv(OUT/dataset/"covariate_thinning_checks.csv",index=False)
    df=pd.concat(rows,ignore_index=True)
    df.to_csv(OUT/dataset/"per_target_all_conditions.csv.gz",index=False)
    pd.DataFrame(null_rows).to_csv(OUT/dataset/"control_null_partitions.csv.gz",index=False)
    return df


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--dataset",choices=["Adamson","Norman"])
    args=parser.parse_args()
    cfg=json.loads(CONFIG.read_text())
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/"manifest.json").exists(): raise SystemExit("Completed pilot exists; refusing overwrite")
    code={p.name:sha256(p) for p in [Path(__file__),RUNNER,HERE/"carf_core.py",HERE/"io_utils.py"]}
    start=datetime.now(timezone.utc).isoformat()
    try:
        for dataset in ([args.dataset] if args.dataset else cfg["datasets"]):
            dataset_run(dataset,cfg)
        all_paths=[OUT/d/"per_target_all_conditions.csv.gz" for d in cfg["datasets"]]
        if not all(p.exists() for p in all_paths):
            print("Requested dataset complete; full pilot awaits remaining dataset",flush=True)
            return
        df=pd.concat([pd.read_csv(p) for p in all_paths],ignore_index=True)
        metrics=["top50_overlap","score_expression_rho","logscore_expression_r2","total_BC","identity_gain","null_response_mean","null_response_abs_mean"]
        seed_keys=["dataset","kind","q","factor","method","target"]
        averaged=df.groupby(seed_keys,dropna=False)[metrics].mean().reset_index()
        averaged.to_csv(OUT/"seed_averaged_per_target.csv.gz",index=False)
        summary=averaged.groupby(seed_keys[:-1],dropna=False)[metrics].agg(["mean","median"])
        summary.columns=["_".join(x) for x in summary.columns]
        summary.reset_index().to_csv(OUT/"condition_summary.csv",index=False)
        # Pair every condition with the original unthinned same-method target.
        base=df[df.kind=="baseline"][["dataset","method","target"]+metrics]
        paired=df.merge(base,on=["dataset","method","target"],suffixes=("","_baseline"),validate="many_to_one")
        for m in metrics: paired[m+"_change"]=paired[m]-paired[m+"_baseline"]
        paired.to_csv(OUT/"paired_changes_all_conditions.csv.gz",index=False)
        changes=[m+"_change" for m in metrics]
        pc=paired.groupby(seed_keys,dropna=False)[changes].mean().reset_index()
        ps=pc.groupby(seed_keys[:-1],dropna=False)[changes].agg(["mean","median"])
        ps.columns=["_".join(x) for x in ps.columns]
        ps.reset_index().to_csv(OUT/"paired_change_summary.csv",index=False)
        for p in [Path(__file__),RUNNER,HERE/"carf_core.py",HERE/"io_utils.py"]:
            if sha256(p)!=code[p.name]: raise RuntimeError("Code changed during analysis")
        outputs={str(p.relative_to(OUT)):sha256(p) for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="manifest.json"}
        write_json(OUT/"manifest.json",dict(status="completed_development_feasibility_not_submission_validation",started_utc=start,
            finished_utc=datetime.now(timezone.utc).isoformat(),protocol=cfg,protocol_sha256=sha256(CONFIG),code_sha256=code,
            output_sha256=outputs,additional_limitation="Historical D-based feature selection includes Dc; only cells and thinning RNG streams are disjoint.",
            uncertainty="All summaries descriptive; seeds averaged within target, no biological CIs or P values."))
        print("Completed measurement diagnostic pilot; all prespecified conditions retained",flush=True)
    except Exception as exc:
        write_json(OUT/"failure.json",dict(time_utc=datetime.now(timezone.utc).isoformat(),error=repr(exc),protocol_sha256=sha256(CONFIG),code_sha256=code))
        raise


if __name__=="__main__": main()
