"""Transparent control-only association comparator on frozen measurement inputs.

This baseline estimates association, not a causal perturbation effect. It uses
no new cells, targets, genes, hyperparameters, or molecular thinning draws.
"""
from __future__ import annotations
from datetime import datetime, timezone
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
from scipy.io import mmread

from carf_core import residualize, topk_weights
from io_utils import HERE, sha256, write_json
from measurement_diagnostic_pilot import evaluate_scores

CONFIG = HERE / "config/pearson_measurement_comparator_v1.json"
SOURCE = HERE / "results/measurement_diagnostic_pilot_v1"
OUT = HERE / "results/pearson_measurement_comparator_v1"


def normalized_model_expression(counts):
    """Normalize cells across the entire fixed model gene universe."""
    counts = np.asarray(counts, float)
    if counts.ndim != 2 or not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all():
        raise ValueError("Expected finite nonnegative integer gene-by-cell counts")
    totals = counts.sum(axis=0)
    if np.any(totals <= 0):
        raise ValueError("Zero-library model cell; protocol forbids reselection")
    return np.log1p(counts / totals[None, :] * 10000), totals


def absolute_pearson(expression, target_indices):
    """Return target-by-gene scores plus explicit zero-variance indicators."""
    expression = np.asarray(expression, float)
    centered = expression - expression.mean(axis=1, keepdims=True)
    squares = np.square(centered).sum(axis=1)
    zero = squares == 0
    target_indices = np.asarray(target_indices, int)
    denominator = np.sqrt(squares[target_indices, None] * squares[None, :])
    signed = np.divide(centered[target_indices] @ centered.T, denominator,
                       out=np.zeros_like(denominator), where=denominator > 0)
    if np.max(np.abs(signed)) > 1 + 1e-12:
        raise ValueError("Correlation outside numerical tolerance")
    scores = np.clip(np.abs(signed), 0, 1)
    return scores, zero[target_indices], zero


def hash_check(path, expected):
    actual = sha256(path)
    if actual != expected: raise ValueError(f"Source hash mismatch: {path}")
    return actual


def frozen_evaluation(dataset, cfg, source_manifest):
    folder = SOURCE / dataset
    src = json.loads((folder / "source_manifest.json").read_text())
    response_dir = HERE / "results/library_calibration_v2" / f"{dataset}_disjoint_library_pools_seed20261024"
    response_manifest = response_dir / "manifest.json"
    hash_check(response_manifest, src["response_manifest_sha256"])
    hash_check(response_dir / "responses.npz", src["responses_sha256"])
    fit = folder / "fits/q1_thin0_model42"
    genes = (fit / "genes.txt").read_text().splitlines()
    targets = (fit / "targets.txt").read_text().splitlines()
    baseline = np.load(folder / "evaluations/baseline_q1_thin0_model42_S0_C0_raw_distance.npz", allow_pickle=False)
    common_genes = baseline["genes"].astype(str)
    if list(baseline["targets"].astype(str)) != targets:
        raise ValueError("Frozen baseline target order mismatch")
    common_indices = np.array([genes.index(g) for g in common_genes], int)
    target_indices = np.array([genes.index(t) for t in targets], int)
    z = np.load(response_dir / "responses.npz", allow_pickle=False)
    response_genes = list(z["genes"].astype(str)); response_tasks = list(z["tasks"].astype(str))
    all_direct = {t for t in response_tasks if not t.startswith("control_vs_control_")}
    if any(g in all_direct for g in common_genes):
        raise ValueError("A directly targeted gene remains in the frozen candidate scope")
    if list(np.asarray(genes)[[g not in all_direct for g in genes]]) != list(common_genes):
        raise ValueError("Candidate universe differs from frozen target exclusion")
    row = np.array([response_tasks.index(t) for t in targets], int)
    col = np.array([response_genes.index(g) for g in common_genes], int)
    null = np.array([i for i,t in enumerate(response_tasks) if t.startswith("control_vs_control_")], int)
    data = dict(dataset=dataset,genes=genes,targets=targets,common_genes=common_genes,
        common_indices=common_indices,target_indices=target_indices,C0=baseline["covariate"],
        B=z["B"][row][:,col],C=z["C"][row][:,col],
        null=z["B"][null][:,col]*z["C"][null][:,col],source=src,
        response_manifest_sha256=sha256(response_manifest),responses_sha256=sha256(response_dir/"responses.npz"))
    return data


def fit_association(data, q, seed, checks):
    fit = SOURCE/data["dataset"]/"fits"/f"q{q:g}_thin{seed}_model42"
    receipt = json.loads((fit/"input_manifest.json").read_text())
    if receipt["dataset"]!=data["dataset"] or receipt["q"]!=q or receipt["thinning_seed"]!=seed:
        raise ValueError("Frozen input condition mismatch")
    for name, expected in receipt["input_sha256"].items(): hash_check(fit/name, expected)
    if (fit/"genes.txt").read_text().splitlines()!=data["genes"] or (fit/"targets.txt").read_text().splitlines()!=data["targets"]:
        raise ValueError("Frozen gene or target order changed across conditions")
    counts = mmread(fit/"control_counts.mtx").toarray()
    if counts.shape!=(len(data["genes"]),400): raise ValueError("Fixed model input dimension changed")
    expression,totals = normalized_model_expression(counts)
    full_scores,zero_targets,zero_genes = absolute_pearson(expression,data["target_indices"])
    independent_error = None
    if q==1:
        with np.errstate(invalid="ignore",divide="ignore"):
            independent = np.abs(np.corrcoef(expression))[data["target_indices"]]
        independent = np.nan_to_num(independent,nan=0.0)
        independent_error = float(np.max(np.abs(independent-full_scores)))
        if not np.allclose(independent,full_scores,rtol=1e-12,atol=1e-12):
            raise ValueError("q1 correlation disagrees with independent numpy.corrcoef")
    checks.append(dict(dataset=data["dataset"],q=q,thinning_seed=seed,n_cells=counts.shape[1],model_genes=counts.shape[0],
        zero_library_cells=int((totals<=0).sum()),input_UMI_mass=int(counts.sum()),
        zero_variance_targets=int(zero_targets.sum()),zero_variance_model_genes=int(zero_genes.sum()),
        zero_variance_candidate_genes=int(zero_genes[data["common_indices"]].sum()),
        numpy_corrcoef_max_absolute_error=independent_error))
    return full_scores[:,data["common_indices"]],zero_targets,zero_genes[data["common_indices"]]


def transformed_scores(raw, cov):
    logs = np.log10(np.maximum(raw,1e-10))
    return {"raw_Pearson":raw,"EDR":np.vstack([residualize(row,cov) for row in logs])},logs


def dataset_run(dataset,cfg,source_manifest):
    data = frozen_evaluation(dataset,cfg,source_manifest)
    folder = OUT/dataset; folder.mkdir(parents=True,exist_ok=True)
    rows=[];null_rows=[];checks=[]
    s0,z0,v0 = fit_association(data,1.0,0,checks)
    c0 = data["C0"]
    baseline_methods,_ = transformed_scores(s0,c0)
    baseline_weights = {m:np.vstack([topk_weights(s,cfg["primary_k"]) for s in scores]) for m,scores in baseline_methods.items()}
    baseline_raw_evaluation = None

    def evaluate(raw,cov,q,seed,factor,kind,zero_targets,zero_genes):
        nonlocal baseline_raw_evaluation
        methods,logs = transformed_scores(raw,cov)
        for method,scores in methods.items():
            result,weights,utility,null = evaluate_scores(scores,cov,logs,data["B"],data["C"],data["null"],baseline_weights[method],data["targets"],cfg["primary_k"])
            result = result.rename(columns={"distance_floor_count":"score_floor_count"})
            result["target_zero_variance"] = zero_targets
            result["zero_variance_candidate_genes"] = int(zero_genes.sum())
            result["score_exact_constant"] = np.ptp(scores,axis=1)==0
            result["all_scores_floored"] = (raw<=1e-10).all(axis=1)
            context=dict(dataset=dataset,q=q,thinning_seed=seed,factor=factor,kind=kind,method=method)
            for key,value in context.items(): result[key]=value
            if kind=="baseline" and method=="raw_Pearson": baseline_raw_evaluation=(scores.copy(),weights.copy(),utility.copy(),null.copy())
            if factor in ["S0_C0","S0_Cq"] and method=="raw_Pearson":
                for actual,expected in zip((scores,weights,utility,null),baseline_raw_evaluation):
                    if not np.array_equal(actual,expected): raise AssertionError("Covariate-only condition changed raw score/evaluation")
            rows.append(result)
            for i in range(null.shape[0]):
                null_rows.append(dict(**context,null_partition=i,mean_over_target_rankings=float(null[i].mean()),
                    q95_over_target_rankings=float(np.quantile(null[i],.95))))
            p=folder/"evaluations"/f"{kind}_q{q:g}_thin{seed}_{factor}_{method}.npz"
            p.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(p,targets=np.asarray(data["targets"]),genes=data["common_genes"],scores=scores,
                raw_correlation=raw,covariate=cov,topk_weights=weights,utility_query_source=utility,null_utility=null,
                target_zero_variance=zero_targets,candidate_zero_variance=zero_genes)

    evaluate(s0,c0,1.0,0,"S0_C0","baseline",z0,v0)
    for q in cfg["thinning_probabilities"]:
        for seed in cfg["thinning_seeds"]:
            sq,zq,vq=fit_association(data,q,seed,checks)
            source_npz=SOURCE/dataset/"evaluations"/f"factorial_q{q:g}_thin{seed}_model42_S0_Cq_raw_distance.npz"
            observed=np.load(source_npz,allow_pickle=False)
            if not np.array_equal(observed["genes"],data["common_genes"]) or list(observed["targets"].astype(str))!=data["targets"]:
                raise ValueError("Cq frozen identity/order mismatch")
            cq=observed["covariate"]
            for factor,raw,c,z,v in [("S0_C0",s0,c0,z0,v0),("S0_Cq",s0,cq,z0,v0),
                                     ("Sq_C0",sq,c0,zq,vq),("Sq_Cq",sq,cq,zq,vq)]:
                evaluate(raw,c,q,seed,factor,"factorial",z,v)
            print(f"Completed association comparator {dataset} q={q} seed={seed}: four arms",flush=True)
    frame=pd.concat(rows,ignore_index=True)
    frame.to_csv(folder/"per_target_all_conditions.csv.gz",index=False)
    pd.DataFrame(null_rows).to_csv(folder/"control_null_partitions.csv.gz",index=False)
    pd.DataFrame(checks).to_csv(folder/"count_and_correlation_checks.csv",index=False)
    write_json(folder/"source_manifest.json",dict(dataset=dataset,source_measurement_manifest_sha256=sha256(SOURCE/"manifest.json"),
        response_manifest_sha256=data["response_manifest_sha256"],responses_sha256=data["responses_sha256"],
        fixed_model_cells=400,targets=len(data["targets"]),candidate_genes=len(data["common_genes"]),
        baseline_role="control-only association, not an intervention prediction model",
        normalization_scope="sum over every fixed model gene; historical native model preprocessing is not asserted to be identical"))
    return frame


def main():
    cfg=json.loads(CONFIG.read_text()); protocol_hash=sha256(CONFIG)
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/"manifest.json").exists(): raise SystemExit("Completed comparator exists; refusing overwrite")
    start=datetime.now(timezone.utc).isoformat()
    code_files=[Path(__file__),HERE/"measurement_diagnostic_pilot.py",HERE/"carf_core.py",HERE/"io_utils.py"]
    code_hash={p.name:sha256(p) for p in code_files}
    try:
        source_manifest=json.loads((SOURCE/"manifest.json").read_text())
        for rel,expected in source_manifest["output_sha256"].items(): hash_check(SOURCE/rel,expected)
        print(f"Verified all {len(source_manifest['output_sha256'])} registered measurement source files",flush=True)
        frame=pd.concat([dataset_run(ds,cfg,source_manifest) for ds in cfg["datasets"]],ignore_index=True)
        metrics=["top50_overlap","score_expression_rho","logscore_expression_r2","total_BC","identity_gain","null_response_mean","null_response_abs_mean"]
        keys=["dataset","kind","q","factor","method","target"]
        averaged=frame.groupby(keys,dropna=False)[metrics].mean().reset_index()
        averaged.to_csv(OUT/"seed_averaged_per_target.csv.gz",index=False)
        def summarize(data,columns,path):
            summary=data.groupby(keys[:-1],dropna=False)[columns].agg(["mean","median"])
            summary.columns=["_".join(c) for c in summary.columns]
            summary.reset_index().to_csv(path,index=False)
        summarize(averaged,metrics,OUT/"condition_summary.csv")
        base=frame[frame.kind=="baseline"][["dataset","method","target"]+metrics]
        paired=frame.merge(base,on=["dataset","method","target"],suffixes=("","_baseline"),validate="many_to_one")
        for m in metrics: paired[m+"_change"]=paired[m]-paired[m+"_baseline"]
        paired.to_csv(OUT/"paired_changes_all_conditions.csv.gz",index=False)
        changes=[m+"_change" for m in metrics]
        summarize(paired.groupby(keys,dropna=False)[changes].mean().reset_index(),changes,OUT/"paired_change_summary.csv")
        native=pd.read_csv(SOURCE/"seed_averaged_per_target.csv.gz")
        native=native[native.kind.isin(["baseline","factorial"])].copy()
        native["method"]=native.method.replace({"raw_distance":"raw_Pearson"})
        comparison=averaged.merge(native,on=keys,suffixes=("_association","_native"),validate="one_to_one")
        for m in metrics: comparison[m+"_association_minus_native"]=comparison[m+"_association"]-comparison[m+"_native"]
        comparison.to_csv(OUT/"native_comparison_per_target.csv.gz",index=False)
        comparison_columns=[m+suffix for m in metrics for suffix in ["_association","_native","_association_minus_native"]]
        comparison.groupby(keys[:-1],dropna=False)[comparison_columns].mean().reset_index().to_csv(OUT/"native_comparison_summary.csv",index=False)
        numerical=frame.groupby(["dataset","kind","q","thinning_seed","factor","method"]).agg(
            n_targets=("target","size"),zero_variance_targets=("target_zero_variance","sum"),
            zero_variance_candidate_genes=("zero_variance_candidate_genes","first"),all_floored_tasks=("all_scores_floored","sum"),
            exactly_constant_score_tasks=("score_exact_constant","sum"),mean_floor_gene_count=("score_floor_count","mean")).reset_index()
        numerical.to_csv(OUT/"numerical_information_audit.csv",index=False)
        for p in code_files:
            if code_hash[p.name]!=sha256(p): raise ValueError("Code changed during comparator run")
        if sha256(CONFIG)!=protocol_hash: raise ValueError("Protocol changed during comparator run")
        # The original source manifest and all its registered files are immutable.
        for rel,expected in source_manifest["output_sha256"].items(): hash_check(SOURCE/rel,expected)
        outputs={str(p.relative_to(OUT)):sha256(p) for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="manifest.json"}
        write_json(OUT/"manifest.json",dict(status="completed_development_association_comparator_not_model_validation",started_utc=start,
            finished_utc=datetime.now(timezone.utc).isoformat(),protocol=cfg,protocol_sha256=protocol_hash,
            code_sha256=code_hash,source_measurement_manifest_sha256=sha256(SOURCE/"manifest.json"),
            registered_source_files_checked_before_and_after=len(source_manifest["output_sha256"]),output_sha256=outputs,
            environment=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),
            comparison_column_note="raw_Pearson rows compare association raw correlation to native raw distance; EDR rows compare respective EDR scores on identical candidates and outcomes",
            limitation="Conditional-on-fixed-scope development baseline. This does not validate a second perturbation prediction model or establish causal effects."))
        print("Completed all prespecified association comparator conditions",flush=True)
    except Exception as exc:
        write_json(OUT/"failure.json",dict(time_utc=datetime.now(timezone.utc).isoformat(),error=repr(exc),protocol_sha256=protocol_hash))
        raise


if __name__=="__main__":main()
