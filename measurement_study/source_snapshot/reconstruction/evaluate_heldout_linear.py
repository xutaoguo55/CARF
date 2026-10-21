"""Published bilinear architecture adapted to strict target-held-out A training.

All source-gene readouts enter fold-local PCA. Only training perturbation A cells
enter model fitting. B/C are evaluation outcomes; empirical A is a labeled anchor.
The control column used in the original publication is deliberately omitted by
the frozen adaptation protocol. Scores are evaluated only within test folds.
"""
from __future__ import annotations
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import h5py
import numpy as np
import pandas as pd

from carf_core import residualize, topk_weights
from io_utils import HERE, column, dense_rows, gene_symbols, input_path, sha256, write_json

CONFIG=HERE/"config/heldout_linear_diagnostic_v1.json"
OUT=HERE/"results/heldout_linear_diagnostic_v1"
MEASUREMENT=HERE/"results/measurement_diagnostic_pilot_v1"
PEARSON=HERE/"results/pearson_measurement_comparator_v1"


def check_hash(path,expected):
    actual=sha256(path)
    if actual!=expected: raise ValueError(f"Source hash mismatch: {path}")
    return actual


def stable_folds(targets,dataset,n_folds=3):
    if len(set(targets))!=len(targets): raise ValueError("Duplicate targets")
    ranked=sorted(targets,key=lambda target:(hashlib.sha256(f"heldout_linear_diagnostic_v1|{dataset}|{target}".encode()).hexdigest(),target))
    assignment={target:i%n_folds for i,target in enumerate(ranked)}
    return np.array([assignment[t] for t in targets],int)


def bilinear_components(Y_train,training_target_gene_indices,k=10,ridge=.1):
    """R prcomp orientation: genes are observations, condition columns centered."""
    Y=np.asarray(Y_train,float)
    if Y.ndim!=2 or not np.isfinite(Y).all(): raise ValueError("Invalid training expression")
    if not 0<k<min(Y.shape): raise ValueError("Fixed PCA rank unsupported by training dimensions")
    center_columns=Y-Y.mean(axis=0,keepdims=True)
    U,s,Vt=np.linalg.svd(center_columns,full_matrices=False)
    G=U[:,:k]*s[:k]
    b=Y.mean(axis=1,keepdims=True)
    P=G[np.asarray(training_target_gene_indices,int)]
    W=bilinear_solve(Y-b,G,P,ridge)
    return G,b,P,W,s[:k]


def bilinear_solve(centered_Y,G,P,ridge=.1):
    left=G.T@G+ridge*np.eye(G.shape[1])
    right=P.T@P+ridge*np.eye(P.shape[1])
    rhs=G.T@centered_Y@P
    return np.linalg.solve(right.T,np.linalg.solve(left,rhs).T).T


def predict_selected_fold(all_post,training_indices,test_indices,target_gene_indices,k=10,ridge=.1):
    """No evaluation response argument; test A columns are never selected."""
    training_indices=np.asarray(training_indices,int);test_indices=np.asarray(test_indices,int)
    if np.intersect1d(training_indices,test_indices).size:raise ValueError("Training/test target overlap")
    G,b,P,W,s=bilinear_components(all_post[:,training_indices],np.asarray(target_gene_indices)[training_indices],k,ridge)
    pred=G@W@G[np.asarray(target_gene_indices)[test_indices]].T+b
    return pred,(G,b,P,W,s)


def within_fold_utilities(scores,B,C,folds,k):
    """Never compare a query with a ranking produced in another training fold."""
    weights=np.vstack([topk_weights(s,k) for s in scores])
    n=len(scores);utility=np.full((n,n),np.nan);same=np.empty(n);other=np.empty(n)
    for fold in np.unique(folds):
        idx=np.flatnonzero(folds==fold)
        if len(idx)<2:raise ValueError("Fewer than two held-out targets in fold")
        u=(B[idx]*C[idx])@weights[idx].T/k
        utility[np.ix_(idx,idx)]=u
        same[idx]=np.diag(u)
        other[idx]=(u.sum(axis=1)-np.diag(u))/(len(idx)-1)
    return weights,utility,same,other


def row_correlation(x,y):
    a=x-x.mean(axis=1,keepdims=True);b=y-y.mean(axis=1,keepdims=True)
    denominator=np.sqrt(np.square(a).sum(axis=1)*np.square(b).sum(axis=1))
    return np.divide((a*b).sum(axis=1),denominator,out=np.full(len(x),np.nan),where=denominator>0)


def source_dataset(dataset,targets,folds,cfg):
    print(f"Reconstructing A-only source means: {dataset}",flush=True)
    folder=OUT/dataset;folder.mkdir(parents=True,exist_ok=True)
    src=json.loads((MEASUREMENT/dataset/"source_manifest.json").read_text())
    run=HERE/"results/library_calibration_v2"/f"{dataset}_disjoint_library_pools_seed{cfg['source_response_seed']}"
    check_hash(run/"manifest.json",src["response_manifest_sha256"])
    receipt=json.loads((run/"manifest.json").read_text())
    for name in ["responses.npz","cell_splits.csv.gz"]:check_hash(run/name,receipt["outputs"][name])
    source=input_path(dataset);check_hash(source,src["source_sha256"])
    frozen=np.load(run/"responses.npz",allow_pickle=False)
    genes=frozen["genes"].astype(str);tasks=list(frozen["tasks"].astype(str))
    if len(set(genes))!=len(genes) or not set(targets)<=set(genes):raise ValueError("Missing/duplicate target or gene")
    target_gene=np.array([list(genes).index(t) for t in targets],int)
    target_row=np.array([tasks.index(t) for t in targets],int)
    splits=pd.read_csv(run/"cell_splits.csv.gz")
    dm=pd.read_csv(MEASUREMENT/dataset/"model_cells.csv")
    if len(dm)!=400 or dm.cell_index.duplicated().any():raise ValueError("Fixed Dm cells changed")
    if set(dm.cell_index)&set(splits.loc[splits.arm!="covariates_D","cell_index"]):raise ValueError("Dm overlaps A/B/C")
    A_indices=[splits.loc[(splits.task==t)&(splits.arm=="pert_A"),"cell_index"].to_numpy(int) for t in targets]
    all_A=np.concatenate(A_indices)
    if len(np.unique(all_A))!=len(all_A):raise ValueError("Pert_A cells overlap target tasks")
    forbidden=set(splits.loc[splits.arm.isin(["pert_B","pert_C","ctrl_B","ctrl_C"]),"cell_index"])
    if set(all_A)&forbidden:raise ValueError("A overlaps held-out B/C")
    cell_audit=[]
    for fold in np.unique(folds):
        train=np.flatnonzero(folds!=fold);test=np.flatnonzero(folds==fold)
        tc=np.concatenate([A_indices[i] for i in train]);vc=np.concatenate([A_indices[i] for i in test])
        if set(tc)&set(vc) or set(tc)&forbidden:raise ValueError("Fold training cell leakage")
        cell_audit.append(dict(fold=int(fold),training_targets=len(train),test_targets=len(test),
            training_A_cells=len(tc),heldout_A_cells=len(vc),overlap_training_test_A=0,overlap_training_BC=0))
    with h5py.File(source) as f:
        if not np.array_equal(gene_symbols(f).astype(str),genes):raise ValueError("Source gene order changed")
        labels=column(f["obs"],"perturbation")
        if not np.all(labels[dm.cell_index.to_numpy(int)]=="control"):raise ValueError("Noncontrol Dm cell")
        dm_mean=dense_rows(f,dm.cell_index.to_numpy(int)).mean(axis=0)
        post=[];reconstruction_error=[]
        for i,target in enumerate(targets):
            idx=A_indices[i]
            if len(idx)!=60 or not np.all(labels[idx]==target):raise ValueError("A cell identity/count changed")
            p=dense_rows(f,idx).mean(axis=0)
            ca=splits.loc[(splits.task==target)&(splits.arm=="ctrl_A"),"cell_index"].to_numpy(int)
            check=p-dense_rows(f,ca).mean(axis=0)
            err=float(np.max(np.abs(check-frozen["A"][target_row[i]])))
            if err>1e-12:raise ValueError("Reconstructed A delta differs from frozen response")
            post.append(p);reconstruction_error.append(err)
    post=np.column_stack(post)
    baseline=np.load(MEASUREMENT/dataset/"evaluations/baseline_q1_thin0_model42_S0_C0_raw_distance.npz",allow_pickle=False)
    common=baseline["genes"].astype(str)
    if list(baseline["targets"].astype(str))!=targets:raise ValueError("Native target order differs")
    columns=np.array([list(genes).index(g) for g in common],int)
    all_direct={t for t in tasks if not t.startswith("control_vs_control_")}
    if set(common)&all_direct or len(set(common))!=len(common):raise ValueError("Candidate exclusion/uniqueness failed")
    pearson=np.load(PEARSON/dataset/"evaluations/baseline_q1_thin0_S0_C0_raw_Pearson.npz",allow_pickle=False)
    if not np.array_equal(pearson["genes"],common) or list(pearson["targets"].astype(str))!=targets:
        raise ValueError("Pearson candidate/target mismatch")
    null_idx=np.array([i for i,t in enumerate(tasks) if t.startswith("control_vs_control_")],int)
    np.savez_compressed(folder/"A_post_and_Dm_reference.npz",genes=genes,targets=np.asarray(targets),A_post=post,
        Dm_mean=dm_mean,target_gene_indices=target_gene,candidate_gene_indices=columns)
    dm.to_csv(folder/"Dm_control_cells.csv",index=False)
    pd.DataFrame([dict(target=t,cell_index=int(c),fold=int(folds[i])) for i,t in enumerate(targets) for c in A_indices[i]]).to_csv(folder/"A_perturbed_cells.csv",index=False)
    write_json(folder/"source_audit.json",dict(source_path=str(source),source_sha256=src["source_sha256"],
        response_manifest_sha256=sha256(run/"manifest.json"),responses_sha256=sha256(run/"responses.npz"),
        splits_sha256=sha256(run/"cell_splits.csv.gz"),target_count=len(targets),PCA_readout_genes=len(genes),
        evaluation_candidates=len(common),model_Dm_cells=len(dm),max_A_reconstruction_error=max(reconstruction_error),
        folds=cell_audit,limitations=["PCA uses all source curated genes; upstream HVG selection used all cells.",
        "Absolute A expression and Dm baseline arise from different technical library pools, so library offsets can enter signed errors.",
        "Original publication includes a ctrl training column; this frozen adaptation deliberately omits it."]))
    return dict(dataset=dataset,targets=targets,folds=folds,genes=genes,target_gene=target_gene,post=post,Dm_mean=dm_mean,
        candidates=common,columns=columns,C0=baseline["covariate"],native_scores=baseline["scores"],pearson_scores=pearson["scores"],
        A=frozen["A"][target_row][:,columns],B=frozen["B"][target_row][:,columns],C=frozen["C"][target_row][:,columns],
        null=frozen["B"][null_idx][:,columns]*frozen["C"][null_idx][:,columns])


def fit_models(data,cfg):
    dataset=data["dataset"];n=len(data["targets"]);g=len(data["candidates"])
    names=["heldout_linear","training_condition_mean"]+[f"training_label_permutation_{s}" for s in cfg["training_label_permutation_seeds"]]
    predictions={name:np.full((n,g),np.nan) for name in names};permutation_records=[]
    for fold in np.unique(data["folds"]):
        train=np.flatnonzero(data["folds"]!=fold);test=np.flatnonzero(data["folds"]==fold)
        predicted,(G,b,P,W,s)=predict_selected_fold(data["post"],train,test,data["target_gene"],cfg["components"],cfg["ridge_penalty_each_side"])
        predictions["heldout_linear"][test]=(predicted-data["Dm_mean"][:,None])[data["columns"]].T
        predictions["training_condition_mean"][test]=np.tile((b[:,0]-data["Dm_mean"])[data["columns"]],(len(test),1))
        perm_W=[];orders=[]
        for seed in cfg["training_label_permutation_seeds"]:
            rng=np.random.default_rng(int.from_bytes(hashlib.sha256(f"heldout_linear_permutation|{dataset}|{fold}|{seed}".encode()).digest()[:16],"little"))
            order=rng.permutation(len(train));pp=P[order]
            wp=bilinear_solve(data["post"][:,train]-b,G,pp,cfg["ridge_penalty_each_side"])
            pred=G@wp@G[data["target_gene"][test]].T+b-data["Dm_mean"][:,None]
            predictions[f"training_label_permutation_{seed}"][test]=pred[data["columns"]].T
            perm_W.append(wp);orders.append(order)
            permutation_records.extend(dict(fold=int(fold),seed=seed,condition_target=data["targets"][train[j]],
                assigned_embedding_target=data["targets"][train[order[j]]]) for j in range(len(train)))
        folder=OUT/dataset/"fits";folder.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(folder/f"fold_{fold}.npz",G=G,b=b,P=P,W=W,singular_values=s,training_indices=train,test_indices=test,
            training_targets=np.asarray(data["targets"])[train],test_targets=np.asarray(data["targets"])[test],
            permutation_W=np.stack(perm_W),permutation_orders=np.stack(orders),seeds=np.asarray(cfg["training_label_permutation_seeds"]))
        print(f"Fitted {dataset} fold {fold}: {len(train)} training / {len(test)} held-out; all five label permutations",flush=True)
    if any(not np.isfinite(p).all() for p in predictions.values()):raise ValueError("Incomplete/nonfinite fold prediction")
    predictions["observed_A_technical_anchor"]=data["A"]
    pd.DataFrame(permutation_records).to_csv(OUT/dataset/"training_label_permutations.csv",index=False)
    np.savez_compressed(OUT/dataset/"signed_predictions.npz",targets=np.asarray(data["targets"]),genes=data["candidates"],folds=data["folds"],**predictions)
    return predictions


def evaluate_dataset(data,predictions,cfg):
    models={name:np.abs(pred) for name,pred in predictions.items()}
    models["native_scTenifoldKnk"]=data["native_scores"]
    models["control_only_Pearson"]=data["pearson_scores"]
    k=cfg["primary_k"];rows=[];null_rows=[];baseline={};raw_invariants={}
    conditions=[(1.0,0,"C0",data["C0"])]
    for q in cfg["thinning_probabilities"]:
        for seed in cfg["thinning_seeds"]:
            p=MEASUREMENT/data["dataset"]/"evaluations"/f"factorial_q{q:g}_thin{seed}_model42_S0_Cq_raw_distance.npz"
            z=np.load(p,allow_pickle=False)
            if not np.array_equal(z["genes"],data["candidates"]) or list(z["targets"].astype(str))!=data["targets"]:raise ValueError("Cq order changed")
            conditions.append((q,seed,"Cq",z["covariate"]))
    for q,seed,condition,cov in conditions:
        for name,raw in models.items():
            logscore=np.log10(np.maximum(raw,1e-10))
            methods={"raw":raw,"EDR":np.vstack([residualize(x,cov) for x in logscore])}
            for adjustment,scores in methods.items():
                weights,utility,same,other=within_fold_utilities(scores,data["B"],data["C"],data["folds"],k)
                null=data["null"]@weights.T/k
                if condition=="C0":baseline[(name,adjustment)]=weights.copy()
                overlap=np.minimum(weights,baseline[(name,adjustment)]).sum(axis=1)/k
                if adjustment=="raw":
                    if condition=="C0":raw_invariants[name]=(scores.copy(),weights.copy(),same.copy(),other.copy(),null.copy())
                    else:
                        for actual,expected in zip((scores,weights,same,other,null),raw_invariants[name]):
                            if not np.array_equal(actual,expected):raise AssertionError("C-only condition changed raw prediction/ranking/evaluation")
                if name=="training_condition_mean" and np.max(np.abs(same-other))>1e-12:
                    raise AssertionError("Same-fold invariant mean predictor has nonzero identity gain")
                pred=predictions.get(name)
                mse_B=np.square(pred-data["B"]).mean(axis=1) if pred is not None else np.full(len(scores),np.nan)
                mse_C=np.square(pred-data["C"]).mean(axis=1) if pred is not None else np.full(len(scores),np.nan)
                corr_B=row_correlation(pred,data["B"]) if pred is not None else np.full(len(scores),np.nan)
                corr_C=row_correlation(pred,data["C"]) if pred is not None else np.full(len(scores),np.nan)
                context=dict(dataset=data["dataset"],model=name,adjustment=adjustment,q=q,thinning_seed=seed,covariate_condition=condition)
                for i,target in enumerate(data["targets"]):
                    rows.append(dict(**context,target=target,fold=int(data["folds"][i]),n_same_fold_sources=int((data["folds"]==data["folds"][i]).sum()),
                        total_BC=float(same[i]),other_target_BC=float(other[i]),identity_gain=float(same[i]-other[i]),top50_overlap=float(overlap[i]),
                        null_response_mean=float(null[:,i].mean()),null_response_abs_mean=float(np.abs(null[:,i]).mean()),
                        signed_prediction_MSE_B=float(mse_B[i]),signed_prediction_MSE_C=float(mse_C[i]),
                        signed_prediction_Pearson_B=float(corr_B[i]),signed_prediction_Pearson_C=float(corr_C[i]),
                        signed_prediction_available=pred is not None,score_floor_count=int((raw[i]<=1e-10).sum()),
                        score_exact_constant=bool(np.ptp(scores[i])==0),all_scores_floored=bool(np.all(raw[i]<=1e-10))))
                for ni in range(null.shape[0]):
                    for fold in np.unique(data["folds"]):
                        idx=np.flatnonzero(data["folds"]==fold)
                        null_rows.append(dict(**context,fold=int(fold),null_partition=ni,mean_over_test_rankings=float(null[ni,idx].mean()),
                            q95_over_test_rankings=float(np.quantile(null[ni,idx],.95))))
                folder=OUT/data["dataset"]/"evaluations";folder.mkdir(parents=True,exist_ok=True)
                np.savez_compressed(folder/f"{name}_{adjustment}_q{q:g}_seed{seed}_{condition}.npz",targets=np.asarray(data["targets"]),
                    genes=data["candidates"],folds=data["folds"],scores=scores,covariate=cov,topk_weights=weights,
                    within_fold_utility=utility,null_utility=null)
    frame=pd.DataFrame(rows);frame.to_csv(OUT/data["dataset"]/"per_target_metrics.csv.gz",index=False)
    pd.DataFrame(null_rows).to_csv(OUT/data["dataset"]/"control_null_partitions.csv.gz",index=False)
    return frame


def main():
    cfg=json.loads(CONFIG.read_text());cfg_hash=sha256(CONFIG);OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/"manifest.json").exists():raise SystemExit("Completed held-out linear analysis exists; refusing overwrite")
    start=datetime.now(timezone.utc).isoformat();files=[Path(__file__),HERE/"carf_core.py",HERE/"io_utils.py"]
    code_hash={p.name:sha256(p) for p in files};source_manifests={}
    try:
        for source in [MEASUREMENT,PEARSON]:
            m=json.loads((source/"manifest.json").read_text())
            for rel,expected in m["output_sha256"].items():check_hash(source/rel,expected)
            source_manifests[str(source.relative_to(HERE))]=sha256(source/"manifest.json")
        design_audit=HERE/"provenance/linear_baseline_design_audit_v1.json"
        if not design_audit.exists():raise FileNotFoundError("Design audit must precede fitting")
        assignments={};fold_rows=[]
        for dataset in cfg["datasets"]:
            targets=(MEASUREMENT/dataset/"fits/q1_thin0_model42/targets.txt").read_text().splitlines()
            folds=stable_folds(targets,dataset,cfg["folds"]);assignments[dataset]=(targets,folds)
            fold_rows.extend(dict(dataset=dataset,target=t,fold=int(folds[i]),sha256_sort_key=hashlib.sha256(f"heldout_linear_diagnostic_v1|{dataset}|{t}".encode()).hexdigest()) for i,t in enumerate(targets))
        pd.DataFrame(fold_rows).to_csv(OUT/"target_folds.csv",index=False)
        fold_hash=sha256(OUT/"target_folds.csv")
        write_json(OUT/"fold_freeze.json",dict(created_before_any_fit_utc=datetime.now(timezone.utc).isoformat(),fold_file_sha256=fold_hash,protocol_sha256=cfg_hash))
        frames=[]
        for dataset in cfg["datasets"]:
            data=source_dataset(dataset,*assignments[dataset],cfg)
            predictions=fit_models(data,cfg)
            frames.append(evaluate_dataset(data,predictions,cfg))
        frame=pd.concat(frames,ignore_index=True)
        metrics=["total_BC","identity_gain","top50_overlap","null_response_mean","null_response_abs_mean", "signed_prediction_MSE_B","signed_prediction_MSE_C","signed_prediction_Pearson_B","signed_prediction_Pearson_C"]
        keys=["dataset","model","adjustment","q","covariate_condition","target","fold"]
        average=frame.groupby(keys,dropna=False)[metrics].mean().reset_index()
        average.to_csv(OUT/"seed_averaged_per_target.csv.gz",index=False)
        summary=average.groupby(keys[:-2],dropna=False)[metrics].agg(["mean","median"])
        summary.columns=["_".join(c) for c in summary.columns]
        summary.reset_index().to_csv(OUT/"condition_summary.csv",index=False)
        average.groupby(keys[:-2]+["fold"],dropna=False)[metrics].mean().reset_index().to_csv(OUT/"fold_summary.csv",index=False)
        mean=frame[frame.model=="training_condition_mean"]
        match=["dataset","adjustment","q","thinning_seed","covariate_condition","target","fold"]
        paired=frame.merge(mean[match+metrics],on=match,suffixes=("","_training_mean"),validate="many_to_one")
        differences=[m+"_minus_training_mean" for m in metrics]
        for m,diff in zip(metrics,differences):paired[diff]=paired[m]-paired[m+"_training_mean"]
        paired.to_csv(OUT/"paired_training_mean_differences.csv.gz",index=False)
        paired.groupby(keys,dropna=False)[differences].mean().reset_index().groupby(keys[:-2],dropna=False)[differences].mean().reset_index().to_csv(OUT/"paired_training_mean_summary.csv",index=False)
        frame.groupby(["dataset","model","adjustment","q","thinning_seed","covariate_condition"]).agg(n_targets=("target","size"),
            exactly_constant_score_tasks=("score_exact_constant","sum"),all_floored_tasks=("all_scores_floored","sum"),mean_gene_floor_count=("score_floor_count","mean")).reset_index().to_csv(OUT/"numerical_information_audit.csv",index=False)
        if sha256(CONFIG)!=cfg_hash or sha256(OUT/"target_folds.csv")!=fold_hash:raise ValueError("Frozen protocol/folds changed")
        for p in files:
            if sha256(p)!=code_hash[p.name]:raise ValueError("Code changed during fitting")
        for source in [MEASUREMENT,PEARSON]:
            m=json.loads((source/"manifest.json").read_text())
            if sha256(source/"manifest.json")!=source_manifests[str(source.relative_to(HERE))]:raise ValueError("Source manifest changed")
            for rel,expected in m["output_sha256"].items():check_hash(source/rel,expected)
        outputs={str(p.relative_to(OUT)):sha256(p) for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="manifest.json"}
        write_json(OUT/"manifest.json",dict(status="completed_development_target_heldout_architecture_adaptation",started_utc=start,
            finished_utc=datetime.now(timezone.utc).isoformat(),protocol=cfg,protocol_sha256=cfg_hash,code_sha256=code_hash,
            source_manifest_sha256=source_manifests,design_audit_sha256=sha256(design_audit),output_sha256=outputs,
            prediction_information_budgets={"heldout_linear_and_permutations":"other training-fold perturbed A cells", "training_condition_mean":"same training-fold perturbed A cells", "native_and_Pearson":"control-only Dm", "observed_A_technical_anchor":"query-specific observed A; not a predictive competitor"},
            signed_metric_note="MSE and Pearson always evaluate the underlying signed predicted delta; EDR is a ranking diagnostic and has no separate signed expression prediction. Native/Pearson unsigned association/distance scores receive missing signed-error metrics, never fabricated deltas.",
            method_adaptation="Full source-gene PCA; no ctrl training column; independent Dm subtraction; technical three-fold target splitting and fixed candidate evaluation, not reproduction of published benchmark numbers.",
            permutation_rng="SHA256(heldout_linear_permutation|dataset|fold|seed), first16bytes little-endian, numpy default_rng",
            environment=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),
            uncertainty="All summaries descriptive; technical splits share study background and controls. No biological-replicate P values or CIs."))
        print("Completed strict target-held-out linear diagnostic; all prespecified models and covariate conditions retained",flush=True)
    except Exception as exc:
        write_json(OUT/"failure.json",dict(time_utc=datetime.now(timezone.utc).isoformat(),error=repr(exc),protocol_sha256=cfg_hash,code_sha256=code_hash))
        raise


if __name__=="__main__":main()
