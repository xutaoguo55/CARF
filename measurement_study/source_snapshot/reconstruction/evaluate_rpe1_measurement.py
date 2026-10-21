"""Frozen RPE1 factorial evaluation; outcomes never enter the model executor."""
from datetime import datetime, timezone
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from carf_core import residualize, topk_weights
from evaluate_heldout_linear import within_fold_utilities, row_correlation
from run_rpe1_measurement_models import rng_stream, thin_full_counts, normalize_full_counts
from io_utils import HERE, sha256, write_json

CONFIG=HERE/'config/rpe1_measurement_validation_v1.json'
BASE=HERE/'results/rpe1_measurement_validation_v1'
PREP=BASE/'preparation'
MODELS=BASE/'models_safe_export'
OUT=BASE/'evaluation'
MAIN=['native_scTenifoldKnk','control_only_Pearson','heldout_linear','training_condition_mean']

def verify_manifest(folder):
    p=folder/'manifest.json'; m=json.loads(p.read_text())
    for name,expected in m['output_sha256'].items():
        if sha256(folder/name)!=expected:raise ValueError(f'Input hash mismatch: {folder/name}')
    return m,sha256(p)

def scale_reference(A,B,C,null_product,folds,k=50):
    w,_,same,other=within_fold_utilities(np.abs(A),B,C,folds,k)
    null=(null_product@w.T/k).mean(axis=1)
    P=float(same.mean()); N=float(null.mean()); delta=P-N
    return dict(empirical_A_total=P,null_partition_mean=N,null_partition_max=float(null.max()),
                total_scale=delta,observed_separation=bool(delta>0 and P>float(null.max())),
                empirical_A_identity=float((same-other).mean()),null_partitions=null.tolist())

def basic_evaluation(raw,cov,B,C,null_product,folds,adjustment,k=50,signed=None):
    logs=np.log10(np.maximum(raw,1e-10))
    scores=raw if adjustment=='raw' else np.vstack([residualize(row,cov) for row in logs])
    w,u,same,other=within_fold_utilities(scores,B,C,folds,k)
    null=null_product@w.T/k
    sst=np.square(logs-logs.mean(axis=1,keepdims=True)).sum(axis=1)
    residual=np.vstack([residualize(row,cov) for row in logs])
    r2=np.divide(np.square(residual).sum(axis=1),sst,out=np.full(len(raw),np.nan),where=sst>1e-20)
    data=dict(total_BC=same,other_target_BC=other,identity_gain=same-other,
        null_response_mean=null.mean(axis=0),null_response_abs_mean=np.abs(null).mean(axis=0),
        score_floor_count=(raw<=1e-10).sum(axis=1),all_scores_floored=(raw<=1e-10).all(axis=1),
        score_exact_constant=np.ptp(scores,axis=1)==0,raw_logscore_expression_R2=1-r2)
    for arm,truth in [('B',B),('C',C)]:
        data['unadjusted_prediction_MSE_'+arm]=np.square(signed-truth).mean(axis=1) if signed is not None else np.full(len(raw),np.nan)
        data['unadjusted_prediction_Pearson_'+arm]=row_correlation(signed,truth) if signed is not None else np.full(len(raw),np.nan)
    return pd.DataFrame(data),scores,w,u,null

def covariate_contrast(left,right,left_weights,right_weights,scale,k=50):
    """Right minus left at fixed model measurement; ranking agreement is symmetric."""
    out=pd.DataFrame({'top50_agreement':np.minimum(left_weights,right_weights).sum(axis=1)/k})
    for name in ['total_BC','identity_gain']:
        delta=right[name].to_numpy()-left[name].to_numpy()
        out[name+'_change']=delta;out[name+'_absolute_change']=abs(delta)
    out['total_absolute_change_over_A_scale']=out.total_BC_absolute_change/scale['total_scale'] if scale['observed_separation'] else np.nan
    out['identity_absolute_change_over_A_identity']=out.identity_gain_absolute_change/scale['empirical_A_identity'] if scale['empirical_A_identity']>0 else np.nan
    out['uninformative_at_either_end']=(left.all_scores_floored|right.all_scores_floored|left.score_exact_constant|right.score_exact_constant).to_numpy()
    return out

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'manifest.json').exists():raise SystemExit('Completed evaluation exists; refusing overwrite')
    start=datetime.now(timezone.utc).isoformat();cfg=json.loads(CONFIG.read_text()); cfg_hash=sha256(CONFIG)
    pm,ph=verify_manifest(PREP); mm,mh=verify_manifest(MODELS)
    if pm['protocol_sha256']!=cfg_hash or mm['protocol_sha256']!=cfg_hash:raise ValueError('Protocol mismatch')
    if pm.get('status')!='completed_frozen_preparation_no_model_or_thinning' or mm.get('status')!='completed_safe_metadata_export_without_numeric_changes':raise ValueError('Incomplete preparation or safe score export')
    if mm.get('preparation_manifest_sha256')!=ph:raise ValueError('Models do not link to this exact preparation')
    original_model_manifest=BASE/'models/manifest.json'
    if mm.get('source_model_manifest_sha256')!=sha256(original_model_manifest):raise ValueError('Safe export does not link to original model execution')
    code_files=[Path(__file__),HERE/'carf_core.py',HERE/'evaluate_heldout_linear.py',HERE/'run_rpe1_measurement_models.py',HERE/'io_utils.py']
    code_hash={p.name:sha256(p) for p in code_files}
    z=np.load(PREP/'evaluation_responses.npz',allow_pickle=False)
    targets=z['targets'].astype(str);genes=z['genes'].astype(str);folds=z['folds'].astype(int)
    A,B,C=[np.asarray(z[k],float) for k in ['A','B','C']]
    null_product=np.asarray(z['null_B'])*np.asarray(z['null_C'])
    if not all(x.shape==(len(targets),len(genes)) for x in [A,B,C]):raise ValueError('Response dimension mismatch')
    if len(set(targets))!=len(targets) or len(set(genes))!=len(genes) or null_product.shape!=(30,len(genes)):raise ValueError('Response identity/null mismatch')
    if not all(np.isfinite(x).all() for x in [A,B,C,null_product]):raise ValueError('Nonfinite reference')
    ei=pd.read_csv(PREP/'evaluation_genes.csv')
    if not np.array_equal(ei.gene_id.astype(str),genes):raise ValueError('Candidate order mismatch')
    source_col=ei.source_gene_index.to_numpy(int)
    dc=load_npz(PREP/'Dc_counts.npz').toarray()
    dc_cells=pd.read_csv(PREP/'Dc_cells.csv')
    small_idx=np.load(PREP/'Dc64_indices.npy',allow_pickle=False).astype(int)
    if len(dc_cells)!=400 or dc_cells.cell_index.duplicated().any() or np.any(small_idx<0) or np.any(small_idx>=400):raise ValueError('Covariate cell identity/index mismatch')
    if dc.shape!=(400,8749) or small_idx.shape!=(64,) or len(set(small_idx))!=64:raise ValueError('Frozen covariate pool mismatch')
    covs={};cov_audit=[]
    for q,seed in [(1.,0)]+[(q,s) for q in cfg['thinning_probabilities'] for s in cfg['thinning_seeds']]:
        description=f'Dc|q={q:g}'
        thinned=thin_full_counts(dc,q,rng_stream(cfg['protocol_id'],description,seed))
        normalized,totals=normalize_full_counts(thinned)
        for n,idx in [(400,np.arange(400)),(64,small_idx)]:
            cov=normalized[idx].mean(axis=0)[source_col]
            covs[(q,seed,n)]=cov
            cov_audit.append(dict(q=q,thinning_seed=seed,n=n,rng_description=description,
                retained_UMI=int(thinned[idx].sum()),min_observed_gene_library=float(totals[idx].min())))
    np.savez_compressed(OUT/'covariates.npz',genes=genes,**{f'q{q:g}_seed{s}_n{n}':v for (q,s,n),v in covs.items()})
    pd.DataFrame(cov_audit).to_csv(OUT/'covariate_measurement_audit.csv',index=False)
    scale=scale_reference(A,B,C,null_product,folds,cfg['primary_k']);write_json(OUT/'fixed_empirical_A_scale.json',scale)
    index=pd.read_csv(MODELS/'model_index.csv');records={}
    for rec in index.to_dict('records'):
        p=MODELS/rec['path']
        if sha256(p)!=rec['sha256']:raise ValueError('Model index hash mismatch')
        f=np.load(p,allow_pickle=False)
        if not np.array_equal(f['targets'].astype(str),targets) or not np.array_equal(f['genes'].astype(str),genes) or not np.array_equal(f['folds'],folds):raise ValueError('Prediction identity/order mismatch')
        records[(rec['model'],float(rec['q']),int(rec['thinning_seed']),int(rec['algorithm_seed']))]=dict(rec=rec,raw=f['raw_scores'],signed=f['signed_predictions'] if 'signed_predictions' in f.files else None)
    if len(records)!=35:raise ValueError('Expected all 35 frozen model conditions')
    for name,raw,signed in [('observed_A_technical_anchor',np.abs(A),A),('fixed_control_expression',np.tile(covs[(1.,0,400)],(len(targets),1)),None),('constant_tie',np.zeros_like(A),None)]:
        records[(name,1.,0,0)]=dict(rec={'path':'fixed_evaluation_anchor'},raw=raw,signed=signed)
    rows=[]; contrasts=[];null_rows=[];baselines={};raw_checks={};invariants=[];serial=0
    eval_dir=OUT/'evaluations';eval_dir.mkdir(exist_ok=True)

    def run(key,q,seed,n,factor,kind,cov_key):
        nonlocal serial
        model,mq,ms,algorithm_seed=key;record=records[key]; raw=record['raw']; cov=covs[cov_key]
        both={}
        for adjustment in ['raw','EDR']:
            f,scores,w,u,null=basic_evaluation(raw,cov,B,C,null_product,folds,adjustment,cfg['primary_k'],record['signed'])
            basekey=(model,adjustment,n)
            is_primary_seed=(algorithm_seed in [0,42])
            if kind=='baseline' and is_primary_seed:baselines[basekey]=w.copy()
            reference=baselines[basekey]
            f['top50_agreement_to_q1_same_n']=np.minimum(w,reference).sum(1)/cfg['primary_k']
            f['top50_agreement_to_q1_n400']=np.minimum(w,baselines[(model,adjustment,400)]).sum(1)/cfg['primary_k']
            if adjustment=='raw':
                if key not in raw_checks:raw_checks[key]=(w.copy(),u.copy(),null.copy())
                else:
                    for actual,expected in zip((w,u,null),raw_checks[key]):
                        if not np.array_equal(actual,expected,equal_nan=True):raise AssertionError('C-only changed raw evaluation')
                    invariants.append(dict(model=model,model_q=mq,model_seed=ms,algorithm_seed=algorithm_seed,q=q,seed=seed,n=n,passed=True))
            if model in ['training_condition_mean','fixed_control_expression','constant_tie'] and np.max(np.abs(f.identity_gain))>1e-12:raise AssertionError('Same-fold constant ranking has nonzero identity')
            context=dict(model=model,q=q,thinning_seed=seed,n_covariate=n,algorithm_seed=algorithm_seed,
                factor=factor,kind=kind,adjustment=adjustment,model_q=mq,model_thinning_seed=ms,covariate_q=cov_key[0])
            f['target']=targets;f['fold']=folds
            for k,v in context.items():f[k]=v
            rows.append(f);both[adjustment]=(f,w)
            for j in range(len(null)):
                null_rows.append(dict(**context,null_partition=j,mean_utility_over_target_rankings=float(null[j].mean())))
            serial+=1
            np.savez_compressed(eval_dir/f'{serial:04d}.npz',targets=targets,genes=genes,folds=folds,scores=scores,
                                topk_weights=w,within_fold_utility=u,null_utility=null)
            f['evaluation_file']=f'evaluations/{serial:04d}.npz'
        return both

    # Primary q1 seeds first, so alternate native algorithm seeds have a fixed reference.
    baseline_keys=sorted([k for k in records if k[1]==1.],key=lambda k:(k[3] not in [0,42],k[0],k[3]))
    for key in baseline_keys:
        for n in [400,64]:run(key,1.,0,n,'S0_C0','baseline' if key[3] in [0,42] else 'algorithm_seed',(1.,0,n))
    for q in cfg['thinning_probabilities']:
        for seed in cfg['thinning_seeds']:
            for model in MAIN:
                algorithm_seed=42 if model=='native_scTenifoldKnk' else 0
                s0=(model,1.,0,algorithm_seed);sq=(model,float(q),seed,algorithm_seed)
                for n in [400,64]:
                    out={}
                    for factor,key,ck in [('S0_C0',s0,(1.,0,n)),('S0_Cq',s0,(q,seed,n)),('Sq_C0',sq,(1.,0,n)),('Sq_Cq',sq,(q,seed,n))]:
                        out[factor]=run(key,q,seed,n,factor,'factorial',ck)
                    for adjustment in ['raw','EDR']:
                        for effect,l,r in [('C_at_S0','S0_C0','S0_Cq'),('C_at_Sq','Sq_C0','Sq_Cq')]:
                            left,lw=out[l][adjustment];right,rw=out[r][adjustment]
                            c=covariate_contrast(left,right,lw,rw,scale,cfg['primary_k'])
                            c['target']=targets;c['fold']=folds
                            for k,v in dict(model=model,q=q,thinning_seed=seed,n_covariate=n,adjustment=adjustment,effect=effect).items():c[k]=v
                            contrasts.append(c)
            for model in ['observed_A_technical_anchor','fixed_control_expression','constant_tie']:
                for n in [400,64]:run((model,1.,0,0),q,seed,n,'S0_Cq','anchor',(q,seed,n))
    frame=pd.concat(rows,ignore_index=True);contrast=pd.concat(contrasts,ignore_index=True)
    frame.to_csv(OUT/'per_target_metrics.csv.gz',index=False);contrast.to_csv(OUT/'covariate_simple_effects.csv.gz',index=False)
    pd.DataFrame(null_rows).to_csv(OUT/'control_null_metrics.csv.gz',index=False)
    pd.DataFrame(invariants).to_csv(OUT/'raw_invariance_checks.csv',index=False)
    metrics=['total_BC','identity_gain','top50_agreement_to_q1_same_n','top50_agreement_to_q1_n400','null_response_mean','null_response_abs_mean',
             'unadjusted_prediction_MSE_B','unadjusted_prediction_MSE_C','unadjusted_prediction_Pearson_B','unadjusted_prediction_Pearson_C']
    keys=['model','q','n_covariate','algorithm_seed','factor','kind','adjustment']
    per_target=frame.groupby(keys+['target','fold'],dropna=False)[metrics].mean().reset_index()
    per_target.groupby(keys,dropna=False)[metrics].mean().reset_index().to_csv(OUT/'condition_summary.csv',index=False)
    per_target.groupby(keys+['fold'],dropna=False)[metrics].mean().reset_index().to_csv(OUT/'fold_summary.csv',index=False)
    frame.groupby(keys+['thinning_seed'],dropna=False)[metrics].mean().reset_index().to_csv(OUT/'seed_summary.csv',index=False)
    # The prediction-error fields refer to unadjusted signed predictions in both rows.
    # Only ranking/response metrics change under EDR.
    comparison_keys=['model','q','thinning_seed','n_covariate','algorithm_seed','factor','kind','target','fold']
    raw=frame[frame.adjustment=='raw']
    paired_adjustment=frame[frame.adjustment=='EDR'].merge(raw[comparison_keys+metrics],on=comparison_keys,
        suffixes=('','_raw'),validate='one_to_one')
    for metric in metrics:paired_adjustment[metric+'_EDR_minus_raw']=paired_adjustment[metric]-paired_adjustment[metric+'_raw']
    paired_adjustment.to_csv(OUT/'paired_adjustment_effects.csv.gz',index=False)
    mean_keys=['q','thinning_seed','n_covariate','factor','kind','adjustment','target','fold']
    comparable=frame[frame.kind.isin(['baseline','factorial'])]
    mean=comparable[comparable.model=='training_condition_mean']
    paired_mean=comparable.merge(mean[mean_keys+metrics],on=mean_keys,suffixes=('','_training_mean'),validate='many_to_one')
    for metric in metrics:paired_mean[metric+'_minus_training_mean']=paired_mean[metric]-paired_mean[metric+'_training_mean']
    paired_mean.to_csv(OUT/'paired_training_mean_effects.csv.gz',index=False)
    budget_keys=['model','q','thinning_seed','algorithm_seed','factor','kind','adjustment','target','fold']
    large=frame[frame.n_covariate==400]
    budget=frame[frame.n_covariate==64].merge(large[budget_keys+metrics],on=budget_keys,suffixes=('','_n400'),validate='one_to_one')
    for metric in metrics:budget[metric+'_n64_minus_n400']=budget[metric]-budget[metric+'_n400']
    budget.to_csv(OUT/'paired_covariate_budget_effects.csv.gz',index=False)
    practical=[];seed_practical=[];ckeys=['model','q','n_covariate','adjustment','effect']
    for key,g in contrast.groupby(ckeys):
        bytarget=g.groupby('target').mean(numeric_only=True)
        mean=float(bytarget.top50_agreement.mean());p10=float(bytarget.top50_agreement.quantile(.1))
        ratio=float(g.total_absolute_change_over_A_scale.mean()) if scale['observed_separation'] else None
        numerical=int(g.groupby('target').uninformative_at_either_end.max().sum())
        seed_means=g.groupby('thinning_seed').top50_agreement.mean()
        rank_pass=bool(mean>=.90 and p10>=.80)
        utility_pass=bool(ratio<=.1) if ratio is not None else None
        seed_rank_pass=[];seed_utility_pass=[]
        for seed,one in g.groupby('thinning_seed'):
            sm=float(one.top50_agreement.mean());sp=float(one.top50_agreement.quantile(.1))
            sr=float(one.total_absolute_change_over_A_scale.mean()) if scale['observed_separation'] else None
            sn=int(one.uninformative_at_either_end.sum());rp=bool(sm>=.90 and sp>=.80)
            up=bool(sr<=.1) if sr is not None else None
            seed_rank_pass.append(rp);seed_utility_pass.append(up)
            seed_practical.append(dict(zip(ckeys,key),thinning_seed=int(seed),agreement_mean=sm,
                agreement_target_p10=sp,mean_absolute_total_change_over_A_scale=sr,
                ranking_threshold_pass=rp,utility_small_change=up,numerically_uninformative_targets=sn,
                ranking_interpretation='information_insufficient_for_full_scope' if sn else ('threshold_met' if rp else 'threshold_not_met'),
                utility_interpretation='information_insufficient_for_scaling' if up is None else ('threshold_met' if up else 'threshold_not_met')))
        practical.append(dict(zip(ckeys,key),agreement_mean=mean,agreement_target_p10=p10,
            targets_below_point8=int((bytarget.top50_agreement<.8).sum()),numerically_uninformative_targets=numerical,
            ranking_threshold_pass=rank_pass,
            ranking_interpretation='information_insufficient_for_full_scope' if numerical else ('threshold_met' if rank_pass else 'threshold_not_met'),
            utility_interpretation='information_insufficient_for_scaling' if utility_pass is None else ('threshold_met' if utility_pass else 'threshold_not_met'),
            ranking_seed_threshold_crossing=len(set(seed_rank_pass))>1,
            utility_seed_threshold_crossing=len(set(seed_utility_pass))>1,
            mean_absolute_total_change=float(g.total_BC_absolute_change.mean()),
            mean_signed_total_change=float(g.total_BC_change.mean()),absolute_total_change_p90=float(g.total_BC_absolute_change.quantile(.9)),
            mean_absolute_total_change_over_A_scale=ratio,utility_small_change=utility_pass,
            seed_agreement_min=float(seed_means.min()),seed_agreement_max=float(seed_means.max()),
            interpretation='operational finite-benchmark criterion; not equivalence or biological inference'))
    pd.DataFrame(practical).to_csv(OUT/'operational_criteria.csv',index=False)
    pd.DataFrame(seed_practical).to_csv(OUT/'operational_criteria_by_seed.csv',index=False)
    # Exact four-condition decomposition, all evaluated on the same frozen outcomes.
    factorial=frame[frame.kind=='factorial'];decomposed=[]
    for metric in ['total_BC','identity_gain']:
        p=factorial.pivot(index=['model','q','thinning_seed','n_covariate','adjustment','target','fold'],columns='factor',values=metric)
        if p.isna().any().any():raise ValueError('Incomplete factorial condition')
        d=p.reset_index()[['model','q','thinning_seed','n_covariate','adjustment','target','fold']].copy()
        d['metric']=metric;d['covariate_at_S0']=(p.S0_Cq-p.S0_C0).to_numpy();d['model_at_C0']=(p.Sq_C0-p.S0_C0).to_numpy()
        d['interaction']=(p.Sq_Cq-p.Sq_C0-p.S0_Cq+p.S0_C0).to_numpy();d['joint']=(p.Sq_Cq-p.S0_C0).to_numpy()
        if not np.allclose(d.covariate_at_S0+d.model_at_C0+d.interaction,d.joint,atol=1e-14,rtol=0):raise AssertionError('Decomposition arithmetic')
        decomposed.append(d)
    pd.concat(decomposed).to_csv(OUT/'factorial_decomposition_per_target.csv.gz',index=False)
    if sha256(CONFIG)!=cfg_hash or sha256(PREP/'manifest.json')!=ph or sha256(MODELS/'manifest.json')!=mh:raise ValueError('Inputs changed during evaluation')
    for p in code_files:
        if sha256(p)!=code_hash[p.name]:raise ValueError('Code changed during evaluation')
    output={str(p.relative_to(OUT)):sha256(p) for p in OUT.rglob('*') if p.is_file()}
    write_json(OUT/'manifest.json',dict(status='completed',started_utc=start,finished_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=cfg_hash,preparation_manifest_sha256=ph,model_manifest_sha256=mh,code_sha256=code_hash,
        n_targets=len(targets),n_candidates=len(genes),target_records=len(frame),output_sha256=output,
        precision='Descriptive finite-benchmark results; seeds, targets and GEM groups are not independent biological replicates'))
    print(f'Completed {len(frame)} target records; {len(output)} registered outputs',flush=True)

if __name__=='__main__':main()
