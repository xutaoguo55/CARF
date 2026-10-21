"""Post-hoc fixed diagnostic comparison. Never refit, edit or select source results."""
from datetime import datetime, timezone
from pathlib import Path
import json
import platform
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from carf_core import residualize, topk_weights
from io_utils import HERE, sha256, write_json

CONFIG = HERE / 'config/posthoc_diagnostic_sensitivity_v1.json'
FREEZE = HERE / 'provenance/posthoc_diagnostic_sensitivity_v1_freeze.json'
OUT = HERE / 'results/posthoc_diagnostic_sensitivity_v1'
RANKINGS = ['raw', 'EDR', 'rank_OLS']
CONTEXT = ['model','q','thinning_seed','n_covariate','algorithm_seed','factor','kind',
           'model_q','model_thinning_seed','covariate_q']
GROUP = ['model','q','n_covariate','algorithm_seed','factor','kind','ranking']
METRICS = ['total_BC','other_target_BC','identity_gain','null_response_mean','null_response_abs_mean',
           'matched_background_BC','matched_total_BC','matched_other_target_BC','matched_identity_gain',
           'matched_null_response_mean','matched_null_response_abs_mean',
           'top50_agreement_to_q1_same_n','top50_agreement_to_q1_n400',
           'expression_mean_difference','expression_mean_absolute_difference',
           'expression_rank_mean_difference','expression_rank_mean_absolute_difference']


def percentile_midrank(x):
    x = np.asarray(x, float)
    if x.ndim != 1 or len(x) == 0 or not np.isfinite(x).all():
        raise ValueError('Expected finite one-dimensional measurements')
    return (rankdata(x, method='average') - .5) / len(x)


def ranking_scores(raw, cov, ranking):
    raw = np.asarray(raw, float)
    if raw.ndim != 2 or not np.isfinite(raw).all() or np.any(raw < 0):
        raise ValueError('Expected finite nonnegative raw scores')
    if len(cov) != raw.shape[1]:
        raise ValueError('Covariate dimension mismatch')
    if ranking == 'raw':
        return raw.copy()
    if ranking == 'EDR':
        return np.vstack([residualize(np.log10(np.maximum(row, 1e-10)), cov) for row in raw])
    if ranking == 'rank_OLS':
        x = percentile_midrank(cov)
        return np.vstack([residualize(percentile_midrank(row), x) for row in raw])
    raise ValueError('Unknown ranking')


def expression_strata(cov, n_strata=5):
    ranks = percentile_midrank(cov)
    return np.minimum(n_strata-1, np.floor(n_strata*ranks).astype(int)), ranks


def match_background(weights, cov, k=50, n_strata=5, capacity_tolerance=1e-12):
    """Deterministic fractional background; unsupported rows remain wholly NA."""
    w = np.asarray(weights, float)
    if w.ndim != 2 or w.shape[1] != len(cov) or not np.isfinite(w).all():
        raise ValueError('Invalid matching dimensions/weights')
    if np.any(w < 0) or np.any(w > 1) or not np.allclose(w.sum(1), k, atol=1e-10, rtol=0):
        raise ValueError('Weights must be inclusion probabilities of total k')
    bins, xr = expression_strata(cov, n_strata)
    v = np.zeros_like(w); supported = np.ones(len(w), bool)
    W = np.zeros((len(w), n_strata)); H = W.copy(); count = np.zeros(n_strata, int)
    for b in range(n_strata):
        idx = bins == b; count[b] = int(idx.sum())
        W[:,b] = w[:,idx].sum(1); H[:,b] = (1-w[:,idx]).sum(1)
        use = W[:,b] > 0
        fail = use & ((H[:,b] <= 0) | (W[:,b]-H[:,b] > capacity_tolerance))
        supported &= ~fail
        ok = use & ~fail
        v[np.ix_(ok,idx)] = (W[ok,b]/H[ok,b])[:,None]*(1-w[np.ix_(ok,idx)])
    v[~supported] = np.nan
    if supported.any():
        if not np.allclose(v[supported].sum(1), k, atol=1e-10, rtol=0):
            raise AssertionError('Matched mass differs from k')
        if np.any(v[supported] < -1e-12) or np.any(w[supported]+v[supported] > 1+1e-12):
            raise AssertionError('Matching exceeds available inclusion capacity')
    mass_error = np.full(len(w), np.nan)
    total_error = np.full(len(w), np.nan)
    expression_delta = np.full(len(w), np.nan)
    rank_delta = np.full(len(w), np.nan)
    if supported.any():
        errors = np.column_stack([v[supported][:,bins==b].sum(1)-W[supported,b] for b in range(n_strata)])
        mass_error[supported] = np.abs(errors).max(1)
        total_error[supported] = v[supported].sum(1)-k
        expression_delta[supported] = (w[supported]-v[supported])@np.asarray(cov)/k
        rank_delta[supported] = (w[supported]-v[supported])@xr/k
    occupied_slack = np.where(W > 0, H-W, np.inf).min(1)
    info = dict(matching_supported=supported, matching_failed_strata=((W>0)&((H<=0)|(W-H>capacity_tolerance))).sum(1),
                matching_min_capacity_slack=occupied_slack, matching_max_capacity_deficit=np.maximum(W-H,0).max(1),
                matching_stratum_mass_max_error=mass_error, matching_total_mass_error=total_error,
                expression_mean_difference=expression_delta, expression_mean_absolute_difference=np.abs(expression_delta),
                expression_rank_mean_difference=rank_delta, expression_rank_mean_absolute_difference=np.abs(rank_delta))
    return v, supported, info, W, H, bins, count


def evaluate_weights(w, v, supported, B, C, null_product, folds, k=50):
    """Source pool is fixed: one unsupported source disables full-fold matched identity."""
    product = np.asarray(B)*np.asarray(C)
    n = len(w)
    if w.shape != product.shape or v.shape != w.shape or np.asarray(folds).shape != (n,):
        raise ValueError('Evaluation identity/dimension mismatch')
    same=np.empty(n);other=np.empty(n);mi=np.full(n,np.nan);mo=mi.copy()
    fold_supported=np.zeros(n,bool)
    own_match = np.full(n,np.nan);background = np.full(n,np.nan)
    own_match[supported] = ((w[supported]-v[supported])*product[supported]).sum(1)/k
    background[supported] = (v[supported]*product[supported]).sum(1)/k
    for fold in np.unique(folds):
        idx = np.flatnonzero(folds == fold)
        if len(idx)<2:raise ValueError('At least two same-fold sources required')
        u=product[idx]@w[idx].T/k
        same[idx]=np.diag(u);other[idx]=(u.sum(1)-np.diag(u))/(len(idx)-1)
        if supported[idx].all():
            mu=product[idx]@(w[idx]-v[idx]).T/k
            mo[idx]=(mu.sum(1)-np.diag(mu))/(len(idx)-1)
            mi[idx]=np.diag(mu)-mo[idx];fold_supported[idx]=True
    null=null_product@w.T/k
    matched_null=np.full((len(null_product),n),np.nan)
    matched_null[:,supported]=null_product@(w[supported]-v[supported]).T/k
    return dict(total_BC=same,other_target_BC=other,identity_gain=same-other,
                null_response_mean=null.mean(0),null_response_abs_mean=np.abs(null).mean(0),
                matched_background_BC=background,matched_total_BC=own_match,matched_other_target_BC=mo,
                matched_identity_gain=mi,matching_fold_supported=fold_supported,
                matched_null_response_mean=matched_null.mean(0),matched_null_response_abs_mean=np.abs(matched_null).mean(0))


def strict_mean(series):
    x=np.asarray(series,float)
    return float(x.mean()) if np.isfinite(x).all() else np.nan


def summarize(frame, keys, metrics):
    # Deliberately retain NA, unlike pandas' default skip-NA reductions.
    values=frame.groupby(keys,dropna=False,sort=True)[metrics].agg(strict_mean).reset_index()
    counts=frame.groupby(keys,dropna=False,sort=True).size().rename('records').reset_index()
    values=values.merge(counts,on=keys,validate='one_to_one')
    for metric in ['matched_total_BC','matched_identity_gain']:
        finite=frame.groupby(keys,dropna=False,sort=True)[metric].count().rename(metric+'_finite_count').reset_index()
        values=values.merge(finite,on=keys,validate='one_to_one')
    return values


def verify_used(folder, manifest, names):
    for name in names:
        if sha256(folder/name)!=manifest['output_sha256'][name]:
            raise ValueError(f'Input hash mismatch: {folder/name}')


def main():
    if OUT.exists():raise SystemExit('Output already exists; refuse overwrite or partial-result reuse')
    cfg=json.loads(CONFIG.read_text());receipt=json.loads(FREEZE.read_text())
    config_hash=sha256(CONFIG)
    if receipt['config_sha256']!=config_hash:raise ValueError('Configuration differs from frozen receipt')
    design=HERE.parent/receipt['design_path']
    if sha256(design)!=receipt['design_sha256']:raise ValueError('Design differs from freeze')
    start=datetime.now(timezone.utc).isoformat();base=HERE/cfg['source_base']
    prep=base/'preparation';models=base/'models_safe_export';evaluation=base/'evaluation'
    inputs={}
    manifests={}
    for folder,key in [(prep,'source_preparation_manifest_sha256'),(models,'source_safe_model_manifest_sha256'),(evaluation,'source_evaluation_manifest_sha256')]:
        actual=sha256(folder/'manifest.json')
        if actual!=cfg[key]:raise ValueError(f'Frozen source mismatch: {folder}')
        inputs[str((folder/'manifest.json').relative_to(HERE))]=actual
        manifests[folder.name]=json.loads((folder/'manifest.json').read_text())
    if sha256(HERE/'config/rpe1_measurement_validation_v1.json')!=cfg['source_protocol_sha256']:raise ValueError('Original protocol changed')
    used_p=['evaluation_responses.npz','evaluation_genes.csv']
    used_e=['per_target_metrics.csv.gz','covariates.npz','fixed_empirical_A_scale.json']
    verify_used(prep,manifests['preparation'],used_p);verify_used(evaluation,manifests['evaluation'],used_e)
    verify_used(models,manifests['models_safe_export'],list(manifests['models_safe_export']['output_sha256']))
    for folder,names in [(prep,used_p),(evaluation,used_e)]:
        for name in names:inputs[str((folder/name).relative_to(HERE))]=sha256(folder/name)
    original=pd.read_csv(evaluation/'per_target_metrics.csv.gz')
    z=np.load(prep/'evaluation_responses.npz',allow_pickle=False)
    targets=z['targets'].astype(str);genes=z['genes'].astype(str);folds=z['folds'].astype(int)
    A,B,C=[np.asarray(z[name],float) for name in ['A','B','C']]
    null=np.asarray(z['null_B'],float)*np.asarray(z['null_C'],float)
    if len(targets)!=cfg['expected_targets'] or len(genes)!=cfg['expected_genes']:raise ValueError('Unexpected frozen scope')
    if len(set(targets))!=len(targets) or len(set(genes))!=len(genes) or set(folds)!={0,1,2}:raise ValueError('Invalid identities/folds')
    if any(x.shape!=(len(targets),len(genes)) or not np.isfinite(x).all() for x in [A,B,C]):raise ValueError('Invalid response data')
    if null.shape!=(30,len(genes)) or not np.isfinite(null).all():raise ValueError('Invalid null products')
    if not np.array_equal(pd.read_csv(prep/'evaluation_genes.csv').gene_id.astype(str),genes):raise ValueError('Candidate order mismatch')
    cov_npz=np.load(evaluation/'covariates.npz',allow_pickle=False)
    if not np.array_equal(cov_npz['genes'].astype(str),genes):raise ValueError('Covariate gene mismatch')
    covs={name:cov_npz[name] for name in cov_npz.files if name!='genes'}
    records={};index=pd.read_csv(models/'model_index.csv')
    for rec in index.to_dict('records'):
        p=models/rec['path']
        if sha256(p)!=rec['sha256']:raise ValueError('Model index mismatch')
        source=np.load(p,allow_pickle=False)
        if not np.array_equal(source['targets'].astype(str),targets) or not np.array_equal(source['genes'].astype(str),genes) or not np.array_equal(source['folds'],folds):raise ValueError('Model identity mismatch')
        key=(rec['model'],float(rec['q']),int(rec['thinning_seed']),int(rec['algorithm_seed']))
        if key in records:raise ValueError('Duplicate model export')
        records[key]=np.asarray(source['raw_scores'],float)
    if len(records)!=cfg['expected_model_exports']:raise ValueError('Missing model exports')
    records[('observed_A_technical_anchor',1.,0,0)]=np.abs(A)
    records[('fixed_control_expression',1.,0,0)]=np.tile(covs['q1_seed0_n400'],(len(targets),1))
    records[('constant_tie',1.,0,0)]=np.zeros_like(A)
    raw_contexts=original[original.adjustment=='raw'][CONTEXT+['evaluation_file']].drop_duplicates().rename(columns={'evaluation_file':'original_raw_file'})
    edr_contexts=original[original.adjustment=='EDR'][CONTEXT+['evaluation_file']].drop_duplicates().rename(columns={'evaluation_file':'original_EDR_file'})
    contexts=raw_contexts.merge(edr_contexts,on=CONTEXT,validate='one_to_one',sort=False)
    contexts['original_order']=contexts.original_raw_file.str.extract(r'(\d+)\.npz')[0].astype(int)
    contexts=contexts.sort_values('original_order').reset_index(drop=True)
    contexts.insert(0,'context_id',[f'context_{i+1:04d}' for i in range(len(contexts))])
    if len(contexts)!=cfg['expected_contexts'] or original.evaluation_file.nunique()!=cfg['expected_original_evaluations']:raise ValueError('Original context coverage differs')
    if len(original)!=cfg['expected_original_evaluations']*len(targets):raise ValueError('Original target row count differs')
    original_groups={name:g.set_index('target').loc[targets] for name,g in original.groupby('evaluation_file',sort=False)}
    if any(len(g)!=len(targets) or not np.array_equal(g.fold.to_numpy(),folds) for g in original_groups.values()):raise ValueError('Original target/fold order invalid')
    code_paths=[Path(__file__),HERE/'carf_core.py',HERE/'io_utils.py',HERE/'tests/test_posthoc_diagnostic_sensitivity.py']
    code_hashes={str(p.relative_to(HERE)):sha256(p) for p in code_paths}
    OUT.mkdir(parents=True);(OUT/'contexts').mkdir()
    write_json(OUT/'execution_start.json',dict(started_utc=start,config_sha256=config_hash,freeze_receipt_sha256=sha256(FREEZE),code_sha256=code_hashes))
    contexts.to_csv(OUT/'context_index.csv',index=False)
    rows=[];strata_rows=[];checks=[];support_rows=[];baseline={};weights_by_context={};frames_by_context={}
    k=cfg['top_k']
    for no,ctx in enumerate(contexts.to_dict('records')):
        key=(ctx['model'],float(ctx['model_q']),int(ctx['model_thinning_seed']),int(ctx['algorithm_seed']))
        raw=records[key]
        cov_seed=0 if ctx['covariate_q']==1 else int(ctx['thinning_seed'])
        cov_key=f"q{ctx['covariate_q']:g}_seed{cov_seed}_n{int(ctx['n_covariate'])}"
        cov=covs[cov_key]
        scores_all=[];weights_all=[];v_all=[];supports=[]
        for ranking in RANKINGS:
            scores=ranking_scores(raw,cov,ranking)
            w=np.vstack([topk_weights(s,k) for s in scores])
            v,supported,info,W,H,bins,counts=match_background(w,cov,k,cfg['n_strata'],cfg['capacity_tolerance'])
            calculated=evaluate_weights(w,v,supported,B,C,null,folds,k)
            f=pd.DataFrame(dict(**calculated,**info))
            f['target']=targets;f['fold']=folds
            f['raw_score_all_floored']=(raw<=1e-10).all(1)
            f['raw_score_exact_constant']=np.ptp(raw,axis=1)==0
            f['ranking_exact_constant']=np.ptp(scores,axis=1)==0
            f['ranking_unique_scores']=[len(np.unique(s)) for s in scores]
            boundaries=np.partition(scores,len(genes)-k,axis=1)[:,len(genes)-k]
            f['boundary_tie_count']=(scores==boundaries[:,None]).sum(1)
            f['fractional_boundary_gene_count']=((w>0)&(w<1)).sum(1)
            f['covariate_unique_values']=len(np.unique(cov));f['covariate_exact_constant']=bool(np.ptp(cov)==0)
            bk=(ctx['model'],ranking,int(ctx['n_covariate']))
            if ctx['kind']=='baseline':baseline[bk]=w.copy()
            if bk not in baseline:raise ValueError('Baseline not established before context')
            f['top50_agreement_to_q1_same_n']=np.minimum(w,baseline[bk]).sum(1)/k
            f['top50_agreement_to_q1_n400']=np.minimum(w,baseline[(ctx['model'],ranking,400)]).sum(1)/k
            for name in CONTEXT+['context_id']:f[name]=ctx[name]
            f['ranking']=ranking;f['covariate_key']=cov_key
            if ctx['model'] in ['training_condition_mean','fixed_control_expression','constant_tie']:
                if np.abs(f.identity_gain).max()>1e-12:raise AssertionError('Shared ranking nonzero identity')
                if f.matched_identity_gain.notna().any() and np.abs(f.matched_identity_gain.dropna()).max()>1e-12:raise AssertionError('Shared ranking nonzero matched identity')
            if ctx['model']=='constant_tie' and np.nanmax(np.abs(f.matched_total_BC))>1e-12:raise AssertionError('Constant matching should equal original weights')
            if ranking in ['raw','EDR']:
                relative=ctx[f'original_{ranking}_file'];path=evaluation/relative
                if sha256(path)!=manifests['evaluation']['output_sha256'][relative]:raise ValueError('Original evaluation file changed')
                old=np.load(path,allow_pickle=False)
                if not np.array_equal(old['targets'].astype(str),targets) or not np.array_equal(old['genes'].astype(str),genes) or not np.array_equal(old['folds'],folds):raise ValueError('Original evaluation identity mismatch')
                exact_weights=np.array_equal(w,old['topk_weights'])
                score_error=float(np.max(np.abs(scores-old['scores'])))
                source_frame=original_groups[relative]
                errors={m:float(np.max(np.abs(f[m].to_numpy()-source_frame[m].to_numpy()))) for m in ['total_BC','other_target_BC','identity_gain','null_response_mean','null_response_abs_mean','top50_agreement_to_q1_same_n','top50_agreement_to_q1_n400']}
                if not exact_weights or score_error>1e-12 or max(errors.values())>1e-12:raise AssertionError('Original raw/EDR reconstruction differs')
                checks.append(dict(context_id=ctx['context_id'],ranking=ranking,original_evaluation_file=relative,exact_top50_weights=exact_weights,score_max_error=score_error,**{m+'_max_error':v for m,v in errors.items()}))
            rows.append(f);frames_by_context[(ctx['context_id'],ranking)]=f;weights_by_context[(ctx['context_id'],ranking)]=w
            scores_all.append(scores);weights_all.append(w);v_all.append(v);supports.append(supported)
            support_rows.append(dict(**{name:ctx[name] for name in CONTEXT+['context_id']},ranking=ranking,
                targets=len(targets),supported_targets=int(supported.sum()),support_rate=float(supported.mean()),
                full_fold_supported_targets=int(f.matching_fold_supported.sum()),full_fold_support_rate=float(f.matching_fold_supported.mean()),
                unsupported_folds='|'.join(str(x) for x in np.unique(folds[~f.matching_fold_supported.to_numpy()])),
                all_floor_targets=int(f.raw_score_all_floored.sum()),constant_rankings=int(f.ranking_exact_constant.sum())))
            for t in range(len(targets)):
                for b in range(cfg['n_strata']):
                    strata_rows.append(dict(context_id=ctx['context_id'],ranking=ranking,target=targets[t],fold=int(folds[t]),stratum=b,
                        stratum_genes=int(counts[b]),selected_mass=W[t,b],available_capacity=H[t,b],capacity_slack=H[t,b]-W[t,b],
                        capacity_deficit=max(0.,W[t,b]-H[t,b]),ranking_supported=bool(supported[t]),
                        matched_mass=float(v[t,bins==b].sum()) if supported[t] else np.nan))
        np.savez_compressed(OUT/'contexts'/f"{ctx['context_id']}.npz",targets=targets,genes=genes,folds=folds,
                            rankings=np.array(RANKINGS),scores=np.stack(scores_all),weights=np.stack(weights_all),
                            matched_weights=np.stack(v_all),matching_supported=np.stack(supports),strata=bins,covariate=cov)
        if (no+1)%32==0:print(f"Evaluated {no+1}/{len(contexts)} full contexts",flush=True)
    frame=pd.concat(rows,ignore_index=True)
    if len(frame)!=cfg['expected_new_evaluations']*len(targets):raise AssertionError('Incomplete evaluation grid')
    frame.to_csv(OUT/'per_target_metrics.csv.gz',index=False)
    pd.DataFrame(strata_rows).to_csv(OUT/'matching_strata_per_target.csv.gz',index=False)
    pd.DataFrame(checks).to_csv(OUT/'original_reconstruction_checks.csv',index=False)
    pd.DataFrame(support_rows).to_csv(OUT/'matching_support_by_context.csv',index=False)
    target_summary=summarize(frame,GROUP+['target','fold'],METRICS)
    target_summary.to_csv(OUT/'per_target_seed_averaged.csv.gz',index=False)
    summarize(target_summary,GROUP,METRICS).to_csv(OUT/'condition_summary.csv',index=False)
    summarize(target_summary,GROUP+['fold'],METRICS).to_csv(OUT/'fold_summary.csv',index=False)
    summarize(frame,GROUP+['thinning_seed'],METRICS).to_csv(OUT/'seed_summary.csv',index=False)
    contrasts=[]
    ckeys=['model','q','thinning_seed','n_covariate','algorithm_seed']
    for key,cs in contexts[contexts.kind=='factorial'].groupby(ckeys):
        context_by_factor=dict(zip(cs.factor,cs.context_id))
        if set(context_by_factor)!={'S0_C0','S0_Cq','Sq_C0','Sq_Cq'}:raise AssertionError('Incomplete four-arm group')
        for ranking in RANKINGS:
            for effect,l,r in [('C_at_S0','S0_C0','S0_Cq'),('C_at_Sq','Sq_C0','Sq_Cq'),('S_at_C0','S0_C0','Sq_C0'),('joint','S0_C0','Sq_Cq')]:
                lc,rc=context_by_factor[l],context_by_factor[r]
                left=frames_by_context[(lc,ranking)];right=frames_by_context[(rc,ranking)]
                d=pd.DataFrame(dict(target=targets,fold=folds,
                    top50_agreement=np.minimum(weights_by_context[(lc,ranking)],weights_by_context[(rc,ranking)]).sum(1)/k))
                for metric in ['total_BC','identity_gain','matched_total_BC','matched_identity_gain']:
                    d[metric+'_change']=right[metric].to_numpy()-left[metric].to_numpy()
                    d[metric+'_absolute_change']=np.abs(d[metric+'_change'])
                d['matching_supported_both']=left.matching_supported.to_numpy()&right.matching_supported.to_numpy()
                d['matching_fold_supported_both']=left.matching_fold_supported.to_numpy()&right.matching_fold_supported.to_numpy()
                d['all_floor_either']=left.raw_score_all_floored.to_numpy()|right.raw_score_all_floored.to_numpy()
                d['constant_either']=left.ranking_exact_constant.to_numpy()|right.ranking_exact_constant.to_numpy()
                for name,val in zip(ckeys,key):d[name]=val
                d['ranking']=ranking;d['effect']=effect;contrasts.append(d)
    contrast=pd.concat(contrasts,ignore_index=True);contrast.to_csv(OUT/'factorial_simple_effects_per_target.csv.gz',index=False)
    cm=['top50_agreement']+[m+s for m in ['total_BC','identity_gain','matched_total_BC','matched_identity_gain'] for s in ['_change','_absolute_change']]
    contrast_keys=['model','q','n_covariate','algorithm_seed','ranking','effect']
    ct=contrast.groupby(contrast_keys+['target','fold'],dropna=False)[cm].agg(strict_mean).reset_index()
    ct.groupby(contrast_keys,dropna=False)[cm].agg(strict_mean).reset_index().to_csv(OUT/'factorial_simple_effects_summary.csv',index=False)
    contrast.groupby(contrast_keys+['thinning_seed'],dropna=False)[cm].agg(strict_mean).reset_index().to_csv(OUT/'factorial_simple_effects_seed_summary.csv',index=False)
    ct.groupby(contrast_keys+['fold'],dropna=False)[cm].agg(strict_mean).reset_index().to_csv(OUT/'factorial_simple_effects_fold_summary.csv',index=False)
    decomp=[]
    for metric in ['total_BC','identity_gain','matched_total_BC','matched_identity_gain']:
        p=frame[frame.kind=='factorial'].pivot(index=ckeys+['ranking','target','fold'],columns='factor',values=metric)
        d=p.reset_index()[ckeys+['ranking','target','fold']].copy();d['metric']=metric
        for name,value in [('covariate_at_S0',p.S0_Cq-p.S0_C0),('model_at_C0',p.Sq_C0-p.S0_C0),
                           ('interaction',p.Sq_Cq-p.Sq_C0-p.S0_Cq+p.S0_C0),('joint',p.Sq_Cq-p.S0_C0)]:d[name]=value.to_numpy()
        decomp.append(d)
    pd.concat(decomp).to_csv(OUT/'factorial_decomposition_per_target.csv.gz',index=False)
    pairs=[]
    for left,right in [('raw','EDR'),('raw','rank_OLS'),('EDR','rank_OLS')]:
        a=frame[frame.ranking==left];b=frame[frame.ranking==right]
        p=b.merge(a[['context_id','target','fold']+METRICS],on=['context_id','target','fold'],suffixes=('','_left'),validate='one_to_one')
        d=p[CONTEXT+['context_id','target','fold']].copy();d['left_ranking']=left;d['right_ranking']=right
        for metric in METRICS:d[metric+'_change']=p[metric]-p[metric+'_left']
        pairs.append(d)
    pd.concat(pairs).to_csv(OUT/'paired_ranking_differences.csv.gz',index=False)
    budget_keys=[x for x in CONTEXT if x!='n_covariate']+['ranking','target','fold']
    large=frame[frame.n_covariate==400];small=frame[frame.n_covariate==64]
    p=small.merge(large[budget_keys+METRICS],on=budget_keys,suffixes=('','_n400'),validate='one_to_one')
    d=p[budget_keys].copy()
    for metric in METRICS:d[metric+'_n64_minus_n400']=p[metric]-p[metric+'_n400']
    d.to_csv(OUT/'paired_covariate_budget_differences.csv.gz',index=False)
    if sha256(CONFIG)!=config_hash or sha256(design)!=receipt['design_sha256']:raise ValueError('Rules changed during execution')
    for p in code_paths:
        if sha256(p)!=code_hashes[str(p.relative_to(HERE))]:raise ValueError('Code changed during execution')
    for name,expected in inputs.items():
        if sha256(HERE/name)!=expected:raise ValueError('Original input changed during execution')
    output={str(p.relative_to(OUT)):sha256(p) for p in OUT.rglob('*') if p.is_file()}
    write_json(OUT/'manifest.json',dict(status='completed_posthoc_diagnostic_sensitivity_no_model_refit',
        started_utc=start,finished_utc=datetime.now(timezone.utc).isoformat(),config_sha256=config_hash,
        freeze_receipt_sha256=sha256(FREEZE),source_input_sha256=inputs,code_sha256=code_hashes,
        contexts=len(contexts),ranking_conditions=len(contexts)*len(RANKINGS),target_records=len(frame),
        original_reconstructed_conditions=len(checks),matching_supported_target_records=int(frame.matching_supported.sum()),
        matching_full_fold_supported_target_records=int(frame.matching_fold_supported.sum()),
        n_targets=len(targets),n_genes=len(genes),rankings=RANKINGS,
        runtime=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),
        output_sha256=output,interpretation=cfg['status']))
    print(f'Completed {len(frame)} target records; {len(output)} output files',flush=True)


if __name__=='__main__':main()
