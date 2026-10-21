"""Independent posthoc diagnostic audit.

Main numerical reconstruction imports no study evaluator/math helpers. Synthetic
branch checks additionally invoke production functions against hand-set cases.

Reads frozen source scores/references and completed diagnostic outputs only.
Writes a separate review/ subdirectory after the main manifest is complete.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE / 'results/rpe1_measurement_validation_v1'
POST = HERE / 'results/posthoc_diagnostic_sensitivity_v1'
OUT = POST / 'review'
CONFIG = HERE / 'config/posthoc_diagnostic_sensitivity_v1.json'
RANKINGS = ['raw', 'EDR', 'rank_OLS']
CONTEXT = ['model', 'q', 'thinning_seed', 'n_covariate', 'algorithm_seed',
           'factor', 'kind', 'model_q', 'model_thinning_seed', 'covariate_q']


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def compare(actual, expected, label, tolerance=1e-12, exact=False):
    a, b = np.asarray(actual), np.asarray(expected)
    if a.shape != b.shape:
        raise AssertionError(f'{label}: shape {a.shape} != {b.shape}')
    if exact:
        if not np.array_equal(a, b, equal_nan=a.dtype.kind in 'fc'):
            raise AssertionError(f'{label}: exact equality failed')
        return 0.0
    if not np.allclose(a, b, atol=tolerance, rtol=0, equal_nan=True):
        finite = np.isfinite(a) & np.isfinite(b)
        delta = np.max(abs(a[finite].astype(float) - b[finite].astype(float))) if finite.any() else None
        raise AssertionError(f'{label}: numeric/NA disagreement (finite max {delta})')
    finite = np.isfinite(a) & np.isfinite(b)
    return float(np.max(abs(a[finite].astype(float) - b[finite].astype(float)))) if finite.any() else 0.0


def normalized_midrank(row):
    """Midranks calculated from unique-group counts, not scipy.rankdata."""
    _, inverse, counts = np.unique(row, return_inverse=True, return_counts=True)
    before = np.cumsum(counts) - counts
    midrank = before + (counts + 1) / 2
    return (midrank[inverse] - .5) / len(row)


def ols_residuals(y, x):
    """Closed-form centered one-regressor OLS, independent of study lstsq."""
    y = np.asarray(y, float)
    centered_y = y - y.mean(axis=1, keepdims=True)
    centered_x = np.asarray(x, float) - np.mean(x)
    if np.std(x) > 1e-12:
        slope = (centered_y @ centered_x) / np.dot(centered_x, centered_x)
        residual = centered_y - slope[:, None] * centered_x[None, :]
    else:
        residual = centered_y
    guard = np.ptp(residual, axis=1) < 1e-12 * np.maximum(1, np.abs(y).max(axis=1))
    residual[guard] = 0
    return residual


def all_scores(raw, cov):
    edr = ols_residuals(np.log10(np.maximum(raw, 1e-10)), cov)
    xr = normalized_midrank(cov)
    yr = np.vstack([normalized_midrank(row) for row in raw])
    rank_ols = ols_residuals(yr, xr)
    return np.stack([raw, edr, rank_ols]), xr


def inclusion(scores, k=50):
    """Descending full sort, preserving exact boundary ties."""
    out = []
    for row in scores:
        boundary = np.sort(row)[-k]
        above, tied = row > boundary, row == boundary
        out.append(above.astype(float) + tied * ((k - above.sum()) / tied.sum()))
    return np.asarray(out)


def matched_background(w, strata, capacity_tolerance=1e-12):
    """Independent whole-row failure propagation from five bin capacities."""
    v = np.zeros_like(w)
    selected = np.zeros((len(w), 5))
    available = np.zeros_like(selected)
    for bin_id in range(5):
        mask = strata == bin_id
        selected[:, bin_id] = w[:, mask].sum(axis=1)
        available[:, bin_id] = (1 - w[:, mask]).sum(axis=1)
        use = (selected[:, bin_id] > 0) & (available[:, bin_id] > 0)
        rows = np.flatnonzero(use)
        v[np.ix_(rows, mask)] = ((selected[rows, bin_id] / available[rows, bin_id])[:, None]
                                * (1 - w[np.ix_(rows, mask)]))
    supported = ~((selected > 0) & (selected - available > capacity_tolerance)).any(axis=1)
    v[~supported] = np.nan
    return v, supported, selected, available


def utilities(weights, response_product, null_product, folds, support=None, k=50):
    """Use complete fold sums; never reduce the fixed mismatched ranking pool."""
    own = np.sum(weights * response_product, axis=1) / k
    other = np.full(len(weights), np.nan)
    fold_support = np.ones(len(weights), bool)
    for fold in np.unique(folds):
        ids = np.flatnonzero(folds == fold)
        if len(ids) < 2:
            raise AssertionError('Fold lacks mismatched source rankings')
        if support is not None and not np.all(support[ids]):
            fold_support[ids] = False
            continue
        summed_sources = weights[ids].sum(axis=0)
        other[ids] = np.sum(response_product[ids] * (summed_sources - weights[ids]), axis=1) / ((len(ids) - 1) * k)
    null = null_product @ weights.T / k
    return dict(total_BC=own, other_target_BC=other, identity_gain=own-other,
                null_response_mean=null.mean(axis=0), null_response_abs_mean=np.abs(null).mean(axis=0),
                full_fold_supported=fold_support), null


def synthetic_checks():
    compare(normalized_midrank(np.array([5, 1, 1, 9.])), [.625, .25, .25, .875], 'hand midranks', exact=True)
    row = np.array([[1, 2, 2, 2, 0.]])
    compare(inclusion(row, 2), [[0, 2/3, 2/3, 2/3, 0]], 'hand fractional ties', exact=True)
    compare(ols_residuals(np.array([[1, 2, 3.], [7, 7, 7]]), [2, 4, 6]), np.zeros((2, 3)), 'affine/constant OLS', exact=True)
    compare(ols_residuals(np.array([[1, 2, 3.]]), [5, 5, 5]), [[-1, 0, 1]], 'constant covariate', exact=True)
    strata = np.array([0, 0, 0, 1, 1, 1, 1, 1])
    w = np.array([[1, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 1, 0, 0]], float)
    v, supported, selected, available = matched_background(w, strata)
    compare(supported, [False, True, True, True], 'capacity deficit support', exact=True)
    assert np.isnan(v[0]).all() and np.isfinite(v[1:]).all()
    compare(v[1], [0, .5, .5, 0, .25, .25, .25, .25], 'hand matching weights', exact=True)
    product = np.arange(32.).reshape(4, 8)
    result, _ = utilities(w-v, product, np.ones((2, 8)), np.zeros(4, int), supported, k=2)
    assert np.isnan(result['identity_gain']).all() and np.isfinite(result['total_BC'][1:]).all()
    assert not result['full_fold_supported'].any()
    constant_w = np.full((4, 8), .25)
    cv, cs, _, _ = matched_background(constant_w, strata)
    compare(cv, constant_w, 'constant background cancellation')
    shared, _ = utilities(constant_w-cv, product, np.ones((2, 8)), np.zeros(4, int), cs, k=2)
    compare(shared['identity_gain'], np.zeros(4), 'constant matched identity')
    # Check the actual production failure branch independently of its own tests.
    import run_posthoc_diagnostic_sensitivity as production
    # Tied covariates place 3 genes in one stratum and 5 in another.
    artificial_cov = np.array([0., 0., 0., 1., 1., 1., 1., 1.])
    actual_v, actual_support, _, _, _, _, _ = production.match_background(w, artificial_cov, k=2)
    compare(actual_v, v, 'production capacity/whole-row NA branch')
    compare(actual_support, supported, 'production matching support', exact=True)
    actual = production.evaluate_weights(w, actual_v, actual_support, product, np.ones_like(product), np.ones((2,8)), np.zeros(4,int), k=2)
    compare(actual['matched_identity_gain'], np.full(4,np.nan), 'production whole-fold NA propagation')
    compare(actual['matched_total_BC'], result['total_BC'], 'production supported own utility retained')
    shared_w = np.tile(w[1], (4,1))
    shared_v, shared_support, _, _, _, _, _ = production.match_background(shared_w, artificial_cov, k=2)
    actual = production.evaluate_weights(shared_w, shared_v, shared_support, product, np.ones_like(product), np.ones((2,8)), np.zeros(4,int), k=2)
    compare(actual['identity_gain'],np.zeros(4),'production shared negative unmatched identity')
    compare(actual['matched_identity_gain'],np.zeros(4),'production shared negative matched identity')
    return dict(independent_cases_passed=10,
                cases=['hand midranks', 'fractional ties', 'affine and constant OLS', 'constant covariate',
                       'matching capacity deficit', 'hand supported matching', 'failure propagates to whole fold',
                       'constant ranking matched cancellation', 'production failure branch independently exercised',
                       'production shared-ranking negative independently exercised'])

METRICS = ['total_BC', 'other_target_BC', 'identity_gain', 'null_response_mean', 'null_response_abs_mean',
           'matched_background_BC', 'matched_total_BC', 'matched_other_target_BC', 'matched_identity_gain',
           'matched_null_response_mean', 'matched_null_response_abs_mean',
           'top50_agreement_to_q1_same_n', 'top50_agreement_to_q1_n400',
           'expression_mean_difference', 'expression_mean_absolute_difference',
           'expression_rank_mean_difference', 'expression_rank_mean_absolute_difference']
GROUP = ['model', 'q', 'n_covariate', 'algorithm_seed', 'factor', 'kind', 'ranking']


def strict_groups(frame, keys, metrics):
    groups = frame.groupby(keys, dropna=False, sort=True)
    mean = groups[metrics].mean()
    complete = groups[metrics].count().eq(groups.size(), axis='index')
    return mean.where(complete).reset_index()


def compare_table(expected, stored, keys, metrics, label):
    if len(expected) != len(stored):
        raise AssertionError(f'{label}: number of rows differs')
    joined = expected.merge(stored, on=keys, suffixes=('_independent', '_saved'), validate='one_to_one')
    if len(joined) != len(expected):
        raise AssertionError(f'{label}: keys differ')
    return max(compare(joined[m+'_independent'], joined[m+'_saved'], f'{label}: {m}') for m in metrics)


def main():
    start = datetime.now(timezone.utc).isoformat()
    if not (POST/'manifest.json').exists():
        raise SystemExit('Completed main manifest required before independent review')
    if OUT.exists():
        raise SystemExit('Review directory already exists; refuse overwriting review evidence')
    manifest_hash = sha(POST/'manifest.json')
    manifest = json.loads((POST/'manifest.json').read_text())
    cfg = json.loads(CONFIG.read_text())
    if manifest['status'] != 'completed_posthoc_diagnostic_sensitivity_no_model_refit':
        raise AssertionError('Incomplete/wrong posthoc status')
    if manifest['config_sha256'] != sha(CONFIG):
        raise AssertionError('Configuration does not match completed run')
    for name, expected in manifest['output_sha256'].items():
        if sha(POST/name) != expected:
            raise AssertionError(f'Posthoc output hash changed: {name}')
    for name, expected in manifest['source_input_sha256'].items():
        if sha(HERE/name) != expected:
            raise AssertionError(f'Original source hash changed: {name}')
    for name, expected in manifest['code_sha256'].items():
        if sha(HERE/name) != expected:raise AssertionError('Executed posthoc code changed: '+name)
    code_hash = sha(__file__)
    synthetic = synthetic_checks()
    refs = np.load(BASE/'preparation/evaluation_responses.npz', allow_pickle=False)
    targets, genes, folds = refs['targets'].astype(str), refs['genes'].astype(str), refs['folds'].astype(int)
    A, product, null = refs['A'], refs['B']*refs['C'], refs['null_B']*refs['null_C']
    if product.shape != (128,665) or null.shape != (30,665):
        raise AssertionError('Reference dimensions differ')
    covs = np.load(BASE/'evaluation/covariates.npz', allow_pickle=False)
    records = {}
    for r in pd.read_csv(BASE/'models_safe_export/model_index.csv').to_dict('records'):
        p = BASE/'models_safe_export'/r['path']
        if sha(p) != r['sha256']:
            raise AssertionError('Safe model score index changed')
        with np.load(p, allow_pickle=False) as z:
            compare(z['targets'].astype(str),targets,'source target order',exact=True)
            compare(z['genes'].astype(str),genes,'source gene order',exact=True)
            compare(z['folds'],folds,'source target folds',exact=True)
            records[(r['model'],r['q'],r['thinning_seed'],r['algorithm_seed'])] = z['raw_scores']
    if len(records) != 35:
        raise AssertionError('Score model coverage changed')
    records[('observed_A_technical_anchor',1.,0,0)] = abs(A)
    records[('fixed_control_expression',1.,0,0)] = np.tile(covs['q1_seed0_n400'],(128,1))
    records[('constant_tie',1.,0,0)] = np.zeros_like(A)
    old = pd.read_csv(BASE/'evaluation/per_target_metrics.csv.gz')
    contexts = pd.read_csv(POST/'context_index.csv')
    frame = pd.read_csv(POST/'per_target_metrics.csv.gz')
    strata_table = pd.read_csv(POST/'matching_strata_per_target.csv.gz')
    support_table = pd.read_csv(POST/'matching_support_by_context.csv')
    if len(contexts)!=256 or len(frame)!=98304 or len(strata_table)!=491520 or len(support_table)!=768:
        raise AssertionError('Posthoc enumeration incomplete')
    if not frame.groupby(['context_id','ranking']).size().eq(128).all():
        raise AssertionError('Missing targets in a context/ranking')
    if set(frame.ranking) != set(RANKINGS) or frame.duplicated(['context_id','ranking','target']).any():
        raise AssertionError('Ranking coverage/target uniqueness differs')
    expected_contexts = old[old.adjustment=='raw'][CONTEXT+['evaluation_file']].drop_duplicates()
    mapped = contexts.merge(expected_contexts,on=CONTEXT,validate='one_to_one')
    if len(mapped)!=256 or not (mapped.original_raw_file==mapped.evaluation_file).all():
        raise AssertionError('Original raw context mapping differs')
    expected_edr = old[old.adjustment=='EDR'][CONTEXT+['evaluation_file']].drop_duplicates()
    mapped = contexts.merge(expected_edr,on=CONTEXT,validate='one_to_one')
    if len(mapped)!=256 or not (mapped.original_EDR_file==mapped.evaluation_file).all():
        raise AssertionError('Original EDR context mapping differs')
    old_groups = {k:g.set_index('target').loc[targets] for k,g in old.groupby('evaluation_file',sort=False)}
    groups = {(c,r):g.set_index('target').loc[targets] for (c,r),g in frame.groupby(['context_id','ranking'],sort=False)}
    strata_groups = {(c,r):g.set_index(['target','stratum']).loc[pd.MultiIndex.from_product([targets,range(5)])]
                     for (c,r),g in strata_table.groupby(['context_id','ranking'],sort=False)}
    supports = support_table.set_index(['context_id','ranking'])
    baselines = {}; cached_weights = {}; summaries = []; errors = {}; raw_conditions = {}
    def record(label, value): errors[label] = max(errors.get(label,0.),float(value))
    for number, ctx in enumerate(contexts.to_dict('records')):
        cid = ctx['context_id']
        raw_key = (ctx['model'],ctx['model_q'],ctx['model_thinning_seed'],ctx['algorithm_seed'])
        raw = records[raw_key]
        cov_seed = 0 if ctx['covariate_q']==1 else int(ctx['thinning_seed'])
        ck = f"q{ctx['covariate_q']:g}_seed{cov_seed}_n{ctx['n_covariate']}"
        cov = covs[ck]; scores, xr = all_scores(raw,cov)
        strata = np.minimum(4,np.floor(5*xr).astype(int))
        with np.load(POST/'contexts'/f'{cid}.npz',allow_pickle=False) as z:
            for name,expected in [('targets',targets),('genes',genes),('folds',folds),('rankings',RANKINGS),('strata',strata),('covariate',cov)]:
                compare(z[name],expected,cid+' '+name,exact=True)
            record('independent_all_scores_max_error',compare(scores,z['scores'],cid+' scores'))
            for j, ranking in enumerate(RANKINGS):
                g = groups[(cid,ranking)]
                compare(g.fold,folds,'metric fold identity',exact=True)
                for name in CONTEXT:
                    if not g[name].eq(ctx[name]).all():raise AssertionError('Metric/context metadata mismatch')
                w = inclusion(scores[j]); v,supported,W,H = matched_background(w,strata)
                compare(w,z['weights'][j],'exact independent Top50 weights',exact=True)
                record('background_weights_max_error',compare(v,z['matched_weights'][j],'independent background'))
                compare(supported,z['matching_supported'][j],'matching support',exact=True)
                if supported.any():
                    if not np.all((v[supported]>=-1e-12)&((v+w)[supported]<=1+1e-12)):
                        raise AssertionError('Matching violates complementary inclusion capacity')
                if not np.isnan(v[~supported]).all():raise AssertionError('Failed matching retained background entries')
                u,nu = utilities(w,product,null,folds)
                mu,mnu = utilities(w-v,product,null,folds,supported)
                expectation = {name:u[name] for name in ['total_BC','other_target_BC','identity_gain','null_response_mean','null_response_abs_mean']}
                expectation.update(matched_background_BC=np.sum(v*product,axis=1)/50)
                expectation.update({'matched_'+name:mu[name] for name in ['total_BC','other_target_BC','identity_gain','null_response_mean','null_response_abs_mean']})
                expectation['matching_fold_supported'] = mu['full_fold_supported']
                expectation['matching_supported'] = supported
                diff = (w-v)@cov/50;rank_diff=(w-v)@xr/50
                expectation.update(expression_mean_difference=diff,expression_mean_absolute_difference=abs(diff),
                    expression_rank_mean_difference=rank_diff,expression_rank_mean_absolute_difference=abs(rank_diff),
                    matching_min_capacity_slack=np.where(W>0,H-W,np.inf).min(1),
                    matching_max_capacity_deficit=np.maximum(W-H,0).max(1),
                    matching_failed_strata=((W>0)&(W-H>1e-12)).sum(1),
                    matching_stratum_mass_max_error=np.abs(np.column_stack([v[:,strata==b].sum(1)-W[:,b] for b in range(5)])).max(1),
                    matching_total_mass_error=v.sum(1)-50,
                    raw_score_all_floored=(raw<=1e-10).all(1),raw_score_exact_constant=np.ptp(raw,axis=1)==0,
                    ranking_exact_constant=np.ptp(scores[j],axis=1)==0,
                    ranking_unique_scores=np.array([len(np.unique(row)) for row in scores[j]]),
                    boundary_tie_count=(scores[j]==np.sort(scores[j],axis=1)[:,-50,None]).sum(1),
                    fractional_boundary_gene_count=((w>0)&(w<1)).sum(1))
                bk = (ctx['model'],ranking,ctx['n_covariate'])
                if ctx['kind']=='baseline':baselines[bk]=w.copy()
                expectation['top50_agreement_to_q1_same_n']=np.minimum(w,baselines[bk]).sum(1)/50
                expectation['top50_agreement_to_q1_n400']=np.minimum(w,baselines[(ctx['model'],ranking,400)]).sum(1)/50
                for name,value in expectation.items():
                    record('per_target_'+name+'_max_error',compare(value,g[name],cid+' '+ranking+' '+name))
                if not g.covariate_unique_values.eq(len(np.unique(cov))).all() or not g.covariate_exact_constant.eq(np.ptp(cov)==0).all():
                    raise AssertionError('Covariate ties/constant flag mismatch')
                if not g.covariate_key.eq(ck).all():raise AssertionError('Wrong covariate key')
                sg=strata_groups[(cid,ranking)]
                for name,value in [('stratum_genes',np.tile(np.bincount(strata,minlength=5),(128,1))),('selected_mass',W),('available_capacity',H),('capacity_slack',H-W),('capacity_deficit',np.maximum(W-H,0)),('matched_mass',np.column_stack([v[:,strata==b].sum(1) for b in range(5)]))]:
                    record('stratum_'+name+'_max_error',compare(value.reshape(-1),sg[name],'all stratum rows '+name))
                compare(sg.ranking_supported,np.repeat(supported,5),'stratum support propagation',exact=True)
                sr=supports.loc[(cid,ranking)]
                for name,value in [('targets',128),('supported_targets',supported.sum()),('support_rate',supported.mean()),('full_fold_supported_targets',mu['full_fold_supported'].sum()),('full_fold_support_rate',mu['full_fold_supported'].mean()),('all_floor_targets',(raw<=1e-10).all(1).sum()),('constant_rankings',(np.ptp(scores[j],axis=1)==0).sum())]:
                    compare(sr[name],value,'context support table '+name)
                if ranking=='raw':
                    if raw_key in raw_conditions:compare(w,raw_conditions[raw_key],'raw C-only invariance',exact=True)
                    else:raw_conditions[raw_key]=w.copy()
                if ranking in ['raw','EDR']:
                    oldfile=ctx['original_'+ranking+'_file'];og=old_groups[oldfile]
                    with np.load(BASE/'evaluation'/oldfile,allow_pickle=False) as original_npz:
                        record('original_raw_EDR_score_max_error',compare(scores[j],original_npz['scores'],'original scores'))
                        compare(w,original_npz['topk_weights'],'original exact Top50 weights',exact=True)
                        record('original_null_array_max_error',compare(nu,original_npz['null_utility'],'original null arrays'))
                        original_u=original_npz['within_fold_utility']
                        for fold in np.unique(folds):
                            ids=np.flatnonzero(folds==fold)
                            expected_u=product[ids]@w[ids].T/50
                            record('original_full_fold_matrix_max_error',compare(expected_u,original_u[np.ix_(ids,ids)],'original full fold utility'))
                        if not np.isnan(original_u[folds[:,None]!=folds[None,:]]).all():raise AssertionError('Original crossfold source included')
                    for name in ['total_BC','other_target_BC','identity_gain','null_response_mean','null_response_abs_mean','top50_agreement_to_q1_same_n','top50_agreement_to_q1_n400']:
                        record('original_'+name+'_max_error',compare(expectation[name],og[name],'original CSV '+name))
                if ctx['model'] in ['training_condition_mean','fixed_control_expression','constant_tie']:
                    compare(u['identity_gain'],np.zeros(128),'shared ranking unmatched identity')
                    compare(mu['identity_gain'][mu['full_fold_supported']],np.zeros(mu['full_fold_supported'].sum()),'shared ranking matched identity')
                if ctx['model']=='constant_tie':compare(mu['total_BC'],np.zeros(128),'constant matched utility')
                cached_weights[(cid,ranking)]=w
                summaries.append(dict(context_id=cid,ranking=ranking,supported_targets=int(supported.sum()),
                    supported_fold_targets=int(mu['full_fold_supported'].sum()),minimum_stratum_genes=int(np.bincount(strata,minlength=5).min()),
                    minimum_occupied_capacity_slack=float(np.where(W>0,H-W,np.inf).min()),
                    all_floor_targets=int((raw<=1e-10).all(1).sum()),constant_rankings=int((np.ptp(scores[j],axis=1)==0).sum())))
        if (number+1)%32==0:print(f'Independently checked {number+1}/256 contexts',flush=True)
    # Strict seed/target/fold aggregation, implemented using completeness masks.
    per_target = strict_groups(frame,GROUP+['target','fold'],METRICS)
    record('seed_averaged_table_max_error',compare_table(per_target,pd.read_csv(POST/'per_target_seed_averaged.csv.gz'),GROUP+['target','fold'],METRICS,'seed means'))
    for filename,source,keys in [('condition_summary.csv',per_target,GROUP),('fold_summary.csv',per_target,GROUP+['fold']),('seed_summary.csv',frame,GROUP+['thinning_seed'])]:
        expected = strict_groups(source,keys,METRICS)
        record(filename+'_max_error',compare_table(expected,pd.read_csv(POST/filename),keys,METRICS,filename))
    # All paired ranking differences; no favorable pair or target subset.
    pairs=pd.read_csv(POST/'paired_ranking_differences.csv.gz')
    if len(pairs)!=98304:raise AssertionError('Paired-ranking table incomplete')
    for (left,right),pg in pairs.groupby(['left_ranking','right_ranking']):
        a=frame[frame.ranking==left].set_index(['context_id','target','fold'])
        b=frame[frame.ranking==right].set_index(['context_id','target','fold'])
        expected=(b[METRICS]-a[METRICS]).rename(columns={m:m+'_change' for m in METRICS}).reset_index()
        record('paired_ranking_max_error',compare_table(expected,pg,['context_id','target','fold'],[m+'_change' for m in METRICS],'ranking pair'))
    budget_keys=[x for x in CONTEXT if x!='n_covariate']+['ranking','target','fold']
    small=frame[frame.n_covariate==64].set_index(budget_keys);large=frame[frame.n_covariate==400].set_index(budget_keys)
    expected=(small[METRICS]-large[METRICS]).rename(columns={m:m+'_n64_minus_n400' for m in METRICS}).reset_index()
    record('paired_budget_max_error',compare_table(expected,pd.read_csv(POST/'paired_covariate_budget_differences.csv.gz'),budget_keys,[m+'_n64_minus_n400' for m in METRICS],'budget pair'))
    # Independently enumerate every four-condition contrast and interaction.
    ckeys=['model','q','thinning_seed','n_covariate','algorithm_seed']
    contrasts=pd.read_csv(POST/'factorial_simple_effects_per_target.csv.gz')
    primary=['total_BC','identity_gain','matched_total_BC','matched_identity_gain']
    expected_contrasts=[]
    for key,cg in contexts[contexts.kind=='factorial'].groupby(ckeys):
        ids=dict(zip(cg.factor,cg.context_id))
        if set(ids)!=set(['S0_C0','S0_Cq','Sq_C0','Sq_Cq']):raise AssertionError('Four-arm scope incomplete')
        for ranking in RANKINGS:
            for effect,left,right in [('C_at_S0','S0_C0','S0_Cq'),('C_at_Sq','Sq_C0','Sq_Cq'),('S_at_C0','S0_C0','Sq_C0'),('joint','S0_C0','Sq_Cq')]:
                a=groups[(ids[left],ranking)];b=groups[(ids[right],ranking)]
                d=pd.DataFrame(dict(target=targets,fold=folds,top50_agreement=np.minimum(cached_weights[(ids[left],ranking)],cached_weights[(ids[right],ranking)]).sum(1)/50))
                for m in primary:
                    delta=b[m].to_numpy()-a[m].to_numpy();d[m+'_change']=delta;d[m+'_absolute_change']=abs(delta)
                for k,value in zip(ckeys,key):d[k]=value
                d['ranking']=ranking;d['effect']=effect;expected_contrasts.append(d)
    ec=pd.concat(expected_contrasts,ignore_index=True);cm=['top50_agreement']+[m+s for m in primary for s in ['_change','_absolute_change']]
    record('factorial_simple_effect_max_error',compare_table(ec,contrasts,ckeys+['ranking','effect','target','fold'],cm,'factorial simple effects'))
    contrasts_keys=['model','q','n_covariate','algorithm_seed','ranking','effect']
    ct=strict_groups(ec,contrasts_keys+['target','fold'],cm)
    for filename,source,keys in [('factorial_simple_effects_summary.csv',ct,contrasts_keys),('factorial_simple_effects_fold_summary.csv',ct,contrasts_keys+['fold']),('factorial_simple_effects_seed_summary.csv',ec,contrasts_keys+['thinning_seed'])]:
        record(filename+'_max_error',compare_table(strict_groups(source,keys,cm),pd.read_csv(POST/filename),keys,cm,filename))
    decomp=pd.read_csv(POST/'factorial_decomposition_per_target.csv.gz')
    dkeys=ckeys+['ranking','target','fold']
    for metric in primary:
        p=frame[frame.kind=='factorial'].pivot(index=dkeys,columns='factor',values=metric)
        d=pd.DataFrame({'covariate_at_S0':p['S0_Cq']-p['S0_C0'],'model_at_C0':p['Sq_C0']-p['S0_C0'],
                        'interaction':p['Sq_Cq']-p['Sq_C0']-p['S0_Cq']+p['S0_C0'],'joint':p['Sq_Cq']-p['S0_C0']}).reset_index()
        record('factorial_decomposition_max_error',compare_table(d,decomp[decomp.metric==metric],dkeys,['covariate_at_S0','model_at_C0','interaction','joint'],'factorial decomposition'))
    if sha(POST/'manifest.json')!=manifest_hash or sha(__file__)!=code_hash:raise AssertionError('Review inputs/code changed')
    summary=dict(status='passed_independent_full_posthoc_audit',started_utc=start,finished_utc=datetime.now(timezone.utc).isoformat(),
        reviewed_manifest_sha256=manifest_hash,config_sha256=sha(CONFIG),review_script_sha256=code_hash,
        verified_posthoc_files=len(manifest['output_sha256']),contexts=256,ranking_conditions=768,target_rows=98304,
        independently_reconstructed_original_raw_EDR_conditions=512,matching_stratum_rows=491520,
        full_fold_source_counts={str(f):int((folds==f).sum()-1) for f in np.unique(folds)},
        unsupported_ranking_target_rows=int((~frame.matching_supported).sum()),
        unsupported_full_fold_target_rows=int((~frame.matching_fold_supported).sum()),
        numerical_max_errors=errors,synthetic_contract_checks=synthetic,
        interpretation_limits=['Posthoc sensitivity on previously inspected RPE1 outcomes, not independent confirmation.',
            'Main reconstruction uses separate midrank, closed-form OLS, sorted ties and complete-fold sum implementations.',
            'Matching capacity failures are retained; zero observed failures do not establish universal support.',
            'Tiny native score differences remain flagged; rank-OLS does not recover biological information.',
            'No new fitting, raw-data download, R provisioning, causal claim or calibrated statistical inference.'])
    OUT.mkdir()
    pd.DataFrame(summaries).to_csv(OUT/'independent_context_checks.csv',index=False)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
