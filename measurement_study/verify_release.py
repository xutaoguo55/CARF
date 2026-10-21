#!/usr/bin/env python3
"""Verify this portable review subset, rerun unit tests, and regenerate figures.

No raw download, R model fit, original workspace, or historical virtual
environment is needed. This does not claim a full raw-to-model reproduction.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.metadata
import json
import platform
import shutil
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
from public_provenance import source_identity_sha256

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def checksum_check():
    manifest=json.loads((ROOT/'release_manifest.json').read_text())
    for name,record in manifest['payload_files'].items():
        path=ROOT/name
        if not path.is_file() or path.stat().st_size!=record['bytes'] or sha(path)!=record['sha256']:
            raise ValueError(f'Package checksum mismatch: {name}')
    return manifest

def numeric_checks(out):
    import numpy as np
    import pandas as pd
    from scipy.sparse import load_npz
    base=ROOT/'data/rpe1';evaluation=base/'evaluation';models=base/'models_safe_export';prep=base/'preparation'
    em=json.loads((evaluation/'manifest.json').read_text());sm=json.loads((models/'manifest.json').read_text())
    if sha(models/'manifest.json')!=em['model_manifest_sha256'] or source_identity_sha256(ROOT,prep/'manifest.json')!=em['preparation_manifest_sha256']:
        raise AssertionError('Frozen manifest identity chain changed (including declared public metadata views)')
    if sm['source_model_manifest_sha256']!=sha(ROOT/'provenance/rpe1_original_models/manifest.json'):
        raise AssertionError('Original model-to-safe-export link changed')
    for p in evaluation.iterdir():
        if p.is_file() and p.name!='manifest.json':
            if sha(p)!=em['output_sha256'][p.name]:raise AssertionError('Frozen table hash changed')
    z=np.load(prep/'evaluation_responses.npz',allow_pickle=False)
    targets=z['targets'].astype(str);genes=z['genes'].astype(str);folds=z['folds'];A=z['A'];B=z['B'];C=z['C'];null=z['null_B']*z['null_C']
    if A.shape!=(128,665) or B.shape!=A.shape or C.shape!=A.shape:raise AssertionError('Reference scope mismatch')
    cov=np.load(evaluation/'covariates.npz',allow_pickle=False)
    scores={};numeric_members=0
    index=pd.read_csv(models/'model_index.csv')
    for r in index.to_dict('records'):
        p=models/r['path']
        if sha(p)!=r['sha256']:raise AssertionError('Score index mismatch')
        f=np.load(p,allow_pickle=False)
        for key in f.files:
            if f[key].dtype.hasobject:raise AssertionError('Unsafe object array in review scores')
            numeric_members+=int(f[key].dtype.kind in 'biufc')
        if not np.array_equal(f['targets'],targets) or not np.array_equal(f['genes'],genes) or not np.array_equal(f['folds'],folds):raise AssertionError('Score identity mismatch')
        scores[(r['model'],r['q'],r['thinning_seed'],r['algorithm_seed'])]=f['raw_scores']
    if len(scores)!=35:raise AssertionError('Missing model condition')
    scores[('observed_A_technical_anchor',1.,0,0)]=abs(A)
    scores[('fixed_control_expression',1.,0,0)]=np.tile(cov['q1_seed0_n400'],(128,1))
    scores[('constant_tie',1.,0,0)]=np.zeros_like(A)

    def equal(a,b,label):
        a=np.asarray(a);b=np.asarray(b)
        if a.shape!=b.shape or not np.allclose(a,b,atol=1e-12,rtol=1e-11,equal_nan=True):raise AssertionError(label)
        good=np.isfinite(a)&np.isfinite(b)
        return float(np.max(abs(a[good]-b[good]))) if good.any() else 0.
    def wgt(values):
        answer=[]
        for row in values:
            boundary=np.sort(row)[-50];above=row>boundary;ties=row==boundary
            answer.append(above.astype(float)+ties*((50-above.sum())/ties.sum()))
        return np.asarray(answer)
    def residual(raw,c):
        y=np.log10(np.maximum(raw,1e-10));yc=y-y.mean(1,keepdims=True);x=c-c.mean()
        result=yc-((yc@x)/(x@x))[:,None]*x[None,:] if c.std()>1e-12 else yc
        result[np.ptp(result,axis=1)<1e-12*np.maximum(1,np.max(abs(y),axis=1))]=0
        return result
    def evaluate(w):
        same=np.empty(128);other=np.empty(128)
        for fold in range(3):
            ids=np.flatnonzero(folds==fold);u=(B[ids]*C[ids])@w[ids].T/50
            same[ids]=np.diag(u);other[ids]=(u.sum(1)-np.diag(u))/(len(ids)-1)
        return same,other,null@w.T/50
    frame=pd.read_csv(evaluation/'per_target_metrics.csv.gz')
    if len(frame)!=65536 or frame.evaluation_file.nunique()!=512:raise AssertionError('Missing evaluation records')
    if not frame.groupby('evaluation_file').size().eq(128).all():raise AssertionError('Targets lost in a condition')
    # Regenerate study means from the per-target, per-seed data.
    metrics=['total_BC','identity_gain','top50_agreement_to_q1_same_n','top50_agreement_to_q1_n400','null_response_mean','null_response_abs_mean','unadjusted_prediction_MSE_B','unadjusted_prediction_MSE_C','unadjusted_prediction_Pearson_B','unadjusted_prediction_Pearson_C']
    keys=['model','q','n_covariate','algorithm_seed','factor','kind','adjustment']
    targetmeans=frame.groupby(keys+['target','fold'],dropna=False)[metrics].mean().reset_index()
    regenerated=targetmeans.groupby(keys,dropna=False)[metrics].mean().reset_index()
    original=pd.read_csv(evaluation/'condition_summary.csv');joined=regenerated.merge(original,on=keys,suffixes=('_new','_old'),validate='one_to_one')
    table_error=max(equal(joined[m+'_new'],joined[m+'_old'],'Condition-summary reproduction') for m in metrics)
    regenerated.to_csv(out/'regenerated_condition_summary.csv',index=False)
    # Recompute all raw rankings plus six representative adjusted conditions.
    raw_error=0.
    for key,g in frame[frame.adjustment=='raw'].groupby(['model','model_q','model_thinning_seed','algorithm_seed']):
        same,other,n=evaluate(wgt(scores[key]));expected=pd.DataFrame(dict(target=targets,total_BC=same,identity_gain=same-other,null_response_mean=n.mean(0)))
        joined=g.merge(expected,on='target',suffixes=('','_new'),validate='many_to_one')
        for m in ['total_BC','identity_gain','null_response_mean']:raw_error=max(raw_error,equal(joined[m],joined[m+'_new'],'Raw ranking reproduction'))
    selected=[('native_scTenifoldKnk',1.,0,400,'baseline','S0_C0',42),('heldout_linear',.25,20261113,64,'factorial','Sq_Cq',0),('control_only_Pearson',.25,20261111,64,'factorial','S0_Cq',0),('training_condition_mean',.5,20261112,400,'factorial','Sq_C0',0),('constant_tie',1.,0,64,'baseline','S0_C0',0),('observed_A_technical_anchor',1.,0,400,'baseline','S0_C0',0)]
    edr_errors=[]
    for model,q,seed,n,kind,factor,algorithm in selected:
        g=frame[(frame.model==model)&(frame.q==q)&(frame.thinning_seed==seed)&(frame.n_covariate==n)&(frame.kind==kind)&(frame.factor==factor)&(frame.algorithm_seed==algorithm)&(frame.adjustment=='EDR')]
        if len(g)!=128:raise AssertionError('Representative EDR condition absent')
        g=g.set_index('target').loc[targets];r=g.iloc[0];cs=0 if r.covariate_q==1 else seed
        c=cov[f'q{r.covariate_q:g}_seed{cs}_n{n}'];raw=scores[(model,r.model_q,r.model_thinning_seed,algorithm)]
        same,other,nul=evaluate(wgt(residual(raw,c)))
        edr_errors.append(equal(np.column_stack([same,same-other,nul.mean(0)]),g[['total_BC','identity_gain','null_response_mean']].to_numpy(),'Independent EDR reproduction'))
    # Recompute all nested Dc covariates from the bundled selected control counts.
    dc=load_npz(prep/'Dc_counts.npz').toarray();small=np.load(prep/'Dc64_indices.npy');cols=pd.read_csv(prep/'evaluation_genes.csv').source_gene_index.to_numpy(int);cov_errors=[]
    for q,seed in [(1.,0)]+[(q,s) for q in [.5,.25] for s in [20261111,20261112,20261113]]:
        rng=np.random.default_rng(int.from_bytes(hashlib.sha256(f'rpe1_measurement_validation_v1|Dc|q={q:g}|{seed}'.encode()).digest()[:16],'little'))
        x=dc.astype(np.int64) if q==1 else rng.binomial(dc.astype(np.int64),q)
        norm=np.log1p(x.astype(float)*(10000/x.sum(1,dtype=float)[:,None]))
        for n,ids in [(400,np.arange(400)),(64,small)]:cov_errors.append(equal(norm[ids].mean(0)[cols],cov[f'q{q:g}_seed{seed}_n{n}'],'Nested full-gene covariance reproduction'))
    return dict(score_exports=35,safe_numeric_members=numeric_members,target_records=65536,summary_table_rows=len(regenerated),
        condition_summary_max_error=table_error,all_raw_metrics_max_error=raw_error,representative_EDR_conditions=len(selected),
        representative_EDR_max_error=max(edr_errors),covariate_vectors=14,covariate_max_error=max(cov_errors),
        full_raw_download_or_model_refit_performed=False)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT.parent/'CARF-measurement-study-2.0.0-check')
    parser.add_argument('--checksums-only',action='store_true')
    args=parser.parse_args();manifest=checksum_check()
    if args.checksums_only:
        print(f"PASS: {len(manifest['payload_files'])} payload hashes and byte counts")
        return
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    start=datetime.now(timezone.utc).isoformat();checks=numeric_checks(out)
    posthoc=[sys.executable,str(ROOT/'verify_posthoc_review.py'),'--output',str(out/'posthoc')]
    result=subprocess.run(posthoc,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'posthoc_reproduction.log').write_text(result.stdout)
    if result.returncode:raise RuntimeError('Posthoc reconstruction failed; see posthoc_reproduction.log')
    posthoc_report=json.loads((out/'posthoc/summary.json').read_text())
    test=[sys.executable,'-m','unittest','discover','-s',str(ROOT/'source_snapshot/reconstruction/tests'),'-v']
    result=subprocess.run(test,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'unit_tests.log').write_text(result.stdout)
    if result.returncode:raise RuntimeError('Unit tests failed; see unit_tests.log')
    public_test=[sys.executable,'-m','unittest','discover','-s',str(ROOT/'public_tests'),'-v']
    result=subprocess.run(public_test,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'public_provenance_tests.log').write_text(result.stdout)
    if result.returncode:raise RuntimeError('Public provenance tests failed')
    figure=[sys.executable,str(ROOT/'source_snapshot/reconstruction/plot_rpe1_measurement.py'),
        '--evaluation',str(ROOT/'data/rpe1/evaluation'),'--review-summary',str(ROOT/'provenance/rpe1_independent_review/summary.json'),
        '--out',str(out/'figures')]
    result=subprocess.run(figure,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'figure_reproduction.log').write_text(result.stdout)
    if result.returncode:raise RuntimeError('Figure reproduction failed; see figure_reproduction.log')
    expected=ROOT/'figures/rpe1';comparisons=[]
    for path in sorted(expected.rglob('*.csv')):
        generated=out/'figures'/path.relative_to(expected)
        if not generated.exists() or sha(generated)!=sha(path):raise AssertionError(f'Figure source table differs: {path.name}')
        comparisons.append(str(path.relative_to(expected)))
    if not comparisons:raise AssertionError('No figure source tables were checked')
    # Run the unchanged journal exporter against an external portable staging tree.
    stage=out/'journal_input'
    stage.mkdir()
    shutil.copytree(ROOT/'data/rpe1/evaluation',stage/'evaluation')
    shutil.copytree(ROOT/'figures/rpe1',stage/'figures')
    journal_code="from pathlib import Path; import sys; sys.path.insert(0,sys.argv[1]); import export_journal_figures as export; export.BASE=Path(sys.argv[2]); export.OUT=Path(sys.argv[3]); export.main()"
    journal=[sys.executable,'-c',journal_code,str(ROOT/'source_snapshot/reconstruction'),str(stage),str(out/'journal_figures')]
    result=subprocess.run(journal,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'journal_figure_reproduction.log').write_text(result.stdout)
    if result.returncode:raise RuntimeError('Journal figure export failed; see journal_figure_reproduction.log')
    journal_tables=[]
    for path in sorted((ROOT/'figures/journal').glob('*.csv')):
        generated=out/'journal_figures'/path.name
        if not generated.exists() or sha(generated)!=sha(path):raise AssertionError('Journal figure table differs: '+path.name)
        journal_tables.append(path.name)
    if len(journal_tables)!=10:raise AssertionError('Unexpected journal figure source table coverage')
    journal_pngs={path.name:sha(path)==sha(out/'journal_figures'/path.name) for path in sorted((ROOT/'figures/journal').glob('*.png'))}
    runtime=dict(executable=sys.executable,prefix=sys.prefix,base_prefix=sys.base_prefix,python=sys.version,
        platform=platform.platform(),packages={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','h5py','matplotlib']})
    report=dict(status='passed_public_package_reproduction',started_utc=start,finished_utc=datetime.now(timezone.utc).isoformat(),
        release_manifest_sha256=sha(ROOT/'release_manifest.json'),payload_files_checked=len(manifest['payload_files']),
        numeric_checks=checks,posthoc_command=posthoc,posthoc_recomputation=posthoc_report,unit_test_command=test,unit_test_return_code=0,
        unit_tests_passed=int(re.search(r'Ran (\d+) tests', (out/'unit_tests.log').read_text()).group(1)),
        public_provenance_tests_passed=int(re.search(r'Ran (\d+) tests',(out/'public_provenance_tests.log').read_text()).group(1)),figure_command=figure,
        journal_figure_command=journal,journal_source_tables_exactly_reproduced=journal_tables,journal_PNG_binary_equality=journal_pngs,
        figure_source_tables_exactly_reproduced=comparisons,runtime=runtime,
        scope='Checksums; original raw metrics and selected EDR; all 768 posthoc ranking/matching conditions from frozen scores; strict summaries and paired/factorial tables; nested Dc normalization; mathematical tests; original and journal figure rendering. Omitted repeated NPZ grids are not compared by the portable workflow; original full-grid independent audits are retained. No full raw/model refit.')
    (out/'verification_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
