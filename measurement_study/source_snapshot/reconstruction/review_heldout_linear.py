"""Independent covariance-eigen fit audit and descriptive development figure."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io_utils import HERE, sha256, write_json

SOURCE=HERE/'results/heldout_linear_diagnostic_v1'
OUT=HERE/'results/heldout_linear_diagnostic_v1_review'

def direct_tie_weights(values,k=50):
    threshold=np.sort(values)[-k]
    above=values>threshold; equal=values==threshold
    return above.astype(float)+equal*((k-above.sum())/equal.sum())

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'review_manifest.json').exists(): raise SystemExit('Review already complete')
    source_manifest=sha256(SOURCE/'manifest.json')
    checks=[]
    for ds in ['Adamson','Norman']:
        d=np.load(SOURCE/ds/'A_post_and_Dm_reference.npz',allow_pickle=False)
        signed=np.load(SOURCE/ds/'signed_predictions.npz',allow_pickle=False)
        for fold in range(3):
            f=np.load(SOURCE/ds/'fits'/f'fold_{fold}.npz',allow_pickle=False)
            train=f['training_indices']; test=f['test_indices']; Y=d['A_post'][:,train]
            X=Y-Y.mean(0,keepdims=True)
            eigenvalues,V=np.linalg.eigh(X.T@X)
            G=X@V[:,np.argsort(eigenvalues)[::-1][:10]]
            b=Y.mean(1,keepdims=True); P=G[d['target_gene_indices'][train]]
            W=np.linalg.inv(G.T@G+0.1*np.eye(10))@G.T@(Y-b)@P@np.linalg.inv(P.T@P+0.1*np.eye(10))
            pred=(G@W@G[d['target_gene_indices'][test]].T+b-d['Dm_mean'][:,None])[d['candidate_gene_indices']].T
            err=float(np.max(np.abs(pred-signed['heldout_linear'][test])))
            if err>1e-7: raise AssertionError('Independent eigen implementation disagrees')
            checks.append(dict(dataset=ds,fold=fold,check='covariance_eigen_predictor',max_absolute_error=err))
        # Recompute OLS without using CARF helper and evaluate with hand-built tie weights.
        p=SOURCE/ds/'evaluations/heldout_linear_EDR_q0.25_seed20261031_Cq.npz'
        z=np.load(p,allow_pickle=False)
        D=np.column_stack([np.ones(len(z['genes'])),z['covariate']])
        L=np.log10(np.maximum(np.abs(signed['heldout_linear']),1e-10))
        residual=L-(D@np.linalg.lstsq(D,L.T,rcond=None)[0]).T
        weights=np.vstack([direct_tie_weights(x) for x in residual])
        if not np.array_equal(weights,z['topk_weights']):raise AssertionError('OLS/tie rankings disagree')
        r=np.load(HERE/'results/library_calibration_v2'/f'{ds}_disjoint_library_pools_seed20261024'/'responses.npz',allow_pickle=False)
        gi={str(g):i for i,g in enumerate(r['genes'])}; ti={str(t):i for i,t in enumerate(r['tasks'])}
        rows=[ti[str(t)] for t in z['targets']]; cols=[gi[str(g)] for g in z['genes']]
        product=r['B'][rows][:,cols]*r['C'][rows][:,cols]
        worst=0.
        for fold in range(3):
            idx=np.flatnonzero(z['folds']==fold)
            actual=product[idx]@weights[idx].T/50
            worst=max(worst,float(np.max(abs(actual-z['within_fold_utility'][np.ix_(idx,idx)]))))
        if worst>1e-12:raise AssertionError('Within-fold utility disagrees')
        checks.append(dict(dataset=ds,check='independent_OLS_ties_within_fold_utility',max_absolute_error=worst))
    pd.DataFrame(checks).to_csv(OUT/'independent_checks.csv',index=False)
    s=pd.read_csv(SOURCE/'condition_summary.csv')
    models=['heldout_linear','training_condition_mean','native_scTenifoldKnk','control_only_Pearson']
    labels=['Linear','Training mean','scTenifoldKnk','Pearson']
    colors=['#225ea8','#757575','#e6550d','#238b45']
    fig,ax=plt.subplots(2,3,figsize=(12,7.2),constrained_layout=True)
    for row,ds in enumerate(['Adamson','Norman']):
        d=s[s.dataset==ds]; baseline=d[d.covariate_condition=='C0']
        x=np.arange(4)
        for j,adj in enumerate(['raw','EDR']):
            vals=[baseline[(baseline.model==m)&(baseline.adjustment==adj)].identity_gain_mean.item() for m in models]
            ax[row,0].bar(x+(j-.5)*.36,vals,.36,label=adj,color=['#636363','#66a61e'][j])
        ax[row,0].axhline(0,color='black',lw=.6);ax[row,0].set_xticks(x,labels,rotation=20,ha='right')
        ax[row,0].set_ylabel(f'{ds}\nWithin-fold identity gain')
        for j,arm in enumerate(['B','C']):
            vals=[baseline[(baseline.model==m)&(baseline.adjustment=='raw')][f'signed_prediction_MSE_{arm}_mean'].item() for m in models[:2]]
            ax[row,1].bar(np.arange(2)+(j-.5)*.36,vals,.36,label=f'MSE vs {arm}',color=['#225ea8','#9ecae1'][j])
        ax[row,1].set_xticks([0,1],labels[:2]);ax[row,1].set_ylabel('Signed prediction MSE')
        for m,label,color in zip(models,labels,colors):
            part=d[(d.model==m)&(d.adjustment=='EDR')].sort_values('q',ascending=False)
            ax[row,2].plot(part.q,part.top50_overlap_mean,'o-',label=label,color=color,ms=4)
        ax[row,2].set_xlim(1.03,.22);ax[row,2].set_xticks([1,.5,.25]);ax[row,2].set_ylim(.7,1.015)
        ax[row,2].set_xlabel('Molecule retention in covariate pool')
        ax[row,2].set_ylabel('EDR Top50 inclusion-weight agreement')
        for a in ax[row]:
            a.spines[['top','right']].set_visible(False)
            a.grid(axis='y',alpha=.15);a.set_axisbelow(True)
    ax[0,0].set_title('Target-specific ranking information')
    ax[0,1].set_title('Accuracy and identity are different')
    ax[0,2].set_title('Frozen prediction; changing covariate')
    ax[0,0].legend(frameon=False,fontsize=8);ax[0,1].legend(frameon=False,fontsize=8)
    ax[1,2].legend(frameon=False,fontsize=8,loc='lower left')
    fig.suptitle('Held-out target development audit — descriptive study means',fontsize=14)
    fig.savefig(OUT/'heldout_linear_overview.png',dpi=180)
    fig.savefig(OUT/'heldout_linear_overview.pdf')
    plt.close(fig)
    if sha256(SOURCE/'manifest.json')!=source_manifest:raise AssertionError('Source manifest changed')
    outputs={p.name:sha256(p) for p in OUT.iterdir() if p.is_file()}
    write_json(OUT/'review_manifest.json',dict(status='completed',source_manifest_sha256=source_manifest,
        code_sha256=sha256(Path(__file__)),checks=checks,output_sha256=outputs,
        scope='Independent eigendecomposition fit and selected OLS/within-fold checks; descriptive figure, no new biological inference'))
    print(pd.DataFrame(checks).to_string(index=False))

if __name__=='__main__':main()
