"""Post-pilot verification, factorial contrasts and a descriptive scientific figure.

This review leaves the declared primary outputs unchanged. Numerical-information
flags are reported rather than silently dropping targets after viewing results.
"""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from io_utils import HERE, sha256, write_json

SOURCE=HERE/"results/measurement_diagnostic_pilot_v1"
OUT=HERE/"results/measurement_diagnostic_pilot_v1_review"
METRICS=["top50_overlap","total_BC","identity_gain","score_expression_rho","logscore_expression_r2"]


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((SOURCE/"manifest.json").read_text())
    for name,expected in manifest["output_sha256"].items():
        if sha256(SOURCE/name)!=expected: raise ValueError(f"Changed primary output: {name}")
    datasets=manifest["protocol"]["datasets"]
    frames=[]; flags=[]; sensitivity=[]; invariants=[]
    for dataset in datasets:
        folder=SOURCE/dataset
        d=pd.read_csv(folder/"per_target_all_conditions.csv.gz")
        frames.append(d)
        base=np.load(folder/"evaluations/baseline_q1_thin0_model42_S0_C0_raw_distance.npz")
        baseline_floor=(base["scores"]<=1e-10).all(axis=1)
        for p in sorted((folder/"evaluations").glob("*.npz")):
            z=np.load(p); scores=z["scores"]
            flags.append(dict(dataset=dataset,file=p.name,targets=len(scores),
                exact_constant_rows=int((np.ptp(scores,axis=1)==0).sum()),
                raw_all_at_or_below_floor_rows=int((scores<=1e-10).all(axis=1).sum()) if p.name.endswith("raw_distance.npz") else None))
            if "_S0_Cq_raw_distance" in p.name:
                ok=all(np.array_equal(z[key],base[key]) for key in ["scores","topk_weights","utility_query_source","null_utility"])
                if not ok: raise AssertionError("Covariate-only intervention changed raw evaluation")
                invariants.append(dict(dataset=dataset,file=p.name,raw_evaluation_exactly_unchanged=True))
            # Explicit post-hoc sensitivity: remove baseline all-floored tasks
            # from both the query and source axes, preserving candidate genes.
            keep=~baseline_floor
            u=z["utility_query_source"][keep][:,keep]
            diag=np.diag(u); other=(u.sum(axis=1)-diag)/(len(diag)-1)
            sensitivity.append(dict(dataset=dataset,file=p.name,n_targets=len(diag),
                excluded_baseline_below_floor_targets=";".join(base["targets"][baseline_floor]),
                total_BC_mean=float(diag.mean()),identity_gain_mean=float((diag-other).mean())))
    pd.DataFrame(flags).to_csv(OUT/"numerical_information_flags.csv",index=False)
    pd.DataFrame(sensitivity).to_csv(OUT/"posthoc_baseline_floor_sensitivity.csv",index=False)
    pd.DataFrame(invariants).to_csv(OUT/"fixed_score_invariants.csv",index=False)
    df=pd.concat(frames,ignore_index=True)
    fact=df[df.kind=="factorial"]
    keys=["dataset","q","thinning_seed","model_seed","method","target"]
    effects=[]
    for metric in METRICS:
        p=fact.pivot(index=keys,columns="factor",values=metric)
        e=p.index.to_frame(index=False)
        e["metric"]=metric
        e["covariate_change"]=(p.S0_Cq-p.S0_C0).to_numpy()
        e["score_change"]=(p.Sq_C0-p.S0_C0).to_numpy()
        e["interaction"]=(p.Sq_Cq-p.Sq_C0-p.S0_Cq+p.S0_C0).to_numpy()
        effects.append(e)
    effect=pd.concat(effects,ignore_index=True)
    effect.to_csv(OUT/"factorial_effects_per_target_seed.csv.gz",index=False)
    group=["dataset","q","method","metric","target"]
    vals=["covariate_change","score_change","interaction"]
    avg=effect.groupby(group)[vals].mean().reset_index()
    avg.groupby(group[:-1])[vals].mean().reset_index().to_csv(OUT/"factorial_effects_summary.csv",index=False)

    summary=pd.read_csv(SOURCE/"condition_summary.csv")
    fig,axes=plt.subplots(2,3,figsize=(13,7.4),constrained_layout=True)
    colors={"S0_Cq":"#20856e","Sq_C0":"#a547b5","Sq_Cq":"#cf6a25"}
    qvalues=[1,.5,.25]
    for ri,dataset in enumerate(datasets):
        sub=summary[summary.dataset==dataset]
        for factor in colors:
            rows=sub[(sub.kind=="factorial")&(sub.factor==factor)&(sub.method=="EDR")].set_index("q")
            vals=[1]+[rows.loc[q,"top50_overlap_mean"] for q in qvalues[1:]]
            axes[ri,0].plot(qvalues,vals,"o-",color=colors[factor],label={"S0_Cq":"Covariate only","Sq_C0":"Model only","Sq_Cq":"Joint"}[factor])
        alg=sub[(sub.kind=="algorithm_seed_control")&(sub.method=="EDR")].top50_overlap_mean.iloc[0]
        axes[ri,0].axhline(alg,color="#777777",ls="--",label="Mean alternate-seed overlap")
        for method,color in [("raw_distance","#296ca0"),("EDR","#cf6a25")]:
            b=sub[(sub.kind=="baseline")&(sub.method==method)].iloc[0]
            rows=sub[(sub.kind=="factorial")&(sub.factor=="Sq_Cq")&(sub.method==method)].set_index("q")
            for col,metric in [(1,"total_BC_mean"),(2,"identity_gain_mean")]:
                vals=[b[metric]]+[rows.loc[q,metric] for q in qvalues[1:]]
                axes[ri,col].plot(qvalues,vals,"o-",color=color,label="Raw distance" if method=="raw_distance" else "EDR")
        for col in range(3):
            ax=axes[ri,col];ax.set_xlim(1.04,.21);ax.set_xticks(qvalues,["100%","50%","25%"])
            ax.grid(alpha=.2);ax.set_xlabel("Retained molecule probability")
            ax.set_title(dataset+" | "+["EDR selection overlap","Joint: total response","Joint: identity increment"][col],fontsize=11)
            if col==0:ax.set_ylim(-.03,1.04)
            if col==2:ax.axhline(0,color="black",lw=.7);ax.ticklabel_format(axis="y",style="sci",scilimits=(-3,-3))
            ax.legend(fontsize=7,loc="best")
    axes[0,0].set_ylabel("Mean overlap of inclusion weights")
    axes[1,0].set_ylabel("Mean overlap of inclusion weights")
    axes[0,1].set_ylabel("Mean held-out B x C")
    axes[1,1].set_ylabel("Mean held-out B x C")
    axes[0,2].set_ylabel("Same target minus other targets")
    axes[1,2].set_ylabel("Same target minus other targets")
    fig.suptitle("Measurement-diagnostic feasibility pilot | descriptive, one model, two K562 studies\nFixed outcomes and gene scope; seeds averaged within targets; no biological inference",fontsize=12)
    fig.savefig(OUT/"measurement_pilot_overview.png",dpi=180)
    fig.savefig(OUT/"measurement_pilot_overview.pdf")
    plt.close(fig)
    write_json(OUT/"review_manifest.json",dict(status="verified_descriptive_postpilot_review",
        primary_manifest_sha256=sha256(SOURCE/"manifest.json"),verified_primary_files=len(manifest["output_sha256"]),
        fixed_score_conditions_checked=len(invariants),script_sha256=sha256(Path(__file__)),
        notes=["Baseline-floor task removal is post-hoc sensitivity, not a change to primary analysis.",
               "Overlap measures agreement of fractional inclusion weights; all-tied scores do not imply an identifiable stable top50.",
               "Causal language applies only to the randomized computational operation conditional on observed data, not biological confounding."],
        outputs={p.name:sha256(p) for p in OUT.iterdir() if p.is_file() and p.name!="review_manifest.json"}))
    print(f"Verified {len(manifest['output_sha256'])} primary files and {len(invariants)} exact raw-invariance conditions; saved figure and contrasts.")


if __name__=="__main__": main()
