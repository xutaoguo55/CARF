"""Publication figures from audited final RPE1 CSVs only; no analysis mutation.

The portable CLI needs only the listed CSVs plus the unchanged evaluation
manifest and independent-review summary. Unused model/raw/NPZ files are not read.
Ranges depict recorded folds/seeds/permutations, never biological confidence.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE / "results/rpe1_measurement_validation_v1"
INPUTS = ["condition_summary.csv", "fold_summary.csv", "seed_summary.csv",
          "per_target_metrics.csv.gz", "paired_adjustment_effects.csv.gz",
          "operational_criteria.csv"]
MODELS = ["native_scTenifoldKnk", "control_only_Pearson", "heldout_linear", "training_condition_mean"]
NAMES = {"native_scTenifoldKnk": "scTenifoldKnk*", "control_only_Pearson": "Pearson association",
         "heldout_linear": "Held-out linear", "training_condition_mean": "Training mean",
         "observed_A_technical_anchor": "Empirical A anchor"}
FACTORS = ["S0_C0", "S0_Cq", "Sq_C0", "Sq_Cq"]
FTICKS = ["S₀C₀", "S₀Cq", "SqC₀", "SqCq"]
RAW, EDR, SMALL = "#285F8F", "#C65E19", "#137F7C"
FOLD_COLORS = ["#7258A6", "#008C91", "#B54F65"]
FOLD_MARKERS = ["o", "s", "^"]
CAPTIONS = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 9, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": .7, "lines.linewidth": 1.6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 240,
        "axes.titlelocation": "left", "figure.facecolor": "white"})


def tidy(ax, axis="y"):
    ax.grid(axis=axis, color="#E7E7E7", lw=.65, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, color="#777777")


def save(fig, out, stem, caption):
    fig.savefig(out / (stem + ".pdf"), bbox_inches="tight")
    fig.savefig(out / (stem + ".png"), bbox_inches="tight")
    plt.close(fig)
    CAPTIONS[stem] = caption


def source(frame, out, name):
    frame.to_csv(out / (name + ".csv"), index=False, float_format="%.17g")


def baseline(d):
    return d[(d.kind == "baseline") & (d.q == 1) & d.algorithm_seed.isin([0, 42])]


def fig1(data, out):
    cond = baseline(data["condition_summary.csv"])
    fold = baseline(data["fold_summary.csv"])
    cond = cond[cond.n_covariate == 400].copy()
    fold = fold[fold.n_covariate == 400].copy()
    keep = MODELS + ["observed_A_technical_anchor"]
    perm = cond[cond.model.str.startswith("training_label_permutation_")].copy()
    source(cond[cond.model.isin(keep)], out, "Figure_1_baseline_means")
    source(fold[fold.model.isin(keep)], out, "Figure_1_baseline_folds")
    source(perm, out, "Figure_1_permutation_draws")
    order = MODELS + ["permuted_linear", "observed_A_technical_anchor"]
    fig, axs = plt.subplots(1, 3, figsize=(12.6, 5.25), gridspec_kw={"width_ratios": [1.15, 1.2, .66]})
    fig.subplots_adjust(left=.165, right=.985, top=.83, bottom=.23, wspace=.36)
    for j, (metric, title) in enumerate([("total_BC", "A  Total response"), ("identity_gain", "B  Target identity: models")]):
        ax = axs[j]
        for i, model in enumerate(order):
            for adj, color, offset in [("raw", RAW, -.13), ("EDR", EDR, .13)]:
                y = i + offset
                if model == "permuted_linear":
                    vals = perm.loc[perm.adjustment == adj, metric].to_numpy() * 1000
                    ax.plot([vals.min(), vals.max()], [y, y], color=color, alpha=.28, lw=6, solid_capstyle="butt")
                    ax.scatter(vals, np.full(len(vals), y), marker="x", s=26, color=color, lw=1.2, zorder=4)
                elif j == 1 and model == "observed_A_technical_anchor":
                    continue
                else:
                    val = cond.loc[(cond.model == model) & (cond.adjustment == adj), metric].item() * 1000
                    ax.scatter(val, y, color=color, s=52, edgecolors="white", linewidth=.7, zorder=5)
                    f = fold[(fold.model == model) & (fold.adjustment == adj)]
                    for row in f.itertuples():
                        ax.scatter(getattr(row, metric) * 1000, y + (row.fold-1)*.045,
                            marker=FOLD_MARKERS[row.fold], s=19, facecolor="none", edgecolor=color, alpha=.7, zorder=3)
        ax.set_ylim(5.6, -.6)
        ax.set_yticks(range(len(order)), [NAMES.get(x, "Linear label permutations") for x in order] if j == 0 else [""]*len(order))
        ax.set_title(title, fontweight="bold", pad=10)
        ax.set_xlabel(("Top-50 B × C utility" if j == 0 else "Within-fold identity gain") + " (×10⁻³)")
        ax.axvline(0, color="#707070", lw=.85, ls="--")
        tidy(ax, "x")
    axs[0].set_xlim(-3, 100)
    axs[1].set_xlim(-4.1, 9.5)
    axs[1].text(1.8, 5, "Technical anchor →", va="center", fontsize=9, color="#555555")
    ax = axs[2]
    anchor = cond[cond.model == "observed_A_technical_anchor"]
    af = fold[fold.model == "observed_A_technical_anchor"]
    for f in range(3):
        vals = [af.loc[(af.fold == f) & (af.adjustment == a), "identity_gain"].item()*1000 for a in ["raw", "EDR"]]
        ax.plot([0, 1], vals, color="#A6A6A6", alpha=.65, lw=1)
        ax.scatter([0, 1], vals, marker=FOLD_MARKERS[f], facecolors="none", edgecolors="#777777", s=30)
    vals = [anchor.loc[anchor.adjustment == a, "identity_gain"].item()*1000 for a in ["raw", "EDR"]]
    ax.plot([0, 1], vals, color="#333333", lw=1.5)
    ax.scatter([0, 1], vals, c=[RAW, EDR], s=65, edgecolors="white", zorder=4)
    ax.set_xticks([0, 1], ["Raw", "EDR"])
    ax.set_xlim(-.4, 1.4)
    ax.set_ylim(0, 70)
    ax.set_ylabel("Within-fold identity (×10⁻³)")
    ax.set_title("C  Empirical A", fontweight="bold", pad=10)
    tidy(ax)
    handles = [Line2D([], [], color=RAW, marker="o", ls="", label="Raw ranking"),
               Line2D([], [], color=EDR, marker="o", ls="", label="EDR ranking"),
               Line2D([], [], color="#777777", marker="o", mfc="none", ls="", label="Small open symbols: fold means"),
               Line2D([], [], color="#777777", marker="x", lw=4, alpha=.5, label="Five permutation draws and range")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.55, .965), ncol=2, frameon=False)
    fig.text(.165, .075, "128 evaluation targets; q = 1; Dc = 400. Filled dots: target-weighted means. No biological confidence intervals.\n"
             "* Native EDR has 2/128 all-floored baseline targets, retained as ties. Empirical A is a technical anchor, not a prediction model.", fontsize=9, color="#4C4C4C")
    save(fig, out, "Figure_1_baseline", "Baseline RPE1 rankings (q=1, Dc400, native algorithm seed42). A: total utility; B: model identity gain; C: empirical-A identity on its own scale. Filled points average 128 targets; open symbols show three technical target folds. Five label-permutation draws are shown individually with their min–max range, not a confidence interval or permutation test. The native model's two all-floored tasks remain included. Model supervision and feature budgets differ.")


def fig2(data, out):
    d = data["condition_summary.csv"]
    d = d[(d.kind == "factorial") & (d.q == .25) & d.model.isin(MODELS) & d.algorithm_seed.isin([0, 42])].copy()
    seeds = data["seed_summary.csv"]
    seeds = seeds[(seeds.kind == "factorial") & (seeds.q == .25) & seeds.model.isin(MODELS) & seeds.algorithm_seed.isin([0, 42])].copy()
    source(d, out, "Figure_2_q025_condition_means")
    source(seeds, out, "Figure_2_q025_seed_means")
    # Raw values are covariate-independent; collapse only after exact verification.
    keys = ["model", "factor", "adjustment"]
    a = d[(d.adjustment == "raw") & (d.n_covariate == 400)].set_index(keys)
    b = d[(d.adjustment == "raw") & (d.n_covariate == 64)].set_index(keys)
    metrics = ["top50_agreement_to_q1_same_n", "total_BC", "identity_gain"]
    np.testing.assert_array_equal(a[metrics].sort_index(), b[metrics].sort_index())
    fig, axs = plt.subplots(3, 4, figsize=(13, 9.1))
    fig.subplots_adjust(left=.075, right=.99, top=.86, bottom=.17, wspace=.24, hspace=.24)
    styles = [("raw", 400, RAW, "o", "-", "Raw (Dc400 = Dc64)"),
              ("EDR", 400, EDR, "s", "-", "EDR · Dc400"),
              ("EDR", 64, SMALL, "^", "--", "EDR · Dc64")]
    for j, model in enumerate(MODELS):
        for i, metric in enumerate(metrics):
            ax = axs[i, j]
            multiplier = 1 if i == 0 else 1000
            for adj, n, color, marker, line, _ in styles:
                vals = d[(d.model == model) & (d.adjustment == adj) & (d.n_covariate == n)].set_index("factor").loc[FACTORS, metric].to_numpy()*multiplier
                small = seeds[(seeds.model == model) & (seeds.adjustment == adj) & (seeds.n_covariate == n)]
                bounds = small.groupby("factor")[metric].agg(["min", "max"]).loc[FACTORS]*multiplier
                ax.vlines(range(4), bounds["min"], bounds["max"], color=color, lw=1, alpha=.45)
                ax.plot(range(4), vals, color=color, marker=marker, ms=4, ls=line, zorder=3)
            ax.axvspan(-.15, 1.15, color="#F3F3F3", zorder=-2)
            ax.axvspan(1.85, 3.15, color="#F3F3F3", zorder=-2)
            ax.set_xlim(-.2, 3.2)
            ax.set_xticks(range(4), FTICKS if i == 2 else [""]*4)
            if i == 0:
                ax.set_ylim(-.03, 1.04)
                ax.set_title(NAMES[model], fontweight="bold", pad=11)
            elif i == 1:
                ax.set_ylim(0, 61)
            else:
                ax.set_ylim(-1.5, 4.8)
                ax.axhline(0, color="#666666", ls="--", lw=.7)
            if j == 0:
                ax.set_ylabel(["Top-50 agreement\nto q1 at the same Dc size", "Total utility (×10⁻³)", "Within-fold identity\ngain (×10⁻³)"][i])
            tidy(ax)
    handles = [Line2D([], [], color=c, marker=m, ls=l, label=label) for _, _, c, m, l, label in styles]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.55, .97), ncol=3, frameon=False)
    fig.text(.075, .91, "Molecule retention q = 0.25  |  Four arms separate covariate measurement from model-input measurement", fontweight="bold", fontsize=12)
    fig.text(.075, .075, "Shaded pairs change only C: S₀C₀ → S₀Cq and SqC₀ → SqCq. The jump between pairs also changes model input S.\n"
             "Points average three seeds within each target; faint vertical ranges span seed means. * Native EDR: 2 baseline and 5 distinct q=.25 all-floor targets retained.\n"
             "Agreement here is to the q1 baseline; covariate-only operational criteria use pairwise agreement (Figure S1). No biological CIs.", fontsize=9, color="#4C4C4C")
    save(fig, out, "Figure_2_factorial_q025", "Four-arm response at q=0.25 for four model types. The raw line is identical at Dc400 and Dc64; EDR lines use each covariate budget. Top-50 agreement is to the same-budget q1 baseline, not the pairwise covariate criterion. Shaded adjacent pairs differ only in measured C. Values average three seeds within targets; vertical ranges span three seed means. Native all-floor EDR tasks remain included and are flagged in source data. Total response can change independently of target identity.")


def fig_s1(data, out):
    d = data["operational_criteria.csv"]
    d = d[(d.q == .25) & (d.adjustment == "EDR") & d.model.isin(MODELS)].copy()
    d["order"] = d.model.map({x:i for i,x in enumerate(MODELS)})*4 + d.n_covariate.map({400:0,64:2}) + d.effect.map({"C_at_S0":0,"C_at_Sq":1})
    d = d.sort_values("order").reset_index(drop=True)
    source(d.drop(columns="order"), out, "Figure_S1_covariate_operational_criteria")
    fig, axs = plt.subplots(1, 3, figsize=(12.4, 8.2), sharey=True)
    fig.subplots_adjust(left=.28, right=.98, top=.88, bottom=.18, wspace=.35)
    specs = [("agreement_mean", .90, "A  Mean pairwise agreement", (.86,1.005)),
             ("agreement_target_p10", .80, "B  Target 10th percentile", (.77,1.005)),
             ("mean_absolute_total_change_over_A_scale", .10, "C  Mean |utility change| / D_A", (0,.115))]
    for ax,(metric,limit,title,xlim) in zip(axs,specs):
        for i,row in d.iterrows():
            color = EDR if row.n_covariate == 400 else SMALL
            marker = "x" if row.numerically_uninformative_targets else ("s" if row.n_covariate == 400 else "^")
            ax.scatter(row[metric],i,color=color,marker=marker,s=38,zorder=3)
            if metric == "agreement_mean":
                ax.plot([row.seed_agreement_min,row.seed_agreement_max],[i,i],color=color,lw=1,alpha=.6)
        ax.axvline(limit,color="#777777",ls="--",lw=1)
        for k in [3.5,7.5,11.5]:ax.axhline(k,color="#D0D0D0",lw=.6)
        ax.set_xlim(*xlim);ax.set_title(title,fontweight="bold",fontsize=10)
        ax.xaxis.set_major_locator(MaxNLocator(4));tidy(ax,"x")
    labels=[]
    for row in d.itertuples():
        flag = "*" if row.numerically_uninformative_targets else ("†" if row.ranking_seed_threshold_crossing else "")
        label = NAMES[row.model].replace("*","") + f" | Dc{row.n_covariate} | " + ("fixed S₀" if row.effect=="C_at_S0" else "fixed Sq") + flag
        labels.append(label)
    axs[0].set_yticks(range(len(d)),labels);axs[0].set_ylim(len(d)-.5,-.5)
    for ax in axs:ax.set_xlabel("Operational readout")
    fig.text(.28,.95,"Covariate-only contrasts at q = 0.25 (EDR)",fontweight="bold",fontsize=13)
    fig.text(.03,.07,"Dashed lines: prespecified references (mean ≥0.90; p10 ≥0.80; scaled change ≤0.10), not equivalence bounds.\n"
             "Mean-agreement whiskers span three seed means. * Native all-floor EDR tasks: full-scope ranking information is insufficient.\n"
             "† At least one seed crosses the ranking threshold. Training mean / Dc64 / fixed Sq is just below 0.90 at full precision; it is not labeled a pass.",fontsize=9,color="#4C4C4C")
    save(fig,out,"Figure_S1_covariate_criteria","Covariate-only EDR contrasts S0_Cq−S0_C0 and Sq_Cq−Sq_C0 at q=0.25. Mean agreement, target tenth-percentile agreement, and absolute utility change scaled by the fixed empirical-A separation are separate readouts. Thin ranges are seed ranges. Native all-floored tasks preclude interpreting threshold agreement as identifiable full-scope stability. Operational thresholds are not biological equivalence criteria. Classification uses unrounded values.")


def fig3(data,out):
    p=data["paired_adjustment_effects.csv.gz"]
    p=p[(p.model=="heldout_linear") & p.kind.isin(["baseline","factorial"])].copy()
    fields=["identity_gain_raw","identity_gain","identity_gain_EDR_minus_raw"]
    keys=["q","n_covariate","factor","kind","target","fold"]
    target=p.groupby(keys,as_index=False)[fields].mean()
    folds=target.groupby(["q","n_covariate","factor","kind","fold"],as_index=False)[fields].mean()
    means=target.groupby(["q","n_covariate","factor","kind"],as_index=False)[fields].mean()
    source(target,out,"Figure_3_linear_target_pairs_seed_averaged")
    source(folds,out,"Figure_3_linear_fold_pairs")
    source(means,out,"Figure_3_linear_target_weighted_means")
    fig,axs=plt.subplots(2,3,figsize=(12.2,7.2))
    fig.subplots_adjust(left=.08,right=.985,top=.84,bottom=.20,wspace=.29,hspace=.32)
    for i,n in enumerate([400,64]):
        ax=axs[i,0]
        f=folds[(folds.q==1)&(folds.n_covariate==n)]
        for row in f.itertuples():
            ax.plot([0,1],[row.identity_gain_raw*1000,row.identity_gain*1000],color=FOLD_COLORS[row.fold],marker=FOLD_MARKERS[row.fold],ms=5)
        v=means[(means.q==1)&(means.n_covariate==n)].iloc[0]
        ax.plot([0,1],[v.identity_gain_raw*1000,v.identity_gain*1000],color="#222222",marker="D",ms=5,lw=2,ls="--")
        ax.set_xlim(-.15,1.15);ax.set_xticks([0,1],["Raw","EDR"]);ax.set_ylim(-1,5.1)
        ax.set_title(f"{'AD'[i]}  q = 1; Dc{n}",fontweight="bold")
        ax.set_ylabel("Within-fold identity (×10⁻³)");ax.axhline(0,color="#888888",lw=.75,ls=":");tidy(ax)
        for j,q in enumerate([.5,.25],1):
            ax=axs[i,j]
            for fold in range(3):
                vals=folds[(folds.q==q)&(folds.n_covariate==n)&(folds.fold==fold)].set_index("factor").loc[FACTORS,"identity_gain_EDR_minus_raw"].to_numpy()*1000
                ax.plot(range(4),vals,color=FOLD_COLORS[fold],marker=FOLD_MARKERS[fold],ms=4,lw=1.3)
            vals=means[(means.q==q)&(means.n_covariate==n)].set_index("factor").loc[FACTORS,"identity_gain_EDR_minus_raw"].to_numpy()*1000
            ax.plot(range(4),vals,color="#222222",marker="D",ms=4,lw=2,ls="--")
            ax.axhline(0,color="#888888",lw=.75,ls=":");ax.set_ylim(-5.8,2.1)
            ax.set_xticks(range(4),FTICKS);ax.set_xlim(-.2,3.2)
            ax.set_title(f"{'BCEF'[i*2+j-1]}  q = {q:g}; Dc{n}",fontweight="bold")
            ax.set_ylabel("Paired identity change\nEDR − raw (×10⁻³)");tidy(ax)
    handles=[Line2D([],[],color=FOLD_COLORS[f],marker=FOLD_MARKERS[f],label=f"Fold {f+1}") for f in range(3)]
    handles.append(Line2D([],[],color="#222222",marker="D",ls="--",label="Target-weighted mean"))
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.54,.965),ncol=4,frameon=False)
    fig.text(.08,.11,"Left: paired raw/EDR identity at q1. Other panels: adjustment-induced identity change, averaged across three thinning seeds.\n"
             "Fold lines describe technical target partitions, not biological confidence intervals. EDR changes rankings only; signed prediction MSE\n"
             "is unchanged by EDR and is shown separately in Figure S2. Linear predictions have no all-floor EDR cases.",fontsize=9,color="#4C4C4C")
    save(fig,out,"Figure_3_linear_identity","Held-out linear identity by target fold, covariate budget and thinning condition. Left: paired raw and EDR baseline values. Middle/right: EDR−raw identity changes at q=0.5/0.25. Three-seed averages are computed within target before fold and overall means. Black dashed lines weight targets equally, not folds equally. Signed predictions are not altered by the ranking adjustment, so MSE is presented separately.")


def fig_s2(data,out):
    d=data["seed_summary.csv"]
    d=d[d.model.isin(["heldout_linear","training_condition_mean"]) & (d.adjustment=="raw") & (d.n_covariate==400) & (((d.kind=="baseline")&(d.q==1))|((d.kind=="factorial")&(d.factor=="Sq_C0")))].copy()
    source(d,out,"Figure_S2_signed_prediction_seed_metrics")
    fig,axs=plt.subplots(1,2,figsize=(8.6,4.6));fig.subplots_adjust(left=.095,right=.98,top=.80,bottom=.26,wspace=.30)
    for ax,arm in zip(axs,["B","C"]):
        metric="unadjusted_prediction_MSE_"+arm
        for model,color,marker in [("heldout_linear",RAW,"o"),("training_condition_mean",EDR,"s")]:
            g=d[d.model==model].groupby("q")[metric].agg(["mean","min","max"]).loc[[1,.5,.25]]
            ax.plot(range(3),g["mean"],marker=marker,color=color,ms=5,label=NAMES[model])
            ax.vlines(range(3),g["min"],g["max"],color=color,lw=1,alpha=.7)
        ax.set_xticks(range(3),["1","0.5","0.25"]);ax.set_xlabel("Model-input molecule retention q");ax.set_ylabel(f"Signed prediction MSE vs {arm}")
        ax.set_title(f"{'A' if arm=='B' else 'B'}  Held-out response arm {arm}",fontweight="bold");tidy(ax)
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.55,.96),ncol=2,frameon=False)
    fig.text(.095,.10,"These are original signed predictions. EDR is a ranking-only transform and does not change MSE.\n"
             "Ranges span recorded thinning seeds, not confidence intervals; q1 has one fixed fit.",fontsize=9,color="#4C4C4C")
    save(fig,out,"Figure_S2_signed_prediction_MSE","Signed prediction MSE against B and C for the held-out linear and same-fold training-mean predictors. q changes both A training measurement and the Dm reference. This is independent of the EDR ranking transform: EDR rows carry the identical unadjusted MSE. Ranges are min–max of three thinning seeds, not biological confidence intervals.")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation",type=Path,default=BASE/"evaluation")
    parser.add_argument("--review-summary",type=Path,default=BASE/"independent_review/summary.json")
    parser.add_argument("--out",type=Path,default=BASE/"figures")
    args=parser.parse_args()
    evaluation_manifest=args.evaluation/"manifest.json";manifest=json.loads(evaluation_manifest.read_text())
    review=json.loads(args.review_summary.read_text())
    if review.get("status")!="passed_independent_read_only_recomputation" or review.get("evaluation_manifest_sha256")!=sha(evaluation_manifest):
        raise ValueError("Figures require the exact independently reviewed evaluation")
    registered=manifest["output_sha256"];checks={}
    for name in INPUTS:
        path=args.evaluation/name
        if sha(path)!=registered.get(name):raise ValueError(f"Final CSV hash mismatch: {name}")
        checks[name]=sha(path)
    data={name:pd.read_csv(args.evaluation/name) for name in INPUTS}
    if data["per_target_metrics.csv.gz"].target.nunique()!=128:raise ValueError("Unexpected target scope")
    args.out.mkdir(parents=True,exist_ok=True)
    style();fig1(data,args.out);fig2(data,args.out);fig_s1(data,args.out);fig3(data,args.out);fig_s2(data,args.out)
    flags=data["per_target_metrics.csv.gz"]
    flags=flags[(flags.model=="native_scTenifoldKnk") & (flags.adjustment=="EDR") & (flags.n_covariate==400)]
    flags=flags.drop_duplicates(["model","model_q","model_thinning_seed","algorithm_seed","target"])
    source(flags[["model_q","model_thinning_seed","algorithm_seed","target","all_scores_floored","score_exact_constant","score_floor_count"]],args.out,"native_numerical_information_flags")
    (args.out/"Figure_legends.md").write_text("# RPE1 scientific figure legends\n\n"+"\n\n".join(f"## {stem}\n\n{text}" for stem,text in CAPTIONS.items())+"\n")
    outputs={p.name:sha(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name!="manifest.json"}
    (args.out/"manifest.json").write_text(json.dumps({"status":"figures_from_audited_final_CSVs","created_utc":datetime.now(timezone.utc).isoformat(),
        "script_sha256":sha(__file__),"evaluation_manifest_sha256":sha(evaluation_manifest),"independent_review_summary_sha256":sha(args.review_summary),
        "input_sha256":checks,"output_sha256":outputs,"versions":{"numpy":np.__version__,"pandas":pd.__version__,"matplotlib":matplotlib.__version__},
        "data_scope":"No source counts, model NPZs or evaluation NPZs read; numeric plots derive only from registered final CSVs.",
        "interval_policy":"Folds, thinning-seed ranges and label-permutation ranges are descriptive; no biological confidence intervals.",
        "visual_inspection":"PNG inspection must be performed after this script completes; not asserted by script."},indent=2)+"\n")
    print(f"Saved {len(CAPTIONS)} PDF/PNG figures and source CSVs to {args.out}")


if __name__=="__main__":main()
