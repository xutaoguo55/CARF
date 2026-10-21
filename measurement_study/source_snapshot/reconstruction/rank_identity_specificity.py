"""Post-hoc, common-scope task-identity audit for unsigned gene rankings.

This is a diagnostic of empirical split rankings, not a held-out model
benchmark. Every query task is scored with its independent B x C response
against rankings from every source task, after excluding all directly targeted
genes from one common candidate universe.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from carf_core import residualize, topk_weights
from io_utils import HERE, ROOT, input_path, dense_rows, sha256, write_json
from run_pilot import coexpression


KS = 50
METHODS = ["repeatability_positive", "legacy_EDR", "legacy_joint_CPS",
           "control_mean_only", "control_variance_only",
           "control_target_coexpression_only", "random_ranking"]


def task_score_matrix(genes, tasks, A, covariates, control_D, seed):
    """Return method -> (task x common-gene) score array."""
    targets = set(map(str, tasks))
    keep = np.array([str(g) not in targets for g in genes], dtype=bool)
    if keep.sum() < KS:
        raise ValueError("Common target-excluded scope is smaller than primary k")
    common_genes = genes[keep]
    cmean = covariates.mean_D.to_numpy(float)[keep]
    cvar = covariates.variance_D.to_numpy(float)[keep]
    out = {name: [] for name in METHODS}
    for i, task in enumerate(tasks):
        a = np.asarray(A[i, keep], float)
        coexpr, _ = coexpression(control_D, genes, str(task))
        coexpr = coexpr[keep]
        logabs = np.log10(np.maximum(np.abs(a), 1e-10))
        rng = np.random.default_rng(np.random.SeedSequence([
            int(seed), int.from_bytes(hashlib.sha256(str(task).encode()).digest()[:4], "little")
        ]))
        out["repeatability_positive"].append(np.abs(a))
        out["legacy_EDR"].append(residualize(logabs, cmean))
        out["legacy_joint_CPS"].append(
            residualize(logabs, np.column_stack([cmean, coexpr]))
        )
        out["control_mean_only"].append(cmean)
        out["control_variance_only"].append(cvar)
        out["control_target_coexpression_only"].append(coexpr)
        out["random_ranking"].append(rng.random(len(common_genes)))
    return common_genes, {key: np.vstack(value) for key, value in out.items()}


def evaluate_cross_task(scores, B, C, k=KS):
    """Evaluate every source-task ranking on every query-task B x C response."""
    n_tasks, n_genes = scores.shape
    if B.shape != scores.shape or C.shape != scores.shape:
        raise ValueError("score and response matrices must have the same shape")
    if k > n_genes:
        raise ValueError("k exceeds the common candidate scope")
    # utility[query, source] uses the same query response and same common scope.
    utility = np.empty((n_tasks, n_tasks), dtype=float)
    for source in range(n_tasks):
        w = topk_weights(scores[source], k)
        for query in range(n_tasks):
            utility[query, source] = float(w @ (B[query] * C[query]) / k)
    rows = []
    for query in range(n_tasks):
        offdiag = np.delete(utility[query], query)
        rows.append({
            "query_index": query,
            "source_task_count": n_tasks,
            "same_task_utility": float(utility[query, query]),
            "all_other_task_mean_utility": float(offdiag.mean()),
            "all_other_task_median_utility": float(np.median(offdiag)),
            "same_minus_other_mean": float(utility[query, query] - offdiag.mean()),
            "same_exceeds_other_mean": bool(utility[query, query] > offdiag.mean()),
            "same_exceeds_other_median": bool(utility[query, query] > np.median(offdiag)),
        })
    return utility, pd.DataFrame(rows)


def run_empirical():
    runs = sorted((HERE / "results/library_calibration_v2").glob("*/manifest.json"))
    if len(runs) != 12:
        raise ValueError(f"Expected 12 frozen runs, found {len(runs)}")
    outdir = HERE / "results/rank_identity_specificity_v1"
    outdir.mkdir(parents=True, exist_ok=True)
    run_summaries, task_rows, pair_rows, inputs = [], [], [], {}
    for manifest_path in runs:
        folder = manifest_path.parent
        manifest = json.loads(manifest_path.read_text())
        response_path = folder / "responses.npz"
        if sha256(response_path) != manifest["outputs"]["responses.npz"]:
            raise ValueError(f"response hash mismatch: {response_path}")
        for path in [manifest_path, response_path]:
            inputs[str(path.relative_to(HERE))] = sha256(path)
        z = np.load(response_path, allow_pickle=False)
        genes = z["genes"].astype(str)
        tasks = z["tasks"].astype(str)
        real = np.array([not t.startswith("control_vs_control_") for t in tasks])
        cov = pd.read_csv(folder / "covariates.csv.gz")
        splits = pd.read_csv(folder / "cell_splits.csv.gz")
        d_idx = splits.loc[splits.arm == "covariates_D", "cell_index"].to_numpy(int)
        with h5py.File(input_path(manifest["dataset"]), "r") as f:
            D = dense_rows(f, d_idx)
        genes_common, all_scores = task_score_matrix(
            genes, tasks, z["A"], cov, D, manifest["seed"]
        )
        if len(genes_common) != len(set(genes_common)):
            raise ValueError("common candidate universe has duplicated gene symbols")
        for task_type, indices in [
            ("real_single_gene_perturbation", np.flatnonzero(real)),
            ("control_split_null", np.flatnonzero(~real)),
        ]:
            if len(indices) < 2:
                continue
            for method, score_matrix in all_scores.items():
                _, per_query = evaluate_cross_task(
                    score_matrix[indices], z["B"][indices][:, np.isin(genes, genes_common)],
                    z["C"][indices][:, np.isin(genes, genes_common)], KS
                )
                # Keep indices and names explicit so every row can be traced.
                per_query["query_task"] = tasks[indices]
                per_query["dataset"] = manifest["dataset"]
                per_query["mode"] = manifest["mode"]
                per_query["seed"] = manifest["seed"]
                per_query["task_type"] = task_type
                per_query["method"] = method
                task_rows.append(per_query)
                utility, _ = evaluate_cross_task(
                    score_matrix[indices], z["B"][indices][:, np.isin(genes, genes_common)],
                    z["C"][indices][:, np.isin(genes, genes_common)], KS
                )
                for qi, query in enumerate(indices):
                    for si, source in enumerate(indices):
                        pair_rows.append({
                            "dataset": manifest["dataset"], "mode": manifest["mode"],
                            "seed": manifest["seed"], "task_type": task_type,
                            "method": method, "query_task": tasks[query],
                            "source_task": tasks[source], "is_identity_match": bool(query == source),
                            "mean_BC_product": float(utility[qi, si]),
                            "k": KS, "common_scope_size": int(len(genes_common)),
                        })
                summary = task_rows[-1]
                run_summaries.append({
                    "dataset": manifest["dataset"], "mode": manifest["mode"],
                    "seed": manifest["seed"], "task_type": task_type, "method": method,
                    "n_query_tasks": len(indices), "k": KS,
                    "common_scope_size": int(len(genes_common)),
                    "mean_same_task_utility": float(summary.same_task_utility.mean()),
                    "mean_all_other_task_utility": float(summary.all_other_task_mean_utility.mean()),
                    "mean_identity_gain": float(summary.same_minus_other_mean.mean()),
                    "median_identity_gain": float(summary.same_minus_other_mean.median()),
                    "fraction_queries_same_above_other_mean": float(summary.same_exceeds_other_mean.mean()),
                    "fraction_queries_same_above_other_median": float(summary.same_exceeds_other_median.mean()),
                })
        print(f"Audited {folder.name}: common genes={len(genes_common)}", flush=True)

    per_run = pd.DataFrame(run_summaries)
    per_task = pd.concat(task_rows, ignore_index=True)
    per_pair = pd.DataFrame(pair_rows)
    per_run.to_csv(outdir / "per_run_summary.csv", index=False)
    per_task.to_csv(outdir / "per_query_task.csv.gz", index=False)
    per_pair.to_csv(outdir / "cross_task_utilities.csv.gz", index=False)
    grouped = (per_run.groupby(["dataset", "mode", "task_type", "method"], as_index=False)
               .agg(n_runs=("seed", "nunique"), n_tasks=("n_query_tasks", "mean"),
                    common_scope_size=("common_scope_size", "mean"),
                    mean_same_task_utility=("mean_same_task_utility", "mean"),
                    mean_all_other_task_utility=("mean_all_other_task_utility", "mean"),
                    mean_identity_gain=("mean_identity_gain", "mean"),
                    median_identity_gain=("median_identity_gain", "median"),
                    fraction_queries_same_above_other_mean=("fraction_queries_same_above_other_mean", "mean"),
                    fraction_queries_same_above_other_median=("fraction_queries_same_above_other_median", "mean")))
    grouped.to_csv(outdir / "grouped_summary.csv", index=False)
    declaration = {
        "status": "posthoc_exploratory_common_scope_identity_audit",
        "declared_utc": datetime.now(timezone.utc).isoformat(),
        "primary_estimator": "Top50 mean(delta_B * delta_C) from an unsigned ranking, with task identity gain equal to same-task utility minus mean utility across all other task rankings",
        "common_scope": "Same retained gene-symbol universe for all queries and source rankings; excludes every real single-gene intervention target in the dataset",
        "ranking_source": "Empirical split-A perturbation response and control-derived baselines; no deployed prediction model is represented",
        "positive_control": "Same-task empirical absolute response ranking in split A; it is a technical repeatability reference, not model prediction",
        "negative_control": "Control-vs-control pseudo-task queries evaluated with all other control pseudo-task rankings; technical shared-control partitions",
        "inference": "Descriptive only; tasks, rankings, cells, controls, seeds, and pairings are dependent; no biological p-values or confirmatory intervals",
        "k": KS,
        "inputs": inputs,
    }
    write_json(outdir / "analysis_declaration.json", declaration)
    write_json(outdir / "manifest.json", {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "code_sha256": {"carf_core.py": sha256(HERE / "carf_core.py"),
                        "io_utils.py": sha256(HERE / "io_utils.py"),
                        "run_pilot.py": sha256(HERE / "run_pilot.py")},
        "inputs": inputs,
        "outputs": {p.name: sha256(p) for p in outdir.iterdir()
                    if p.is_file() and p.name != "manifest.json"},
        "limitations": [
            "Post-hoc analysis of existing technical split runs.",
            "Empirical repeatability scores are not model predictions.",
            "Task identity gain is descriptive and does not make shared biology artifactual.",
        ],
    })
    report = [
        "# CARF 共同候选范围下的扰动身份特异性复核（事后探索）", "",
        "本复核使用既有 v2 真实 K562 扰动任务，重新构造严格相同的候选基因范围：所有在数据中被直接扰动的基因均从每个查询与排名中剔除。每个查询任务用独立 B、C 细胞组的 signed delta 乘积评价；每一种排名都与同一数据集、同一拆分设计和同一任务类型的所有排名交叉比较。", "",
        "主要量为 Top50 平均 delta_B×delta_C；身份增量为同任务排名效用减去所有其他任务排名的平均效用。正值表示该经验排名比错配排名更贴合本任务可重复响应。结果只反映技术拆分；来源排名是 split A 的经验响应，不是模型预测。", "",
        "| 数据 | 拆分 | 查询类型 | 排名方法 | 任务数 | 共同候选基因数 | 同任务效用 | 错配效用 | 身份增量 | 同任务高于错配均值比例 |", "|---|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in grouped.itertuples():
        report.append(
            f"| {row.dataset} | {'随机细胞' if row.mode.startswith('random') else '建库组分离'} | "
            f"{'真实扰动' if row.task_type == 'real_single_gene_perturbation' else '对照伪扰动'} | "
            f"{row.method} | {row.n_tasks:.0f} | {row.common_scope_size:.0f} | "
            f"{row.mean_same_task_utility:.5f} | {row.mean_all_other_task_utility:.5f} | "
            f"{row.mean_identity_gain:.5f} | {row.fraction_queries_same_above_other_mean:.2f} |"
        )
    report += ["", "所有逐任务、逐配对分数及哈希清单位于 `reconstruction/results/rank_identity_specificity_v1/`。本复核增加了共同范围和全部错配比较，但不能替代独立生物重复、外部细胞背景或真实模型的留出验证，也不能单独构成新方法贡献。"]
    report_path = ROOT / "review/CARF_共同范围身份特异性复核_20261021.md"
    report_path.write_text("\n".join(report) + "\n")
    manifest = json.loads((outdir / "manifest.json").read_text())
    manifest["outputs"]["report_sha256"] = sha256(report_path)
    write_json(outdir / "manifest.json", manifest)
    print(grouped.to_string(index=False), flush=True)


if __name__ == "__main__":
    run_empirical()
