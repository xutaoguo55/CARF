"""Frozen exploratory cell-split calibration, using verified real H5AD inputs.

This evaluates an empirical repeatability positive control, NOT a deployable
prediction model. No cell split is reported as an independent biological repeat.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from carf_core import (rank_metrics, topk_weights, conditional_moments,
                       quantile_strata, residualize, welch_reference)
from io_utils import HERE, input_path, column, gene_symbols, dense_rows, sha256, write_json


def rng_for(seed, name):
    n = int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "little")
    return np.random.default_rng(np.random.SeedSequence([seed, n]))


def coexpression(control, genes, target):
    idx = np.flatnonzero(genes == target)
    if not len(idx):
        return np.zeros(len(genes)), False
    centered = control - control.mean(axis=0)
    t = centered[:, idx[0]]
    denom = np.sqrt((centered * centered).sum(axis=0) * (t * t).sum())
    r = np.divide(t @ centered, denom, out=np.zeros(len(genes)), where=denom > 0)
    return np.abs(r), True


def continuous_excess(w, response, strata):
    expected = 0.0
    for group in np.unique(strata):
        keep = strata == group
        expected += w[keep].sum() * response[keep].mean()
    return float((w @ response - expected) / w.sum())


def evaluate(task, genes, cov, cfg, seed, dataset, wrong_delta):
    rng = rng_for(seed, task["name"] + "_scores")
    valid = genes != task["name"]
    cmean, cvar, cdetect, control = cov
    coexpr, target_found = coexpression(control, genes, task["name"])
    log_score = np.log10(np.maximum(np.abs(task["A"]), 1e-10))
    edr = residualize(log_score, cmean)
    cps = residualize(edr, coexpr)
    joint = residualize(log_score, np.column_stack([cmean, coexpr]))
    methods = {
        "independent_repeatability_positive": np.abs(task["A"]),
        "legacy_EDR_positive": edr,
        "legacy_two_stage_CPS_positive": cps,
        "joint_CPS_positive": joint,
        "wrong_perturbation": np.abs(wrong_delta),
        "control_mean_only": cmean,
        "control_variance_only": cvar,
        "control_detection_only": cdetect,
        "control_target_coexpression_only": coexpr,
        "random_ranking": rng.random(len(genes)),
        "constant_zero_delta": np.zeros(len(genes)),
    }
    if task["type"] == "control_split_null":
        # Another null partition may reuse this task's evaluation cells.
        # It is not an independent mismatched prediction and must not be scored.
        methods.pop("wrong_perturbation")
    eligible = np.flatnonzero(valid)
    scopes = {"all_curated_genes": eligible}
    # Ties in defining a scope are resolved by stable gene order, never labels.
    expr_order = np.lexsort((genes[eligible], -cmean[eligible]))
    scopes["control_expression_upper_half"] = eligible[expr_order[:len(eligible) // 2]]
    scopes["random_half"] = np.sort(rng.choice(eligible, len(eligible) // 2, replace=False))
    M = int(task["supported"][valid].sum())
    rows = []
    for scope_name, idx in scopes.items():
        expr_strata = quantile_strata(cmean[idx], cfg["covariate_quantile_bins"])
        joint_strata = quantile_strata(np.column_stack([cmean[idx], cvar[idx]]), cfg["covariate_quantile_bins"])
        sizes = np.bincount(joint_strata)
        y = task["supported"][idx]
        b_abs = np.abs(task["B"][idx])
        for method, values in methods.items():
            # Refit residualization inside each candidate scope, matching its estimand.
            if method == "legacy_EDR_positive":
                scores = residualize(log_score[idx], cmean[idx])
            elif method == "legacy_two_stage_CPS_positive":
                scores = residualize(residualize(log_score[idx], cmean[idx]), coexpr[idx])
            elif method == "joint_CPS_positive":
                scores = residualize(log_score[idx], np.column_stack([cmean[idx], coexpr[idx]]))
            else:
                scores = values[idx]
            for k in [cfg["primary_k"]] + cfg["secondary_k"]:
                m = rank_metrics(scores, y, k, M)
                w = topk_weights(scores, k)
                global_m = conditional_moments(w, y)
                expr_m = conditional_moments(w, y, expr_strata)
                joint_m = conditional_moments(w, y, joint_strata)
                m.update(
                    dataset=dataset, seed=seed, perturbation=task["name"], task_type=task["type"],
                    method=method, scope=scope_name, cells_per_arm=task["n"],
                    global_excess_hits=global_m["excess_hits"],
                    expr_excess_hits=expr_m["excess_hits"],
                    expr_standardized_excess=expr_m["standardized_excess"],
                    joint_excess_hits=joint_m["excess_hits"],
                    joint_standardized_excess=joint_m["standardized_excess"],
                    joint_information=joint_m["within_stratum_information"],
                    joint_min_stratum_size=int(sizes.min()),
                    joint_small_stratum_gene_fraction=float(sizes[sizes < cfg["small_stratum_warning_below"]].sum() / len(idx)),
                    mean_abs_response_B=float(w @ b_abs / k),
                    joint_excess_abs_response_B=continuous_excess(w, b_abs, joint_strata),
                    target_present_for_coexpression=target_found,
                )
                rows.append(m)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["Adamson", "Norman"], required=True)
    parser.add_argument("--seed", type=int, default=20261021)
    args = parser.parse_args()
    cfg_path = HERE / "config/pilot_v1.json"
    cfg = json.loads(cfg_path.read_text())
    if args.seed not in [cfg["primary_seed"]] + cfg["sensitivity_seeds"]:
        raise ValueError("Seed not in frozen pilot specification")
    verification = HERE / "provenance" / (args.dataset + "_file_verification.json")
    receipt = json.loads(verification.read_text())
    if not receipt["public_bytes_match"]:
        raise ValueError("Curated file has not passed public-byte verification")
    source = input_path(args.dataset)
    if sha256(source) != receipt["local_sha256"]:
        raise ValueError("Local file has changed since provenance verification")
    out = HERE / "results/pilot_v1" / f"{args.dataset}_seed{args.seed}"
    if (out / "manifest.json").exists():
        raise SystemExit("Completed result already exists; do not silently overwrite")
    out.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    code_hashes = {p.name: sha256(p) for p in [Path(__file__), HERE / "carf_core.py", HERE / "io_utils.py"]}
    with h5py.File(source, "r") as f:
        genes = gene_symbols(f)
        if len(np.unique(genes)) != len(genes):
            raise ValueError("Duplicate gene symbols require an explicit aggregation policy")
        labels = column(f["obs"], "perturbation")
        barcodes = column(f["obs"], f["obs"].attrs["_index"])
        controls = np.flatnonzero(labels == cfg["control_label_exact"])
        rng = rng_for(args.seed, args.dataset + "_controls")
        pools = np.array_split(rng.permutation(controls), 4)
        for i in range(4):
            for j in range(i):
                assert not np.intersect1d(pools[i], pools[j]).size
        D = dense_rows(f, pools[3])
        cov = D.mean(axis=0), D.var(axis=0, ddof=1), (D > 0).mean(axis=0), D
        pd.DataFrame({"gene_symbol": genes, "control_mean_D": cov[0],
                      "control_variance_D": cov[1], "control_detection_D": cov[2]}).to_csv(out / "covariates.csv.gz", index=False)
        counts = pd.Series(labels).value_counts().sort_index()
        tasks, split_rows, inventory = [], [], []
        for condition, n_cells in counts.items():
            reason = ("control" if condition == "control" else "combination" if "+" in condition
                      else "fewer_than_180_cells" if n_cells < cfg["min_perturbation_cells"] else "included")
            inventory.append({"condition": condition, "cells": int(n_cells), "decision": reason})
            if reason != "included":
                continue
            rng = rng_for(args.seed, args.dataset + "_" + condition)
            pert = rng.permutation(np.flatnonzero(labels == condition))
            n = min(len(pert) // 3, cfg["max_cells_per_arm"], min(len(p) for p in pools[:3]))
            pa, pb, pc = pert[:n], pert[n:2*n], pert[2*n:3*n]
            ca, cb, cc = [rng.choice(p, n, replace=False) for p in pools[:3]]
            groups = [pa, pb, pc, ca, cb, cc, pools[3]]
            flat = np.concatenate(groups)
            assert len(np.unique(flat)) == len(flat), "overlapping prediction/evaluation cells"
            A = dense_rows(f, pa).mean(axis=0) - dense_rows(f, ca).mean(axis=0)
            B = dense_rows(f, pb).mean(axis=0) - dense_rows(f, cb).mean(axis=0)
            C, q, supported = welch_reference(dense_rows(f, pc), dense_rows(f, cc),
                                               cfg["reference_fdr"], cfg["reference_min_abs_log_normalized_delta"])
            tasks.append({"name": condition, "n": n, "A": A, "B": B, "C": C, "q": q,
                          "supported": supported, "type": "real_single_gene_perturbation"})
            for arm, group in zip(["pert_A", "pert_B", "pert_C", "ctrl_A", "ctrl_B", "ctrl_C"], groups[:6]):
                split_rows.extend({"task": condition, "arm": arm, "cell_index": int(i), "barcode": barcodes[i]} for i in group)
        # Null cells drawn from controls outside D, all six groups disjoint.
        null_pool = np.concatenate(pools[:3])
        for i in range(cfg["real_null_partitions"]):
            rng = rng_for(args.seed, args.dataset + "_null_" + str(i))
            n = cfg["real_null_cells_per_arm"]
            groups = np.array_split(rng.choice(null_pool, 6*n, replace=False), 6)
            assert len(np.unique(np.concatenate(groups))) == 6*n
            A = dense_rows(f, groups[0]).mean(axis=0) - dense_rows(f, groups[1]).mean(axis=0)
            B = dense_rows(f, groups[2]).mean(axis=0) - dense_rows(f, groups[3]).mean(axis=0)
            C, q, supported = welch_reference(dense_rows(f, groups[4]), dense_rows(f, groups[5]),
                                               cfg["reference_fdr"], cfg["reference_min_abs_log_normalized_delta"])
            name = f"control_vs_control_{i:03d}"
            tasks.append({"name": name, "n": n, "A": A, "B": B, "C": C, "q": q,
                          "supported": supported, "type": "control_split_null"})
            for arm, group in zip(["pert_A", "ctrl_A", "pert_B", "ctrl_B", "pert_C", "ctrl_C"], groups):
                split_rows.extend({"task": name, "arm": arm, "cell_index": int(j), "barcode": barcodes[j]} for j in group)
    split_rows.extend({"task": "all_tasks", "arm": "covariates_D", "cell_index": int(i), "barcode": barcodes[i]} for i in pools[3])
    pd.DataFrame(split_rows).to_csv(out / "cell_splits.csv.gz", index=False)
    pd.DataFrame(inventory).to_csv(out / "task_inclusion.csv", index=False)
    np.savez_compressed(out / "responses_and_references.npz",
                        genes=genes, tasks=np.array([t["name"] for t in tasks]),
                        **{key: np.stack([t[key] for t in tasks]) for key in ["A", "B", "C", "q", "supported"]})
    print(f"{args.dataset}: prepared {sum(t['type'].startswith('real') for t in tasks)} real perturbations and {cfg['real_null_partitions']} null splits", flush=True)
    real = [t for t in tasks if t["type"].startswith("real")]
    wrong = {t["name"]: real[(i+1) % len(real)]["A"] for i, t in enumerate(real)}
    null = [t for t in tasks if t["type"] == "control_split_null"]
    wrong.update({t["name"]: null[(i+1) % len(null)]["A"] for i, t in enumerate(null)})
    rows, quality = [], []
    for i, task in enumerate(tasks):
        rows.extend(evaluate(task, genes, cov, cfg, args.seed, args.dataset, wrong[task["name"]]))
        valid = genes != task["name"]
        quality.append({
            "task": task["name"], "task_type": task["type"], "cells_per_arm": task["n"],
            "supported_genes_C": int(task["supported"][valid].sum()),
            "pearson_delta_A_B_independent_controls": float(np.corrcoef(task["A"][valid], task["B"][valid])[0, 1]),
            "pearson_abs_delta_A_B_independent_controls": float(np.corrcoef(np.abs(task["A"][valid]), np.abs(task["B"][valid]))[0, 1]),
        })
        if (i+1) % 20 == 0:
            print(f"{args.dataset}: scored {i+1}/{len(tasks)} tasks", flush=True)
    pd.DataFrame(rows).to_csv(out / "metrics.csv.gz", index=False)
    pd.DataFrame(quality).to_csv(out / "reference_quality.csv", index=False)
    write_json(out / "manifest.json", {
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": cfg, "protocol_sha256": sha256(cfg_path), "code_sha256": code_hashes,
        "source_sha256": receipt["local_sha256"], "file_verification_receipt_sha256": sha256(verification),
        "dataset": args.dataset, "seed": args.seed, "real_tasks": len(real), "null_tasks": len(null),
        "n_metric_rows": len(rows), "claim_status": "exploratory_technical_calibration_only",
        "limitations": [
            "No biological replicate inference; random cells may share guides, batches and experimental history.",
            "D and control pools reused across tasks; tasks and split seeds are not independent.",
            "Upstream selected HVGs on full dataset; confirmation requires independent scope construction.",
            "Conditional moments assume within-stratum exchangeability and are diagnostics, not calibrated biological p-values.",
            "Reference non-members are unsupported/unknown, not validated biological negatives.",
            "Repeatability positive control uses treated cells and is not a held-out perturbation prediction model.",
        ],
        "outputs": {p.name: sha256(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"},
    })
    print(f"Completed {out}", flush=True)


if __name__ == "__main__":
    main()
