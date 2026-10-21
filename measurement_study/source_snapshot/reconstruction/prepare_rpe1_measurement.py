"""Prepare the frozen RPE1 measurement protocol, without fitting or thinning.

The ordered gates are source verification, metadata-only cell/fold freezing,
Dm-only gene/evaluation-target freezing, and then response construction.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz

from verify_rpe1_contents import read_column, source_gate, file_hashes

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "config/rpe1_measurement_validation_v1.json"
OUT = HERE / "results/rpe1_measurement_validation_v1/preparation"
EXPECTED_PROTOCOL_SHA256 = "ed682799610d7c32338b991f437c86dff6f873d2efe2db556f484b44be473913"
PROTOCOL = "rpe1_measurement_validation_v1"


def sha(path):
    return file_hashes(Path(path))["sha256"]


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def rng_for(description, seed, ledger=None):
    if ledger is not None:
        ledger.add(description)
    material = f"{PROTOCOL}|{description}|{seed}".encode()
    return np.random.default_rng(int.from_bytes(hashlib.sha256(material).digest()[:16], "little"))


def largest_remainder(counts, size):
    """Exact integer remainders; caller supplies keys in the prescribed tie order."""
    counts = np.asarray(counts, dtype=np.int64)
    if counts.ndim != 1 or (counts < 0).any() or size < 0 or size > int(counts.sum()):
        raise ValueError("Impossible largest-remainder sample")
    if size == 0:
        return np.zeros_like(counts)
    total = int(counts.sum())
    quota = (counts * size) // total
    remainder = (counts * size) % total
    remaining = size - int(quota.sum())
    order = np.lexsort((np.arange(len(counts)), -remainder))
    quota[order[:remaining]] += 1
    if int(quota.sum()) != size or np.any(quota > counts):
        raise AssertionError("Largest-remainder capacity failure")
    return quota


def proportional_sample(cells, strata, size, description, seed, ledger=None):
    cells, strata = np.asarray(cells, dtype=np.int64), np.asarray(strata)
    if len(cells) != len(strata) or len(np.unique(cells)) != len(cells):
        raise ValueError("Misaligned or duplicated sampling population")
    keys = sorted(set(strata.tolist()))
    groups = [np.sort(cells[strata == key]) for key in keys]
    quota = largest_remainder([len(x) for x in groups], size)
    parts = []
    for key, group, n in zip(keys, groups, quota):
        if n:
            stream = f"{description}|stratum={key}"
            parts.append(rng_for(stream, seed, ledger).permutation(group)[:int(n)])
    result = np.sort(np.concatenate(parts)) if parts else np.array([], dtype=np.int64)
    if len(result) != size:
        raise AssertionError("Incorrect selected-cell count")
    return result


def matched_controls(pert_cells, control_pool, gems, description, seed, ledger=None):
    parts = []
    for gem in sorted(set(map(int, gems[pert_cells]))):
        n = int(np.sum(gems[pert_cells] == gem))
        choices = np.sort(control_pool[gems[control_pool] == gem])
        if len(choices) < n:
            raise ValueError(f"Insufficient controls before expression access: {description}, GEM {gem}, need {n}, have {len(choices)}")
        parts.append(rng_for(f"{description}|gem={gem}", seed, ledger).permutation(choices)[:n])
    result = np.sort(np.concatenate(parts))
    if len(np.unique(result)) != len(pert_cells):
        raise AssertionError("Duplicated matched controls")
    return result


def target_folds(targets):
    ordered = sorted(targets, key=lambda t: (hashlib.sha256(f"{PROTOCOL}|fold|{t}".encode()).hexdigest(), t))
    return {target: i % 3 for i, target in enumerate(ordered)}


def normalized_counts(counts):
    counts = np.asarray(counts)
    if not np.isfinite(counts).all() or np.any(counts < 0) or np.any(counts != np.floor(counts)):
        raise ValueError("Invalid raw counts in selected cells")
    totals = counts.sum(axis=1, dtype=np.float64)
    if np.any(totals <= 0):
        raise ValueError("Frozen selected cell has zero observed-gene library; no replacement permitted")
    return np.log1p(counts.astype(np.float64) * (10000.0 / totals[:, None])), totals


def read_rows(dataset, cells, checked, total_by_cell):
    """Read only requested sorted rows, at most 256 per HDF5 request."""
    cells = np.asarray(cells, dtype=np.int64)
    if len(np.unique(cells)) != len(cells):
        raise ValueError("Repeated cells in an expression arm")
    order = np.argsort(cells)
    sorted_cells = cells[order]
    blocks = []
    for start in range(0, len(cells), 256):
        selected = sorted_cells[start:start + 256]
        x = dataset[selected, :]
        if np.any(x > np.iinfo(np.uint32).max):
            raise ValueError("Raw count exceeds exact uint32 cache capacity")
        _, totals = normalized_counts(x)
        checked.update(map(int, selected))
        total_by_cell.update({int(i): float(n) for i, n in zip(selected, totals)})
        blocks.append(x.astype(np.uint32))
    return np.vstack(blocks)[np.argsort(order)]


def cell_rows(cells, metadata, role, target=None, fold=None, partition=None):
    frame = metadata.loc[cells, ["cell_index", "cell_barcode", "gem_group", "gene_id", "gene_transcript", "sgID_AB"]].copy()
    frame["role"] = role
    if target is not None: frame["target_id"] = target
    if fold is not None: frame["fold"] = fold
    if partition is not None: frame["null_partition"] = partition
    return frame


def metadata_stage(f, cfg, source):
    print("Stage 1: freeze cells and folds using metadata only", flush=True)
    seed = cfg["split_seed"]
    ledger = set()
    genes = read_column(f["var"], "gene_id").astype(str)
    names = read_column(f["var"], "gene_name").astype(str)
    ids = read_column(f["obs"], "gene_id").astype(str)
    labels = read_column(f["obs"], "gene_transcript").astype(str)
    guides = read_column(f["obs"], "sgID_AB").astype(str)
    gems = np.asarray(read_column(f["obs"], "gem_group"), dtype=np.int64)
    bars = read_column(f["obs"], "cell_barcode").astype(str)
    if len(genes) != 8749 or len(set(genes)) != len(genes) or len(set(bars)) != len(bars):
        raise ValueError("Unexpected source gene or cell identity")
    metadata = pd.DataFrame({"cell_index": np.arange(len(ids)), "cell_barcode": bars, "gem_group": gems,
                             "gene_id": ids, "gene_transcript": labels, "sgID_AB": guides})
    control = np.flatnonzero(np.char.find(labels, "non-targeting") >= 0)
    stable = np.array([re.fullmatch(r"ENSG[0-9]{11}", g) is not None for g in ids])
    if any(label.rsplit("_", 1)[-1] != gene for label, gene in zip(labels[stable], ids[stable])):
        raise ValueError("Stable target IDs and perturbation labels disagree")
    all_targets = sorted(set(ids[stable]))
    by_target = {t: np.flatnonzero(ids == t) for t in all_targets}
    arm_ids = (gems - 1) % 3
    eligibility = []
    measured = set(genes)
    for target, cells in by_target.items():
        n = [int(np.sum(arm_ids[cells] == i)) for i in range(3)]
        reasons = (["target_not_in_measured_genes"] if target not in measured else [])
        reasons += [f"{arm}_fewer_than_{need}" for arm, need, count in zip("ABC", [40, 30, 30], n) if count < need]
        eligibility.append({"target_id": target, "total_cells": len(cells), "cells_A": n[0], "cells_B": n[1],
                            "cells_C": n[2], "target_in_measured_genes": target in measured,
                            "eligible": not reasons, "exclusion_reasons": ";".join(reasons)})
    eligibility = pd.DataFrame(eligibility)
    targets = sorted(eligibility.loc[eligibility.eligible, "target_id"])
    if len(targets) != 505 or len(control) != 11485 or set(gems) != set(range(1, 57)):
        raise ValueError("Metadata preflight does not reproduce frozen 505-target/control/GEM scope")
    folds = target_folds(targets)
    dm = proportional_sample(control, gems[control], 400, "control_pool|Dm", seed, ledger)
    remaining = np.setdiff1d(control, dm)
    dc = proportional_sample(remaining, gems[remaining], 400, "control_pool|Dc", seed, ledger)
    remaining = np.setdiff1d(remaining, dc)
    pools = {arm: remaining[arm_ids[remaining] == i] for i, arm in enumerate("ABC")}
    dc64_indices = np.sort(rng_for("covariate_subset|Dc64", seed, ledger).choice(400, 64, replace=False))
    dc64 = dc[dc64_indices]
    if len(set(dm) | set(dc) | set(remaining)) != len(control):
        raise AssertionError("Control partition overlap")
    sampled = {}
    split_frames = []
    for target in targets:
        cells = by_target[target]
        sampled[target] = {}
        for i, (arm, n) in enumerate(zip("ABC", [40, 30, 30])):
            available = cells[arm_ids[cells] == i]
            strata = np.asarray([json.dumps([guides[c], int(gems[c])], separators=(",", ":")) for c in available])
            selected = proportional_sample(available, strata, n, f"target|{target}|arm={arm}", seed, ledger)
            controls = matched_controls(selected, pools[arm], gems, f"matched_control|{target}|arm={arm}", seed, ledger)
            sampled[target][arm] = (selected, controls)
            split_frames.extend([cell_rows(selected, metadata, f"{arm}_pert", target, folds[target]),
                                 cell_rows(controls, metadata, f"{arm}_ctrl", target, folds[target])])
        all_six = np.concatenate([x for arm in "ABC" for x in sampled[target][arm]])
        if len(np.unique(all_six)) != 200 or np.intersect1d(all_six, np.r_[dm, dc]).size:
            raise AssertionError("Target arms overlap each other or D controls")
    nulls, null_frames = [], []
    for partition in range(cfg["null_partitions"]):
        arms = {}
        for arm in "BC":
            pseudo = np.sort(rng_for(f"null|{partition}|arm={arm}|pseudo", seed, ledger).choice(pools[arm], 30, replace=False))
            pool = np.setdiff1d(pools[arm], pseudo)
            matched = matched_controls(pseudo, pool, gems, f"null|{partition}|arm={arm}|matched", seed, ledger)
            arms[arm] = (pseudo, matched)
            null_frames.extend([cell_rows(pseudo, metadata, f"{arm}_pseudo", partition=partition),
                                cell_rows(matched, metadata, f"{arm}_ctrl", partition=partition)])
        if len(np.unique(np.concatenate([v for arm in "BC" for v in arms[arm]]))) != 120:
            raise AssertionError("Null arms overlap")
        nulls.append(arms)
    genes_frame = pd.DataFrame({"gene_id": genes, "gene_name": names, "source_gene_index": np.arange(len(genes))})
    genes_frame.to_csv(OUT / "genes_all.csv", index=False)
    eligibility.to_csv(OUT / "target_eligibility.csv", index=False)
    ambiguous = metadata.loc[(~stable) & (~metadata.cell_index.isin(control)), "gene_transcript"].value_counts()
    ambiguous.rename_axis("gene_transcript").reset_index(name="cell_count").to_csv(OUT / "excluded_noncontrol_labels.csv", index=False)
    targets_frame = eligibility.set_index("target_id").loc[targets].reset_index()
    targets_frame.insert(1, "fold", [folds[t] for t in targets])
    targets_frame.insert(2, "all_target_index", np.arange(len(targets)))
    targets_frame.to_csv(OUT / "targets_all.csv", index=False)
    pd.concat(split_frames, ignore_index=True).to_csv(OUT / "target_cell_splits.csv.gz", index=False)
    pd.concat(null_frames, ignore_index=True).to_csv(OUT / "null_cell_splits.csv.gz", index=False)
    cell_rows(dm, metadata, "Dm").to_csv(OUT / "Dm_cells.csv", index=False)
    cell_rows(dc, metadata, "Dc").assign(in_Dc64=np.isin(dc, dc64)).to_csv(OUT / "Dc_cells.csv", index=False)
    np.save(OUT / "Dc64_indices.npy", dc64_indices, allow_pickle=False)
    pd.concat([cell_rows(dm, metadata, "Dm"), cell_rows(dc, metadata, "Dc"),
               *[cell_rows(pools[a], metadata, a) for a in "ABC"]], ignore_index=True).to_csv(OUT / "control_pools.csv.gz", index=False)
    with gzip.open(OUT / "rng_stream_descriptions.txt.gz", "wt") as handle:
        handle.write("\n".join(sorted(ledger)) + "\n")
    write_json(OUT / "metadata_manifest.json", {"stage": "metadata_only_before_any_new_X_read", "protocol_sha256": sha(CONFIG),
        "source_sha256": source["sha256"], "targets": len(targets), "fold_counts": {str(k): sum(f == k for f in folds.values()) for k in range(3)},
        "Dm_cells": len(dm), "Dc_cells": len(dc), "Dc64_cells": len(dc64),
        "remaining_control_pools": {a: len(pools[a]) for a in "ABC"}, "null_partitions": len(nulls),
        "guide_GEM_stratum_encoding": "compact JSON array [sgID_AB, integer gem_group]; lexical string tie order",
        "GEM_quota_tie_order": "ascending numeric GEM ID", "all_matched_control_capacities_checked": True,
        "independent_rng_streams": len(ledger), "output_sha256": {p.name: sha(p) for p in OUT.iterdir() if p.is_file()}})
    return dict(genes=genes, genes_frame=genes_frame, targets=targets, folds=folds, all_source_targets=all_targets,
                metadata=metadata, dm=dm, dc=dc, dc64_indices=dc64_indices, sampled=sampled, nulls=nulls)


def dm_freeze_stage(f, data, checked, totals):
    print("Stage 2: Dm-only feature and evaluation-target selection", flush=True)
    raw = read_rows(f["X"], data["dm"], checked, totals)
    normalized, _ = normalized_counts(raw)
    detection, variance = (raw > 0).mean(axis=0), normalized.var(axis=0)
    selectable = (detection >= 0.05) & (variance > 0)
    gene_to_col = {g: i for i, g in enumerate(data["genes"])}
    candidates = np.flatnonzero(selectable)
    selected = candidates[np.lexsort((data["genes"][candidates], -variance[candidates]))][:1000]
    eligible_targets = [t for t in data["targets"] if selectable[gene_to_col[t]]]
    eval_targets = sorted(eligible_targets, key=lambda t: (hashlib.sha256(f"{PROTOCOL}|evaluation_target|{t}".encode()).hexdigest(), t))[:128]
    fold_counts = {k: sum(data["folds"][t] == k for t in eval_targets) for k in range(3)}
    if len(eval_targets) < 64 or min(fold_counts.values()) < 10:
        raise ValueError("Dm-only evaluation target selection fails minimum number/fold support")
    model_cols = np.sort(np.unique(np.r_[selected, [gene_to_col[t] for t in eval_targets]]))
    evaluation_cols = np.sort([i for i in selected if data["genes"][i] not in set(data["all_source_targets"])])
    if len(evaluation_cols) <= 50:
        raise ValueError("Fixed evaluation gene universe is not larger than Top50")
    data.update(dm_raw=raw, dm_normalized_mean=normalized.mean(axis=0), evaluation_targets=eval_targets,
                evaluation_cols=evaluation_cols, model_cols=model_cols)
    gene_table = data["genes_frame"].copy().assign(Dm_detection=detection, Dm_variance=variance, selectable=selectable,
               selected_top_variance=np.isin(np.arange(len(variance)), selected), in_model_scope=np.isin(np.arange(len(variance)), model_cols),
               in_evaluation_scope=np.isin(np.arange(len(variance)), evaluation_cols))
    gene_table.to_csv(OUT / "Dm_gene_selection.csv", index=False)
    data["genes_frame"].iloc[model_cols].to_csv(OUT / "model_genes.csv", index=False)
    data["genes_frame"].iloc[evaluation_cols].to_csv(OUT / "evaluation_genes.csv", index=False)
    pd.DataFrame({"target_id": eval_targets, "fold": [data["folds"][t] for t in eval_targets],
                  "all_target_index": [data["targets"].index(t) for t in eval_targets],
                  "source_gene_index": [gene_to_col[t] for t in eval_targets]}).to_csv(OUT / "evaluation_targets.csv", index=False)
    pd.DataFrame({"target_id": data["targets"], "Dm_selectable": [selectable[gene_to_col[t]] for t in data["targets"]],
                  "evaluation_selected": [t in set(eval_targets) for t in data["targets"]],
                  "selection_reason": ["selected_hash_order" if t in set(eval_targets) else "hash_rank_beyond_128" if t in eligible_targets else "Dm_detection_or_variance" for t in data["targets"]]}).to_csv(OUT / "evaluation_target_selection.csv", index=False)
    save_npz(OUT / "Dm_counts.npz", csr_matrix(raw), compressed=True)
    np.save(OUT / "Dm_q1_normalized_mean.npy", data["dm_normalized_mean"], allow_pickle=False)
    write_json(OUT / "gene_freeze_manifest.json", {"stage": "Dm_only_before_target_or_Dc_expression_read", "protocol_sha256": sha(CONFIG),
        "metadata_manifest_sha256": sha(OUT / "metadata_manifest.json"), "Dm_cells": len(raw),
        "selectable_genes": int(selectable.sum()), "top_variance_genes": len(selected), "model_genes": len(model_cols),
        "evaluation_genes": len(evaluation_cols), "Dm_selectable_targets": len(eligible_targets), "evaluation_targets": len(eval_targets),
        "evaluation_targets_per_fold": {str(k): v for k, v in fold_counts.items()},
        "output_sha256": {name: sha(OUT / name) for name in ["Dm_gene_selection.csv", "model_genes.csv", "evaluation_genes.csv",
                            "evaluation_targets.csv", "evaluation_target_selection.csv", "Dm_counts.npz", "Dm_q1_normalized_mean.npy"]}})


def response_stage(f, data, checked, totals):
    print("Stage 3: frozen response construction; no model fitting or thinning", flush=True)
    raw_dc = read_rows(f["X"], data["dc"], checked, totals)
    norm_dc, _ = normalized_counts(raw_dc)
    save_npz(OUT / "Dc_counts.npz", csr_matrix(raw_dc), compressed=True)
    pd.DataFrame({"gene_id": data["genes"], "mean_Dc400": norm_dc.mean(axis=0),
                  "mean_Dc64": norm_dc[data["dc64_indices"]].mean(axis=0)}).to_csv(OUT / "q1_covariates.csv", index=False)
    n_targets, n_genes = len(data["targets"]), len(data["genes"])
    a_post, a_control = np.zeros((n_targets, n_genes)), np.zeros((n_targets, n_genes))
    a_indices = np.stack([data["sampled"][t]["A"][0] for t in data["targets"]])
    with h5py.File(OUT / "A_raw_counts.h5", "w") as cache:
        cache.attrs["source_sha256"] = "9b05ef1f81526216fa008d677e9e0d03dce9a2f7a95499a4fb81e505e9d88ef1"
        cache.attrs["normalization"] = "none; exact observed raw counts"
        cache.attrs["source_cell_order"] = "ascending source row within each frozen target A arm"
        ds = cache.create_dataset("counts", shape=(n_targets, 40, n_genes), dtype="uint32", chunks=(1, 40, n_genes), compression="gzip", compression_opts=1, shuffle=True)
        cache.create_dataset("targets", data=np.asarray(data["targets"], dtype=object), dtype=h5py.string_dtype())
        cache.create_dataset("genes", data=np.asarray(data["genes"], dtype=object), dtype=h5py.string_dtype())
        cache.create_dataset("cell_indices", data=a_indices)
        for i, target in enumerate(data["targets"]):
            pert, ctrl = data["sampled"][target]["A"]
            raw = read_rows(f["X"], pert, checked, totals)
            ds[i] = raw
            a_post[i] = normalized_counts(raw)[0].mean(axis=0)
            a_control[i] = normalized_counts(read_rows(f["X"], ctrl, checked, totals))[0].mean(axis=0)
            if (i + 1) % 100 == 0 or i + 1 == n_targets:
                print(f"Frozen A means and raw cache: {i+1}/{n_targets} targets", flush=True)
    np.save(OUT / "A_post_mean.npy", a_post, allow_pickle=False)
    np.save(OUT / "A_control_mean.npy", a_control, allow_pickle=False)
    m = len(data["evaluation_targets"])
    deltas = {arm: np.zeros((m, n_genes)) for arm in "BC"}
    for i, target in enumerate(data["evaluation_targets"]):
        for arm in "BC":
            pert, ctrl = data["sampled"][target][arm]
            deltas[arm][i] = normalized_counts(read_rows(f["X"], pert, checked, totals))[0].mean(axis=0) - normalized_counts(read_rows(f["X"], ctrl, checked, totals))[0].mean(axis=0)
    null_delta = {arm: np.zeros((len(data["nulls"]), n_genes)) for arm in "BC"}
    for i, task in enumerate(data["nulls"]):
        for arm in "BC":
            pseudo, ctrl = task[arm]
            null_delta[arm][i] = normalized_counts(read_rows(f["X"], pseudo, checked, totals))[0].mean(axis=0) - normalized_counts(read_rows(f["X"], ctrl, checked, totals))[0].mean(axis=0)
    # The zero-library rule covers every frozen cell, including B/C cells of
    # targets not selected for evaluation. Only library QC is retained for those.
    all_frozen = set(map(int, np.r_[data["dm"], data["dc"]]))
    for task in list(data["sampled"].values()) + data["nulls"]:
        for arms in task.values():
            for cells in arms: all_frozen.update(map(int, cells))
    unchecked = np.array(sorted(all_frozen - checked), dtype=np.int64)
    for start in range(0, len(unchecked), 256):
        read_rows(f["X"], unchecked[start:start+256], checked, totals)
    if checked != all_frozen:
        raise AssertionError("Selected-cell library QC coverage mismatch")
    pd.DataFrame({"cell_index": sorted(checked), "observed_gene_UMI_sum": [totals[i] for i in sorted(checked)]}).to_csv(OUT / "selected_cell_library_checks.csv.gz", index=False)
    eval_indices = np.array([data["targets"].index(t) for t in data["evaluation_targets"]])
    cols = data["evaluation_cols"]
    np.savez_compressed(OUT / "evaluation_responses.npz", targets=np.asarray(data["evaluation_targets"]),
        folds=np.array([data["folds"][t] for t in data["evaluation_targets"]]), genes=data["genes"][cols],
        source_gene_indices=cols, all_target_indices=eval_indices, A=(a_post[eval_indices]-a_control[eval_indices])[:, cols],
        B=deltas["B"][:, cols], C=deltas["C"][:, cols], null_B=null_delta["B"][:, cols], null_C=null_delta["C"][:, cols])
    return {"all_frozen_cells_checked": len(checked), "zero_library_cells": 0, "A_targets": n_targets,
            "evaluation_targets": m, "null_partitions": len(data["nulls"]), "observed_genes": n_genes}


def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise SystemExit("Preparation directory is nonempty; refusing to overwrite frozen stages")
    if sha(CONFIG) != EXPECTED_PROTOCOL_SHA256:
        raise SystemExit("Frozen protocol SHA256 changed; explicit review required")
    cfg = json.loads(CONFIG.read_text())
    start = datetime.now(timezone.utc).isoformat()
    audit_dir = HERE / "provenance/rpe1_contents_v1"
    audit = json.loads((audit_dir / "summary.json").read_text())
    audit_manifest = json.loads((audit_dir / "manifest.json").read_text())
    if not audit["passed"] or sha(audit_dir / "summary.json") != audit_manifest["output_sha256"]["summary.json"]:
        raise SystemExit("Verified contents audit is absent or has changed")
    if audit_manifest["sources"][0]["sha256"] != cfg["source"]["sha256"]:
        raise SystemExit("Contents audit identifies a different singlecell source")
    source_path = HERE / cfg["source"]["file"]
    source = source_gate(source_path, HERE / "provenance/download_receipts/Replogle_RPE1_raw_singlecell.json",
                         cfg["source"]["bytes"], cfg["source"]["md5"], 35775606)
    if source["sha256"] != cfg["source"]["sha256"]:
        raise SystemExit("Source SHA256 does not match frozen protocol")
    code_paths = [Path(__file__), HERE / "verify_rpe1_contents.py"]
    code_hashes = {p.name: sha(p) for p in code_paths}
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        checked, totals = set(), {}
        with h5py.File(source_path, "r") as f:
            data = metadata_stage(f, cfg, source)
            dm_freeze_stage(f, data, checked, totals)
            summary = response_stage(f, data, checked, totals)
        if sha(CONFIG) != EXPECTED_PROTOCOL_SHA256 or any(sha(p) != code_hashes[p.name] for p in code_paths):
            raise ValueError("Protocol or preparation code changed during execution")
        stat = source_path.stat()
        if (stat.st_size, stat.st_mtime_ns) != (source["bytes"], source["mtime_ns"]):
            raise ValueError("Source file changed during preparation")
        write_json(OUT / "manifest.json", {"status": "completed_frozen_preparation_no_model_or_thinning",
            "started_utc": start, "completed_utc": datetime.now(timezone.utc).isoformat(),
            "protocol_sha256": EXPECTED_PROTOCOL_SHA256, "protocol": cfg, "source_sha256": source["sha256"],
            "source_verification": source, "contents_audit_manifest_sha256": sha(audit_dir / "manifest.json"),
            "code_sha256": code_hashes, "runtime": {"python": sys.version, "numpy": np.__version__, "h5py": h5py.__version__},
            "summary": summary, "output_sha256": {str(p.relative_to(OUT)): sha(p) for p in sorted(OUT.iterdir()) if p.is_file()},
            "stage_order": ["metadata IDs and folds saved", "Dm-only genes/evaluation targets saved", "Dc and frozen A/B/C/null responses read"],
            "model_access_policy": "Model fitting may use Dm raw counts, A raw counts, A post means, frozen IDs/folds and gene lists. No B/C/null outcome file may be used for model fitting or selection.",
            "limits": ["GEM-disjoint technical evaluation, not biological replicates", "Guide-pair mixtures remain", "Observed gene universe was author-filtered using the original study"]})
        print(json.dumps(summary, indent=2), flush=True)
    except Exception as exc:
        write_json(OUT / "failure.json", {"status": "preparation_failed_without_replacement", "error": repr(exc),
            "protocol_sha256": EXPECTED_PROTOCOL_SHA256, "code_sha256": code_hashes})
        raise


if __name__ == "__main__":
    main()
