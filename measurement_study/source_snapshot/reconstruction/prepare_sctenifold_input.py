"""Prepare control-only raw-count matrices for a scTenifoldKnk pilot.

Feature filtering and target eligibility use only the independent D control
pool. Held-out perturbation and evaluation cells are not used to select genes.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

from io_utils import HERE, column, gene_symbols, input_path, sha256, write_json


DATASETS = ["Adamson", "Norman"]
SEED = 20261024
N_HVG = 1000
MIN_TARGET_DETECTION = 0.05
N_TRAIN_CELLS = 400


def decode(values):
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values])


def load_counts_subset(h5path, row_indices):
    with h5py.File(h5path, "r") as f:
        group = f["layers/counts"]
        shape = tuple(int(x) for x in group.attrs["shape"])
        full = csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)
        return full[np.asarray(row_indices, dtype=int), :].tocsr()


def prepare_dataset(dataset):
    run_dir = HERE / "results/library_calibration_v2" / f"{dataset}_disjoint_library_pools_seed{SEED}"
    response_path = run_dir / "responses.npz"
    response_manifest = run_dir / "manifest.json"
    split_path = run_dir / "cell_splits.csv.gz"
    source = input_path(dataset)
    provenance_path = HERE / f"provenance/{dataset}_file_verification.json"
    receipt = json.loads(provenance_path.read_text())
    run_receipt = json.loads(response_manifest.read_text())
    if not receipt.get("public_bytes_match") or sha256(source) != receipt.get("local_sha256"):
        raise ValueError(f"Source file no longer matches the verified public input: {dataset}")
    if sha256(response_path) != run_receipt.get("outputs", {}).get("responses.npz"):
        raise ValueError(f"Frozen response input hash mismatch: {dataset}")
    out = HERE / "results/model_rankings_sctenifoldknk_pilot" / dataset
    if (out / "preparation_manifest.json").exists():
        raise SystemExit(f"Refusing to overwrite prepared inputs: {out}")
    out.mkdir(parents=True, exist_ok=True)
    split = pd.read_csv(split_path)
    d = split.loc[split.arm == "covariates_D"].copy()
    if d.cell_index.duplicated().any():
        raise ValueError("D control cells are not unique")
    z = np.load(response_path, allow_pickle=False)
    genes = z["genes"].astype(str)
    tasks = z["tasks"].astype(str)
    targets_all = tasks[~np.char.startswith(tasks, "control_vs_control_")]
    if len(set(genes)) != len(genes):
        raise ValueError("Duplicate gene symbols in the frozen source gene universe")
    if not set(targets_all).issubset(set(genes)):
        raise ValueError("At least one perturbation target is absent from the gene universe")

    counts_d = load_counts_subset(source, d.cell_index.to_numpy(int))
    with h5py.File(source, "r") as f:
        if not np.array_equal(gene_symbols(f).astype(str), genes):
            raise ValueError("Frozen response gene order no longer matches source H5AD")
        if len(d.cell_index) and int(d.cell_index.max()) >= int(f["X"].shape[0]):
            raise ValueError("D control cell index exceeds source H5AD")
    detected = np.asarray((counts_d > 0).mean(axis=0)).ravel()
    eligible_target = np.isin(genes, targets_all) & (detected >= MIN_TARGET_DETECTION)
    eligible_targets = genes[eligible_target]
    if len(eligible_targets) < 20:
        raise ValueError("Too few targets pass the control-only detection threshold")

    # Rank variability on log1p CPM computed from D controls only.
    library_size = np.asarray(counts_d.sum(axis=1)).ravel()
    normalized = counts_d.multiply(1e4 / np.maximum(library_size, 1)[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    variance = np.asarray(normalized.power(2).mean(axis=0)).ravel() - np.square(
        np.asarray(normalized.mean(axis=0)).ravel()
    )
    eligible_gene = detected >= MIN_TARGET_DETECTION
    candidates = np.flatnonzero(eligible_gene)
    # Stable tie breaking by symbol is pre-outcome and reproducible.
    order = candidates[np.lexsort((genes[candidates], -variance[candidates]))]
    top_hvg = order[:min(N_HVG, len(order))]
    selected = np.unique(np.concatenate([top_hvg, np.flatnonzero(eligible_target)]))
    selected_genes = genes[selected]

    # Keep all eligible D controls when possible; otherwise sample 400 per GEM
    # group with fixed seed. This step uses no perturbation/effect outcomes.
    rng = np.random.default_rng(SEED + (1 if dataset == "Adamson" else 2))
    if len(d) <= N_TRAIN_CELLS:
        chosen = np.arange(len(d))
    else:
        groups = []
        levels = sorted(d.gemgroup.unique())
        allocation = {g: N_TRAIN_CELLS // len(levels) for g in levels}
        for g in levels[:N_TRAIN_CELLS % len(levels)]:
            allocation[g] += 1
        for g in levels:
            local = np.flatnonzero(d.gemgroup.to_numpy() == g)
            n = min(allocation[g], len(local))
            groups.extend(rng.choice(local, n, replace=False).tolist())
        chosen = np.asarray(sorted(groups), dtype=int)
    if len(chosen) < 100:
        raise ValueError("Too few independent D controls for model fitting")

    train_cell_indices = d.cell_index.to_numpy(int)[chosen]
    X = counts_d[chosen, :][:, selected].T.tocsc()
    if X.shape != (len(selected_genes), len(train_cell_indices)):
        raise AssertionError("Model input dimensions are misaligned")
    mmwrite(out / "control_counts.mtx", X)
    (out / "genes.txt").write_text("\n".join(selected_genes) + "\n")
    (out / "targets.txt").write_text("\n".join(eligible_targets) + "\n")
    (out / "all_real_targets.txt").write_text("\n".join(targets_all) + "\n")
    pd.DataFrame({"cell_index": train_cell_indices,
                  "gemgroup": d.gemgroup.to_numpy()[chosen]}).to_csv(out / "training_cells.csv", index=False)
    pd.DataFrame({"gene_symbol": genes, "detection_fraction_D": detected,
                  "variance_log1p_cpm_D": variance,
                  "selected_for_model": np.isin(np.arange(len(genes)), selected),
                  "target_eligible": eligible_target}).to_csv(out / "gene_selection.csv", index=False)
    config = {
        "dataset": dataset, "analysis_type": "model_output_pilot_development_only",
        "source_file": str(source), "source_sha256": sha256(source),
        "source_receipt_sha256": sha256(provenance_path),
        "response_manifest_sha256": sha256(response_manifest),
        "split_sha256": sha256(split_path),
        "raw_counts_layer": "layers/counts (original sparse integer UMI counts)",
        "model_cells": int(len(train_cell_indices)), "available_D_controls": int(len(d)),
        "training_gemgroups": sorted(set(map(int, d.gemgroup.to_numpy()[chosen]))),
        "feature_rule": f"Top {N_HVG} log1p(CPM) variance genes in D controls, plus eligible perturbation targets; all genes require >=5% detection in D",
        "minimum_target_detection_fraction": MIN_TARGET_DETECTION,
        "eligible_targets": int(len(eligible_targets)), "all_real_targets": int(len(targets_all)),
        "selected_gene_count": int(len(selected_genes)), "primary_k": 50,
        "seed": SEED,
        "preparation_script_sha256": sha256(Path(__file__)),
        "limitations": [
            "Only control D cells fit the model; B and C perturbation outcomes are held out for evaluation.",
            "This is a technical single-cell split from one study, not an independent biological replicate.",
            "Feature scope is a D-only variance subset plus targets and does not represent whole-transcriptome coverage.",
            "Transcriptome-wide heat distance is a native scTenifoldKnk output but differs from its default single-target manifold mode.",
        ],
    }
    config["input_sha256"] = {name: sha256(out / name) for name in [
        "control_counts.mtx", "genes.txt", "targets.txt", "all_real_targets.txt",
        "training_cells.csv", "gene_selection.csv"]}
    write_json(out / "preparation_manifest.json", config)
    print(dataset, "selected genes", len(selected_genes), "eligible targets", len(eligible_targets),
          "/", len(targets_all), "training cells", len(train_cell_indices), flush=True)


def main():
    for dataset in DATASETS:
        prepare_dataset(dataset)


if __name__ == "__main__":
    main()
