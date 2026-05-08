#!/usr/bin/env python3
"""Download and standardize public Perturb-seq datasets for CARF-Benchmark.

The script has two layers:

* source activation: create source manifests from stable public h5ad records;
* data preparation: convert a downloaded h5ad into CARF inputs.

Large h5ad files are not downloaded in CI by default. A dataset becomes a full
leaderboard run when its h5ad has been downloaded and this script has written
`inputs/{ground_truth,covariates,method_scores}.csv`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


CARF_ROOT = Path(__file__).resolve().parents[1]
SOURCE_MANIFEST = CARF_ROOT / "registry" / "public_dataset_sources.csv"
DEFAULT_RAW_DIR = CARF_ROOT / "raw"
DEFAULT_RUNS_DIR = CARF_ROOT / "runs"
PERTURBATION_KEY_CANDIDATES = [
    "condition",
    "perturbation",
    "perturbation_name",
    "perturbation_id",
    "gene",
    "target",
    "target_gene",
    "guide_target",
    "sgRNA",
    "sgRNA_gene",
    "cov_drug_dose_name",
    "treatment",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare public Perturb-seq CARF inputs.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("list-sources", help="Print configured public data sources.")

    materialize = subparsers.add_parser(
        "materialize-sources",
        help="Write source_manifest.json files under carf_benchmark/runs/*.",
    )
    materialize.add_argument("--dataset-id", default="all")
    materialize.add_argument("--runs-dir", type=Path, default=DEFAULT_RUNS_DIR)

    download = subparsers.add_parser("download", help="Download one or all source h5ad files.")
    download.add_argument("--dataset-id", required=True)
    download.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    download.add_argument("--decompress", action="store_true", default=True)
    download.add_argument("--no-decompress", action="store_false", dest="decompress")

    inspect = subparsers.add_parser("inspect", help="Inspect h5ad obs/var fields.")
    inspect.add_argument("--h5ad", required=True, type=Path)

    prepare = subparsers.add_parser("prepare", help="Convert a downloaded h5ad into CARF inputs.")
    prepare.add_argument("--dataset-id", required=True)
    prepare.add_argument("--h5ad", required=True, type=Path)
    prepare.add_argument("--output-run-dir", type=Path)
    prepare.add_argument("--perturbation-key")
    prepare.add_argument("--gene-symbol-key")
    prepare.add_argument("--control-regex")
    prepare.add_argument("--perturbations")
    prepare.add_argument("--max-perturbations", type=int, default=3)
    prepare.add_argument("--min-cells", type=int, default=25)
    prepare.add_argument("--max-cells-per-group", type=int, default=1000)
    prepare.add_argument("--max-genes", type=int)
    prepare.add_argument("--fdr-alpha", type=float, default=0.05)
    prepare.add_argument("--min-abs-delta", type=float, default=0.0)
    prepare.add_argument("--seed", type=int, default=1)

    return parser.parse_args()


def load_sources() -> pd.DataFrame:
    return pd.read_csv(SOURCE_MANIFEST)


def optional_string(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value)
    if text.lower() in {"nan", "none", "null"}:
        return ""
    return text


def select_sources(dataset_id: str) -> pd.DataFrame:
    sources = load_sources()
    if dataset_id == "all":
        return sources
    selected = sources[sources["dataset_id"] == dataset_id]
    if selected.empty:
        raise SystemExit(f"Unknown dataset_id: {dataset_id}")
    return selected


def list_sources() -> None:
    cols = ["dataset_id", "file_name", "compressed_bytes", "source_doi", "primary_accession"]
    print(load_sources()[cols].to_string(index=False))


def materialize_sources(dataset_id: str, runs_dir: Path) -> None:
    for _, row in select_sources(dataset_id).iterrows():
        run_dir = runs_dir / row["dataset_id"]
        run_dir.mkdir(parents=True, exist_ok=True)
        manifest = {
            "dataset_id": row["dataset_id"],
            "source_dataset_name": row["source_dataset_name"],
            "source_record": row["source_record"],
            "source_doi": row["source_doi"],
            "primary_accession": row["primary_accession"],
            "file_name": row["file_name"],
            "url": row["url"],
            "compressed_md5": row["compressed_md5"],
            "compressed_bytes": int(row["compressed_bytes"]),
            "license": row["license"],
            "activation_state": "source_manifest_active_inputs_pending",
            "notes": row["notes"],
        }
        path = run_dir / "source_manifest.json"
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote {path}")


def download_sources(dataset_id: str, raw_dir: Path, decompress: bool) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    for _, row in select_sources(dataset_id).iterrows():
        dataset_dir = raw_dir / row["dataset_id"]
        dataset_dir.mkdir(parents=True, exist_ok=True)
        compressed_path = dataset_dir / row["file_name"]
        if not compressed_path.exists():
            print(f"Downloading {row['dataset_id']} from {row['source_doi']}...")
            subprocess.run(
                ["curl", "-L", "--fail", "--retry", "3", "-o", str(compressed_path), row["url"]],
                check=True,
            )
        verify_md5(compressed_path, row["compressed_md5"])
        if decompress and compressed_path.name.endswith(".gz"):
            decompressed_path = compressed_path.with_suffix("")
            if not decompressed_path.exists():
                print(f"Decompressing {compressed_path.name}...")
                with gzip.open(compressed_path, "rb") as source, decompressed_path.open("wb") as target:
                    shutil.copyfileobj(source, target)
            print(f"Ready: {decompressed_path}")
        else:
            print(f"Ready: {compressed_path}")


def verify_md5(path: Path, expected: str) -> None:
    expected = expected.replace("md5:", "")
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    actual = h.hexdigest()
    if actual != expected:
        raise SystemExit(f"MD5 mismatch for {path}: expected {expected}, got {actual}")


def read_h5ad(path: Path):
    try:
        import anndata as ad
    except ImportError as exc:
        raise SystemExit("Preparing public h5ad datasets requires anndata.") from exc
    return ad.read_h5ad(path)


def inspect_h5ad(path: Path) -> None:
    adata = read_h5ad(path)
    print(f"shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print("obs columns:")
    for col in adata.obs.columns:
        nunique = adata.obs[col].nunique(dropna=True)
        print(f"  {col}: {nunique} unique")
    print("var columns:")
    for col in adata.var.columns:
        print(f"  {col}")


def infer_perturbation_key(obs: pd.DataFrame) -> str:
    lower_to_actual = {col.lower(): col for col in obs.columns}
    for candidate in PERTURBATION_KEY_CANDIDATES:
        if candidate.lower() in lower_to_actual:
            return lower_to_actual[candidate.lower()]
    scored: list[tuple[int, str]] = []
    for col in obs.columns:
        values = obs[col].astype(str)
        nunique = values.nunique(dropna=True)
        if 1 < nunique < max(5000, len(values) * 0.8):
            score = int(values.str.contains("control|ctrl|non-target|non_target|nt", case=False, regex=True).any())
            scored.append((score, col))
    if scored:
        return sorted(scored, reverse=True)[0][1]
    raise SystemExit("Could not infer perturbation key; pass --perturbation-key.")


def infer_gene_symbols(adata, gene_symbol_key: str | None) -> np.ndarray:
    if gene_symbol_key and gene_symbol_key in adata.var.columns:
        return adata.var[gene_symbol_key].astype(str).to_numpy()
    for candidate in ["gene_symbol", "gene_name", "symbol", "feature_name", "genes"]:
        if candidate in adata.var.columns:
            values = adata.var[candidate].astype(str).to_numpy()
            if len(set(values)) == len(values):
                return values
    return np.asarray(adata.var_names.astype(str))


def dense(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def sample_indices(indices: np.ndarray, max_cells: int, rng: np.random.Generator) -> np.ndarray:
    if max_cells and len(indices) > max_cells:
        return np.sort(rng.choice(indices, size=max_cells, replace=False))
    return indices


def select_gene_indices(X, genes: np.ndarray, max_genes: int | None) -> np.ndarray:
    valid = np.array([bool(g) and g.lower() not in {"nan", "none"} for g in genes])
    indices = np.flatnonzero(valid)
    if max_genes and len(indices) > max_genes:
        means = np.asarray(X[:, indices].mean(axis=0)).ravel()
        top = np.argsort(means)[::-1][:max_genes]
        indices = np.sort(indices[top])
    return indices


def mean_var(X) -> tuple[np.ndarray, np.ndarray]:
    arr = dense(X).astype(float)
    return np.nanmean(arr, axis=0), np.nanvar(arr, axis=0, ddof=1)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    idx = np.where(valid)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1.0)
    return q


def clean_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.+-]+", "_", str(value).strip())
    return cleaned.strip("_") or "perturbation"


def coexpression_with_target(control_matrix: np.ndarray, genes: np.ndarray, target: str) -> np.ndarray:
    corr = signed_coexpression_with_target(control_matrix, genes, target)
    return np.abs(corr)


def signed_coexpression_with_target(control_matrix: np.ndarray, genes: np.ndarray, target: str) -> np.ndarray:
    target_parts = re.split(r"[+|,;]", target)
    target_parts = [part.strip() for part in target_parts if part.strip()]
    target_indices: list[int] = []
    for part in target_parts:
        matches = np.where(genes == part)[0]
        if len(matches):
            target_indices.append(int(matches[0]))
    if not target_indices or control_matrix.shape[0] < 3:
        return np.full(len(genes), np.nan)
    arr = control_matrix.astype(float)
    gene_sd = np.nanstd(arr, axis=0, ddof=1)
    centered = arr - np.nanmean(arr, axis=0)
    corrs = []
    for target_idx in target_indices:
        target_vec = arr[:, target_idx]
        target_sd = np.nanstd(target_vec, ddof=1)
        if not np.isfinite(target_sd) or target_sd == 0:
            continue
        centered_target = target_vec - np.nanmean(target_vec)
        cov = np.nanmean(centered * centered_target[:, None], axis=0)
        denom = target_sd * gene_sd
        with np.errstate(divide="ignore", invalid="ignore"):
            corrs.append(cov / denom)
    if not corrs:
        return np.full(len(genes), np.nan)
    corr_mat = np.vstack(corrs)
    best = np.nanargmax(np.abs(np.nan_to_num(corr_mat, nan=0.0)), axis=0)
    return corr_mat[best, np.arange(corr_mat.shape[1])]


def make_baseline_scores(
    dataset_id: str,
    perturbation_id: str,
    genes: np.ndarray,
    scores: np.ndarray,
    model_id: str,
    model_name: str,
    model_version: str,
    score_type: str,
    source_file: str,
) -> pd.DataFrame:
    clean_scores = np.nan_to_num(scores.astype(float), nan=0.0, posinf=0.0, neginf=0.0)
    df = pd.DataFrame(
        {
            "dataset_id": dataset_id,
            "perturbation_id": perturbation_id,
            "model_id": model_id,
            "model_name": model_name,
            "model_family": "statistical_baseline",
            "model_version": model_version,
            "gene_id": np.nan,
            "gene_symbol": genes,
            "score": clean_scores,
            "score_abs": np.abs(clean_scores),
            "direction": np.where(clean_scores > 0, "positive", np.where(clean_scores < 0, "negative", "zero")),
            "scope_included": True,
            "score_type": score_type,
            "source_file": source_file,
            "adapter_version": "carf_public_perturbseq_v1",
        }
    )
    df = df.sort_values(["score_abs", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1, dtype=int)
    columns = [
        "dataset_id",
        "perturbation_id",
        "model_id",
        "model_name",
        "model_family",
        "model_version",
        "gene_id",
        "gene_symbol",
        "score",
        "score_abs",
        "rank",
        "direction",
        "scope_included",
        "score_type",
        "source_file",
        "adapter_version",
    ]
    return df[columns]


def prepare_dataset(args: argparse.Namespace) -> None:
    source_rows = select_sources(args.dataset_id)
    source = source_rows.iloc[0].to_dict()
    control_regex = args.control_regex or optional_string(source.get("default_control_regex")) or "control|ctrl|non-target|nt"
    rng = np.random.default_rng(args.seed)

    adata = read_h5ad(args.h5ad)
    perturbation_key = args.perturbation_key or optional_string(source.get("default_perturbation_key"))
    if not perturbation_key:
        perturbation_key = infer_perturbation_key(adata.obs)
    if perturbation_key not in adata.obs.columns:
        raise SystemExit(f"Perturbation key {perturbation_key!r} not found in obs.")

    labels = adata.obs[perturbation_key].astype(str)
    control_mask = labels.str.contains(control_regex, case=False, regex=True, na=False).to_numpy()
    if control_mask.sum() < args.min_cells:
        raise SystemExit(
            f"Only {control_mask.sum()} control cells detected with regex {control_regex!r}."
        )

    genes_all = infer_gene_symbols(adata, args.gene_symbol_key)
    gene_idx = select_gene_indices(adata.X, genes_all, args.max_genes)
    genes = genes_all[gene_idx]

    counts = labels[~control_mask].value_counts()
    counts = counts[counts >= args.min_cells]
    if args.perturbations:
        wanted = {item.strip() for item in args.perturbations.split(",") if item.strip()}
        perturbations = [p for p in counts.index if p in wanted]
    else:
        perturbations = counts.head(args.max_perturbations).index.tolist()
    if not perturbations:
        raise SystemExit("No perturbations met the cell-count filter.")

    run_dir = args.output_run_dir or (DEFAULT_RUNS_DIR / args.dataset_id)
    input_dir = run_dir / "inputs"
    native_dir = run_dir / "native"
    input_dir.mkdir(parents=True, exist_ok=True)
    native_dir.mkdir(parents=True, exist_ok=True)

    control_idx = sample_indices(np.flatnonzero(control_mask), args.max_cells_per_group, rng)
    control_X = dense(adata.X[control_idx, :][:, gene_idx]).astype(float)
    ctrl_mean = np.nanmean(control_X, axis=0)
    ctrl_var = np.nanvar(control_X, axis=0, ddof=1)
    ctrl_sd = np.sqrt(ctrl_var)

    all_ground_truth: list[pd.DataFrame] = []
    all_covariates: list[pd.DataFrame] = []
    all_method_scores: list[pd.DataFrame] = []

    for perturbation in perturbations:
        pert_mask = (labels == perturbation).to_numpy()
        pert_idx = sample_indices(np.flatnonzero(pert_mask), args.max_cells_per_group, rng)
        pert_X = dense(adata.X[pert_idx, :][:, gene_idx]).astype(float)
        pert_mean, pert_var = mean_var(pert_X)

        delta = pert_mean - ctrl_mean
        se = np.sqrt((pert_var / max(len(pert_idx), 1)) + (ctrl_var / max(len(control_idx), 1)))
        with np.errstate(divide="ignore", invalid="ignore"):
            z = delta / se
        p = 2 * stats.norm.sf(np.abs(z))
        fdr = bh_fdr(p)
        is_positive = (fdr <= args.fdr_alpha) & (np.abs(delta) >= args.min_abs_delta)
        perturbation_id = clean_id(perturbation)

        native_delta = pd.DataFrame(
            {
                "gene_symbol": genes,
                "observed_delta": delta,
                "z_score": z,
                "p_value": p,
                "fdr": fdr,
            }
        )
        native_path = native_dir / f"{perturbation_id}_observed_delta.csv"
        native_delta.to_csv(native_path, index=False)

        all_ground_truth.append(
            pd.DataFrame(
                {
                    "dataset_id": args.dataset_id,
                    "perturbation_id": perturbation_id,
                    "gene_id": np.nan,
                    "gene_symbol": genes,
                    "is_positive": is_positive,
                    "label_source": f"{args.dataset_id}_within_dataset_differential_expression",
                    "evidence_type": "normal_approximation_fdr",
                    "effect_direction": np.where(delta > 0, "positive", np.where(delta < 0, "negative", "zero")),
                    "fdr": fdr,
                    "split": "validation",
                    "notes": f"Derived from {len(pert_idx)} perturbation cells and {len(control_idx)} sampled control cells.",
                }
            )
        )

        abs_coexpr = coexpression_with_target(control_X, genes, perturbation)
        all_covariates.append(
            pd.DataFrame(
                {
                    "dataset_id": args.dataset_id,
                    "perturbation_id": perturbation_id,
                    "gene_id": np.nan,
                    "gene_symbol": genes,
                    "mean_expression": ctrl_mean,
                    "sd_expression": ctrl_sd,
                    "cv_expression": np.divide(
                        ctrl_sd,
                        ctrl_mean,
                        out=np.full_like(ctrl_sd, np.nan),
                        where=ctrl_mean != 0,
                    ),
                    "abs_coexpression_with_perturbed_gene": abs_coexpr,
                    "gene_family": "not_annotated",
                    "platform_gene_in_scope": True,
                }
            )
        )

        signed_coexpr = signed_coexpression_with_target(control_X, genes, perturbation)
        all_method_scores.extend(
            [
                make_baseline_scores(
                    args.dataset_id,
                    perturbation_id,
                    genes,
                    np.zeros(len(genes)),
                    "mean_baseline",
                    "Mean prediction baseline",
                    "control_mean_delta_zero",
                    "predicted_delta",
                    "generated_control_mean_delta_zero",
                ),
                make_baseline_scores(
                    args.dataset_id,
                    perturbation_id,
                    genes,
                    ctrl_mean,
                    "mean_expression_baseline",
                    "Mean expression baseline",
                    "control_mean_expression",
                    "mean_expression",
                    "generated_control_mean_expression",
                ),
                make_baseline_scores(
                    args.dataset_id,
                    perturbation_id,
                    genes,
                    signed_coexpr,
                    "control_coexpression_baseline",
                    "Control co-expression baseline",
                    "control_target_coexpression",
                    "correlation",
                    "generated_control_target_coexpression",
                ),
            ]
        )

        print(
            f"Prepared {args.dataset_id}/{perturbation_id}: "
            f"{len(pert_idx)} perturbed cells, {int(is_positive.sum())} positives"
        )

    ground_truth = pd.concat(all_ground_truth, ignore_index=True)
    covariates = pd.concat(all_covariates, ignore_index=True)
    method_scores = pd.concat(all_method_scores, ignore_index=True)

    write_csv(ground_truth, input_dir / "ground_truth.csv")
    write_csv(covariates, input_dir / "covariates.csv")
    write_csv(method_scores, input_dir / "method_scores.csv")

    run_manifest = {
        "benchmark_id": "carf_benchmark_v1",
        "benchmark_version": "1.1.0",
        "dataset_id": args.dataset_id,
        "perturbation_id": "multiple",
        "created_by": "CARF public Perturb-seq activation script",
        "inputs": {
            "method_scores": str(input_dir / "method_scores.csv"),
            "ground_truth": str(input_dir / "ground_truth.csv"),
            "covariates": str(input_dir / "covariates.csv"),
        },
        "source": source,
        "preparation": {
            "h5ad": str(args.h5ad),
            "perturbation_key": perturbation_key,
            "control_regex": control_regex,
            "max_perturbations": args.max_perturbations,
            "min_cells": args.min_cells,
            "max_cells_per_group": args.max_cells_per_group,
            "max_genes": args.max_genes,
            "fdr_alpha": args.fdr_alpha,
        },
        "metrics": {
            "k_values": [10, 25, 50, 75, 100, 200, 500, 1000],
            "pce_alpha": 0.5,
            "threshold_status": "heuristic_not_pass_fail",
        },
    }
    (run_dir / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote run manifest: {run_dir / 'run_manifest.json'}")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"Wrote {len(df):,} rows: {path}")


def main() -> int:
    args = parse_args()
    if args.command == "list-sources":
        list_sources()
    elif args.command == "materialize-sources":
        materialize_sources(args.dataset_id, args.runs_dir)
    elif args.command == "download":
        download_sources(args.dataset_id, args.raw_dir, args.decompress)
    elif args.command == "inspect":
        inspect_h5ad(args.h5ad)
    elif args.command == "prepare":
        prepare_dataset(args)
    else:
        raise SystemExit(f"Unsupported command: {args.command}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
