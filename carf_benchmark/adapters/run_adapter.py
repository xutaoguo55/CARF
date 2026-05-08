#!/usr/bin/env python3
"""Convert native perturbation-model outputs into CARF method_scores.csv.

This adapter layer intentionally separates "running the model" from
"standardizing the model output". Foundation-model repositories change quickly,
but their benchmark outputs usually reduce to one of two native artifacts:

1. a gene-level score table, or
2. a predicted post-perturbation expression matrix.

Both are converted here into the CARF method_scores schema without fabricating
scores when native outputs are unavailable.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


MODEL_SPECS: dict[str, dict[str, str]] = {
    "scgpt": {
        "model_name": "scGPT",
        "model_family": "foundation_model",
        "score_type": "predicted_delta",
    },
    "scfoundation": {
        "model_name": "scFoundation",
        "model_family": "foundation_model",
        "score_type": "predicted_delta",
    },
    "uce": {
        "model_name": "UCE",
        "model_family": "foundation_model",
        "score_type": "predicted_delta",
    },
    "scbert": {
        "model_name": "scBERT",
        "model_family": "foundation_model",
        "score_type": "predicted_delta",
    },
    "gears": {
        "model_name": "GEARS",
        "model_family": "deep_perturbation_model",
        "score_type": "predicted_delta",
    },
    "cpa": {
        "model_name": "CPA",
        "model_family": "deep_perturbation_model",
        "score_type": "predicted_delta",
    },
    "mean_baseline": {
        "model_name": "Mean prediction baseline",
        "model_family": "statistical_baseline",
        "score_type": "predicted_delta",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert native model outputs into CARF method_scores.csv."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    score_csv = subparsers.add_parser(
        "gene-score-csv",
        help="Convert a native gene-level score table into CARF method_scores.",
    )
    add_common_args(score_csv)
    score_csv.add_argument("--native-output", required=True, type=Path)
    score_csv.add_argument("--gene-column", default="gene_symbol")
    score_csv.add_argument("--score-column", default="predicted_delta")
    score_csv.add_argument("--gene-id-column")
    score_csv.add_argument("--direction-column")

    expr = subparsers.add_parser(
        "expression-delta",
        help=(
            "Convert predicted and baseline expression matrices into gene-level "
            "predicted deltas. Inputs may be CSV or h5ad."
        ),
    )
    add_common_args(expr)
    expr.add_argument("--predicted-expression", required=True, type=Path)
    expr.add_argument("--baseline-expression", required=True, type=Path)
    expr.add_argument(
        "--layer",
        default=None,
        help="Optional AnnData layer to read from h5ad files instead of X.",
    )

    append = subparsers.add_parser(
        "append",
        help="Append one or more method_scores files into a run input file.",
    )
    append.add_argument("--run-method-scores", required=True, type=Path)
    append.add_argument("--adapter-output", nargs="+", required=True, type=Path)
    append.add_argument("--output", required=True, type=Path)

    return parser.parse_args()


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--model-id", required=True, choices=sorted(MODEL_SPECS))
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--perturbation-id", required=True)
    parser.add_argument("--model-version", default="external_native_output")
    parser.add_argument("--adapter-version", default="carf_adapter_v1")
    parser.add_argument("--source-file")
    parser.add_argument("--output", required=True, type=Path)


def read_gene_score_csv(path: Path, gene_column: str, score_column: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [col for col in (gene_column, score_column) if col not in df.columns]
    if missing:
        raise SystemExit(
            f"{path} is missing required native columns: {', '.join(missing)}"
        )
    out = pd.DataFrame(
        {
            "gene_symbol": df[gene_column].astype(str),
            "score": pd.to_numeric(df[score_column], errors="coerce"),
        }
    )
    return out.join(df.drop(columns=[gene_column, score_column], errors="ignore"))


def read_expression_matrix(path: Path, layer: str | None = None) -> tuple[np.ndarray, list[str]]:
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".h5ad"):
        try:
            import anndata as ad
        except ImportError as exc:
            raise SystemExit("Reading h5ad native outputs requires anndata.") from exc
        adata = ad.read_h5ad(path)
        matrix = adata.layers[layer] if layer else adata.X
        genes = [str(g) for g in adata.var_names]
        return dense_array(matrix), genes

    df = pd.read_csv(path)
    if "gene_symbol" in df.columns and "predicted_delta" in df.columns:
        matrix = np.asarray([pd.to_numeric(df["predicted_delta"], errors="coerce")])
        return matrix, df["gene_symbol"].astype(str).tolist()
    first = df.columns[0]
    if first.lower() in {"cell", "cell_id", "obs_name", "index"}:
        df = df.drop(columns=[first])
    genes = [str(col) for col in df.columns]
    matrix = df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    return matrix, genes


def dense_array(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def expression_delta(predicted_path: Path, baseline_path: Path, layer: str | None) -> pd.DataFrame:
    predicted, predicted_genes = read_expression_matrix(predicted_path, layer)
    baseline, baseline_genes = read_expression_matrix(baseline_path, layer)

    shared = [gene for gene in predicted_genes if gene in set(baseline_genes)]
    if not shared:
        raise SystemExit("Predicted and baseline expression matrices share no genes.")

    pred_index = {gene: i for i, gene in enumerate(predicted_genes)}
    base_index = {gene: i for i, gene in enumerate(baseline_genes)}
    pred_mean = np.nanmean(predicted[:, [pred_index[g] for g in shared]], axis=0)
    base_mean = np.nanmean(baseline[:, [base_index[g] for g in shared]], axis=0)
    return pd.DataFrame({"gene_symbol": shared, "score": pred_mean - base_mean})


def standardize_scores(
    score_df: pd.DataFrame,
    args: argparse.Namespace,
    source_file: str,
) -> pd.DataFrame:
    spec = MODEL_SPECS[args.model_id]
    df = score_df.copy()
    df = df.dropna(subset=["gene_symbol", "score"])
    df = df[df["gene_symbol"].astype(str).str.len() > 0]
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df[np.isfinite(df["score"])]
    df["score_abs"] = df["score"].abs()
    df = df.sort_values(["score_abs", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1, dtype=int)

    gene_id = None
    if getattr(args, "gene_id_column", None) and args.gene_id_column in df.columns:
        gene_id = df[args.gene_id_column].astype(str)

    direction = np.where(df["score"] > 0, "positive", np.where(df["score"] < 0, "negative", "zero"))
    if getattr(args, "direction_column", None) and args.direction_column in df.columns:
        direction = df[args.direction_column].astype(str).to_numpy()

    out = pd.DataFrame(
        {
            "dataset_id": args.dataset_id,
            "perturbation_id": args.perturbation_id,
            "model_id": args.model_id,
            "model_name": spec["model_name"],
            "model_family": spec["model_family"],
            "model_version": args.model_version,
            "gene_id": gene_id if gene_id is not None else np.nan,
            "gene_symbol": df["gene_symbol"].astype(str),
            "score": df["score"].astype(float),
            "score_abs": df["score_abs"].astype(float),
            "rank": df["rank"].astype(int),
            "direction": direction,
            "scope_included": True,
            "score_type": spec["score_type"],
            "source_file": source_file,
            "adapter_version": args.adapter_version,
        }
    )
    return out


def append_method_scores(run_method_scores: Path, adapter_outputs: Iterable[Path], output: Path) -> None:
    frames = [pd.read_csv(run_method_scores)]
    frames.extend(pd.read_csv(path) for path in adapter_outputs)
    merged = pd.concat(frames, ignore_index=True)
    key_cols = ["dataset_id", "perturbation_id", "model_id", "gene_symbol"]
    duplicated = merged.duplicated(key_cols, keep=False)
    if duplicated.any():
        dupes = merged.loc[duplicated, key_cols].head(10).to_dict(orient="records")
        raise SystemExit(f"Duplicate CARF method-score keys after append: {dupes}")
    write_csv(merged, output)


def write_csv(df: pd.DataFrame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"Wrote {len(df):,} rows: {output}")


def main() -> int:
    args = parse_args()
    if args.command == "gene-score-csv":
        source_file = args.source_file or str(args.native_output)
        native = read_gene_score_csv(args.native_output, args.gene_column, args.score_column)
        write_csv(standardize_scores(native, args, source_file), args.output)
        return 0

    if args.command == "expression-delta":
        source_file = args.source_file or str(args.predicted_expression)
        native = expression_delta(args.predicted_expression, args.baseline_expression, args.layer)
        write_csv(standardize_scores(native, args, source_file), args.output)
        return 0

    if args.command == "append":
        append_method_scores(args.run_method_scores, args.adapter_output, args.output)
        return 0

    raise SystemExit(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
