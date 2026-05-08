#!/usr/bin/env python3
"""Raw Geneformer embedding-density and attention diagnostics.

This script uses raw model hidden states and attention weights from a
tokenized Geneformer dataset. It replaces the earlier expression-space proxy
with a direct diagnostic of the model representation used for inference.

Full 3,992-token attention tensors are O(sequence_length^2) and are too large
to retain across all layers/cells. Therefore the default workflow computes
embedding density on full sequences and computes raw attention in a fixed
rank-order window around the WWOX token. The window size and sampled cells are
reported in the metadata and can be increased when memory permits.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch
from datasets import load_from_disk
from transformers import AutoModel

try:
    from scipy.stats import spearmanr
except Exception:  # pragma: no cover - scipy is expected but not required
    spearmanr = None


SPECIAL_TOKEN_IDS = {0, 1, 2, 3}
WWOX_ENSEMBL = "ENSG00000186153"


def first_existing(paths: Iterable[str | Path]) -> Path | None:
    for path in paths:
        path = Path(path).expanduser()
        if path.exists():
            return path
    return None


def default_project_dir() -> Path:
    return Path(__file__).resolve().parents[1]


def is_lfs_pointer(path: Path) -> bool:
    if not path.exists() or not path.is_file():
        return False
    try:
        head = path.read_bytes()[:128]
    except OSError:
        return False
    return head.startswith(b"version https://git-lfs.github.com/spec")


def default_model_dir() -> Path | None:
    env = os.environ.get("GENEFORMER_MODEL_DIR")
    candidates = [env] if env else []
    candidates.extend([
        "/private/tmp/Geneformer",
        "/tmp/Geneformer",
        "/private/tmp/Geneformer/Geneformer-V2-104M",
        "/tmp/Geneformer/Geneformer-V2-104M",
    ])
    return first_existing([p for p in candidates if p])


def default_dataset_path() -> Path | None:
    env = os.environ.get("GENEFORMER_DATASET_DIR")
    candidates = [env] if env else []
    candidates.extend([
        "/private/tmp/geneformer_test/tokenized_final/GSE10846_filtered.dataset",
        "/tmp/geneformer_test/tokenized_final/GSE10846_filtered.dataset",
        "/private/tmp/geneformer_test/tokenized_filtered/GSE10846_filtered.dataset",
        "/tmp/geneformer_test/tokenized_filtered/GSE10846_filtered.dataset",
    ])
    return first_existing([p for p in candidates if p])


def default_token_dictionary(model_dir: Path | None) -> Path | None:
    env = os.environ.get("GENEFORMER_TOKEN_DICTIONARY")
    candidates = [env] if env else []
    if model_dir is not None:
        candidates.extend([
            model_dir / "geneformer" / "token_dictionary_gc104M.pkl",
            model_dir / "geneformer" / "gene_dictionaries_30m" / "token_dictionary_gc30M.pkl",
        ])
    candidates.extend([
        "/private/tmp/Geneformer/geneformer/token_dictionary_gc104M.pkl",
        "/tmp/Geneformer/geneformer/token_dictionary_gc104M.pkl",
    ])
    return first_existing([p for p in candidates if p])


def default_gene_name_dictionary(model_dir: Path | None) -> Path | None:
    env = os.environ.get("GENEFORMER_GENE_NAME_DICTIONARY")
    candidates = [env] if env else []
    if model_dir is not None:
        candidates.extend([
            model_dir / "geneformer" / "gene_name_id_dict_gc104M.pkl",
            model_dir / "geneformer" / "gene_dictionaries_30m" / "gene_name_id_dict_gc30M.pkl",
        ])
    candidates.extend([
        "/private/tmp/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
        "/tmp/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
    ])
    return first_existing([p for p in candidates if p])


def parse_args() -> argparse.Namespace:
    project_dir = default_project_dir()
    model_dir = default_model_dir()
    dataset_path = default_dataset_path()
    token_dictionary = default_token_dictionary(model_dir)
    gene_name_dictionary = default_gene_name_dictionary(model_dir)

    parser = argparse.ArgumentParser(
        description="Compute raw Geneformer hidden-state density and WWOX attention diagnostics."
    )
    parser.add_argument("--project-dir", type=Path, default=project_dir)
    parser.add_argument("--model-dir", type=Path, default=model_dir)
    parser.add_argument("--dataset-path", type=Path, default=dataset_path)
    parser.add_argument("--token-dictionary", type=Path, default=token_dictionary)
    parser.add_argument("--gene-name-dictionary", type=Path, default=gene_name_dictionary)
    parser.add_argument("--wwox-ensembl", default=WWOX_ENSEMBL)
    parser.add_argument("--max-cells", type=int, default=8)
    parser.add_argument("--attention-cells", type=int, default=4)
    parser.add_argument("--attention-window", type=int, default=512)
    parser.add_argument("--density-k", type=int, default=20)
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda", "mps", "auto"])
    parser.add_argument("--preflight", action="store_true")
    return parser.parse_args()


def load_pickle(path: Path) -> dict:
    with path.open("rb") as handle:
        return pickle.load(handle)


def resolve_device(requested: str) -> torch.device:
    if requested == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is not available.")
    if requested == "mps":
        if getattr(torch.backends, "mps", None) is None or not torch.backends.mps.is_available():
            raise RuntimeError("MPS was requested but is not available.")
    return torch.device(requested)


def check_inputs(args: argparse.Namespace) -> tuple[bool, list[str]]:
    errors: list[str] = []
    for label, path in [
        ("project_dir", args.project_dir),
        ("model_dir", args.model_dir),
        ("dataset_path", args.dataset_path),
        ("token_dictionary", args.token_dictionary),
        ("gene_name_dictionary", args.gene_name_dictionary),
    ]:
        if path is None:
            errors.append(f"{label} was not provided and no default was found.")
        elif not Path(path).exists():
            errors.append(f"{label} does not exist: {path}")

    if args.model_dir is not None and Path(args.model_dir).exists():
        model_files = [
            Path(args.model_dir) / "model.safetensors",
            Path(args.model_dir) / "pytorch_model.bin",
        ]
        usable_weights = [p for p in model_files if p.exists() and not is_lfs_pointer(p)]
        if len(usable_weights) == 0:
            errors.append(
                f"model_dir has no usable local weights (Git LFS pointer files do not count): {args.model_dir}"
            )

    if args.max_cells <= 0:
        errors.append("--max-cells must be positive.")
    if args.attention_cells < 0:
        errors.append("--attention-cells cannot be negative.")
    if args.attention_window <= 8:
        errors.append("--attention-window must be > 8.")
    if args.density_k <= 0:
        errors.append("--density-k must be positive.")

    return len(errors) == 0, errors


def load_model(model_dir: Path, device: torch.device, eager_attention: bool = False) -> AutoModel:
    kwargs = {"trust_remote_code": True, "add_pooling_layer": False}
    if eager_attention:
        kwargs["attn_implementation"] = "eager"
    try:
        model = AutoModel.from_pretrained(model_dir, **kwargs)
    except TypeError:
        kwargs.pop("attn_implementation", None)
        model = AutoModel.from_pretrained(model_dir, **kwargs)
    model.eval()
    model.to(device)
    return model


def select_cell_indices(n_rows: int, n_cells: int) -> list[int]:
    n_cells = min(n_rows, n_cells)
    if n_cells <= 0:
        return []
    return sorted(set(np.linspace(0, n_rows - 1, n_cells, dtype=int).tolist()))


def find_token_positions(input_ids: list[int], token_id: int) -> list[int]:
    return [i for i, value in enumerate(input_ids) if value == token_id]


def run_embedding_pass(
    model: AutoModel,
    dataset,
    cell_indices: list[int],
    device: torch.device,
) -> tuple[dict[int, np.ndarray], dict[int, int]]:
    token_sums: dict[int, np.ndarray] = {}
    token_counts: dict[int, int] = {}

    for run_idx, cell_idx in enumerate(cell_indices, start=1):
        input_ids = list(map(int, dataset[int(cell_idx)]["input_ids"]))
        ids_tensor = torch.tensor([input_ids], dtype=torch.long, device=device)
        mask = torch.ones_like(ids_tensor)

        start = time.time()
        with torch.no_grad():
            output = model(input_ids=ids_tensor, attention_mask=mask, return_dict=True)
        hidden = output.last_hidden_state[0].detach().cpu().numpy().astype(np.float32)

        for position, token_id in enumerate(input_ids):
            if token_id in SPECIAL_TOKEN_IDS:
                continue
            if token_id not in token_sums:
                token_sums[token_id] = np.zeros(hidden.shape[1], dtype=np.float64)
                token_counts[token_id] = 0
            token_sums[token_id] += hidden[position].astype(np.float64)
            token_counts[token_id] += 1

        elapsed = time.time() - start
        print(f"Embedding pass {run_idx}/{len(cell_indices)}: cell={cell_idx}, tokens={len(input_ids)}, {elapsed:.1f}s")

        del output, hidden, ids_tensor, mask
        gc.collect()

    token_means = {
        token_id: (token_sums[token_id] / token_counts[token_id]).astype(np.float32)
        for token_id in token_sums
    }
    return token_means, token_counts


def compute_density(
    token_means: dict[int, np.ndarray],
    token_counts: dict[int, int],
    density_k: int,
) -> pd.DataFrame:
    token_ids = np.array(sorted(token_means), dtype=np.int64)
    matrix = np.vstack([token_means[int(token_id)] for token_id in token_ids]).astype(np.float32)
    norms = np.linalg.norm(matrix, axis=1)
    safe_norms = np.where(norms == 0, 1, norms)
    normalized = matrix / safe_norms[:, None]

    similarity = normalized @ normalized.T
    np.fill_diagonal(similarity, -np.inf)
    k = min(density_k, max(1, similarity.shape[0] - 1))
    topk = np.partition(similarity, kth=similarity.shape[1] - k, axis=1)[:, -k:]
    local_density = topk.mean(axis=1)

    return pd.DataFrame({
        "token_id": token_ids,
        "raw_embedding_norm": norms,
        "raw_embedding_cells_observed": [token_counts[int(token_id)] for token_id in token_ids],
        f"raw_embedding_density_k{k}": local_density,
        f"raw_embedding_sparsity_k{k}": 1 - local_density,
    })


def attention_window(input_ids: list[int], target_position: int, window_size: int) -> tuple[list[int], int, int, int]:
    if window_size >= len(input_ids):
        return input_ids, target_position, 0, len(input_ids)
    start = max(0, min(target_position - window_size // 2, len(input_ids) - window_size))
    end = start + window_size
    return input_ids[start:end], target_position - start, start, end


def run_attention_pass(
    model: AutoModel,
    dataset,
    cell_indices: list[int],
    device: torch.device,
    wwox_token_id: int,
    window_size: int,
) -> tuple[pd.DataFrame, list[dict]]:
    records: dict[int, dict[str, float]] = {}
    windows: list[dict] = []

    def update(token_id: int, values: dict[str, float]) -> None:
        if token_id in SPECIAL_TOKEN_IDS:
            return
        rec = records.setdefault(token_id, {
            "token_id": token_id,
            "raw_attention_cells_observed": 0,
            "raw_attention_wwox_to_gene_all_layers": 0.0,
            "raw_attention_gene_to_wwox_all_layers": 0.0,
            "raw_attention_wwox_to_gene_last_layer": 0.0,
            "raw_attention_gene_to_wwox_last_layer": 0.0,
        })
        rec["raw_attention_cells_observed"] += 1
        for key, value in values.items():
            rec[key] += value

    for run_idx, cell_idx in enumerate(cell_indices, start=1):
        full_ids = list(map(int, dataset[int(cell_idx)]["input_ids"]))
        positions = find_token_positions(full_ids, wwox_token_id)
        if not positions:
            print(f"Attention pass {run_idx}/{len(cell_indices)}: cell={cell_idx}, WWOX token absent; skipped")
            continue

        window_ids, wwox_local_pos, start, end = attention_window(full_ids, positions[0], window_size)
        ids_tensor = torch.tensor([window_ids], dtype=torch.long, device=device)
        mask = torch.ones_like(ids_tensor)

        start_time = time.time()
        with torch.no_grad():
            output = model(
                input_ids=ids_tensor,
                attention_mask=mask,
                output_attentions=True,
                return_dict=True,
            )

        layer_from = []
        layer_to = []
        for attention in output.attentions:
            # Shape: batch x heads x query_token x key_token.
            mean_heads = attention[0].mean(dim=0)
            layer_from.append(mean_heads[wwox_local_pos, :].detach().cpu().numpy())
            layer_to.append(mean_heads[:, wwox_local_pos].detach().cpu().numpy())

        all_from = np.vstack(layer_from).mean(axis=0)
        all_to = np.vstack(layer_to).mean(axis=0)
        last_from = layer_from[-1]
        last_to = layer_to[-1]

        for local_pos, token_id in enumerate(window_ids):
            update(token_id, {
                "raw_attention_wwox_to_gene_all_layers": float(all_from[local_pos]),
                "raw_attention_gene_to_wwox_all_layers": float(all_to[local_pos]),
                "raw_attention_wwox_to_gene_last_layer": float(last_from[local_pos]),
                "raw_attention_gene_to_wwox_last_layer": float(last_to[local_pos]),
            })

        elapsed = time.time() - start_time
        print(
            f"Attention pass {run_idx}/{len(cell_indices)}: cell={cell_idx}, "
            f"window={start}:{end}, wwox_local_pos={wwox_local_pos}, {elapsed:.1f}s"
        )
        windows.append({
            "cell_index": int(cell_idx),
            "sequence_length": len(full_ids),
            "window_start": int(start),
            "window_end": int(end),
            "wwox_position_full": int(positions[0]),
            "wwox_position_window": int(wwox_local_pos),
        })

        del output, ids_tensor, mask
        gc.collect()

    rows = []
    for rec in records.values():
        count = rec["raw_attention_cells_observed"]
        row = {"token_id": rec["token_id"], "raw_attention_cells_observed": count}
        for key, value in rec.items():
            if key in {"token_id", "raw_attention_cells_observed"}:
                continue
            row[key] = value / count if count else np.nan
        rows.append(row)

    return pd.DataFrame(rows), windows


def build_gene_maps(project_dir: Path, token_dictionary: dict, gene_name_dictionary: dict) -> tuple[dict[int, str], dict[str, str], dict[str, int]]:
    token_to_ensembl = {
        int(token_id): str(ensembl)
        for ensembl, token_id in token_dictionary.items()
        if isinstance(token_id, (int, np.integer)) and str(ensembl).startswith("ENSG")
    }

    ensembl_to_symbol = {
        str(ensembl): str(symbol)
        for symbol, ensembl in gene_name_dictionary.items()
        if isinstance(symbol, str) and isinstance(ensembl, str)
    }

    symbol_to_ensembl = {
        str(symbol): str(ensembl)
        for symbol, ensembl in gene_name_dictionary.items()
        if isinstance(symbol, str) and isinstance(ensembl, str)
    }

    mapping_file = project_dir / "benchmark_results" / "gene_symbol_to_ensembl.csv"
    if mapping_file.exists():
        mapping = pd.read_csv(mapping_file)
        for _, row in mapping.dropna(subset=["SYMBOL", "ENSEMBL"]).iterrows():
            symbol = str(row["SYMBOL"])
            ensembl = str(row["ENSEMBL"])
            symbol_to_ensembl.setdefault(symbol, ensembl)
            ensembl_to_symbol.setdefault(ensembl, symbol)

    symbol_to_token = {
        symbol: int(token_dictionary[ensembl])
        for symbol, ensembl in symbol_to_ensembl.items()
        if ensembl in token_dictionary
    }
    return token_to_ensembl, ensembl_to_symbol, symbol_to_token


def enrich_results(
    project_dir: Path,
    density_df: pd.DataFrame,
    attention_df: pd.DataFrame,
    token_to_ensembl: dict[int, str],
    ensembl_to_symbol: dict[str, str],
    symbol_to_token: dict[str, int],
    density_k: int,
) -> pd.DataFrame:
    results = density_df.copy()
    if not attention_df.empty:
        results = results.merge(attention_df, on="token_id", how="left")

    results["ensembl_id"] = results["token_id"].map(token_to_ensembl)
    results["gene_symbol"] = results["ensembl_id"].map(ensembl_to_symbol)

    gf_file = project_dir / "benchmark_results" / "benchmark_geneformer_50cell.csv"
    gf = pd.read_csv(gf_file)
    gf = gf.rename(columns={
        "rank": "geneformer_rank",
        "method": "geneformer_method",
    })
    gf["geneformer_token_id_from_dictionary"] = gf["gene_symbol"].map(symbol_to_token)

    results = results.merge(gf, on="gene_symbol", how="outer")
    results["token_id"] = results["token_id"].fillna(results["geneformer_token_id_from_dictionary"])

    expr_file = project_dir / "benchmark_results" / "GSE10846_gene_expression_log2.csv"
    expr = pd.read_csv(expr_file, index_col=0)
    mean_expr = expr.mean(axis=1, skipna=True)
    results["mean_expr_log2"] = results["gene_symbol"].map(mean_expr)

    gt_file = project_dir / "benchmark_results" / "ground_truth_29.csv"
    if gt_file.exists():
        gt = pd.read_csv(gt_file)
        status = dict(zip(gt["gene"], gt["status"]))
        results["is_validated_target"] = results["gene_symbol"].map(status).eq(True)
    else:
        results["is_validated_target"] = False

    density_col = f"raw_embedding_density_k{min(density_k, max(1, len(density_df) - 1))}"
    results["raw_mapping_status"] = np.select(
        [
            results[density_col].notna() & results["geneformer_rank"].notna(),
            results[density_col].notna() & results["geneformer_rank"].isna(),
            results[density_col].isna() & results["geneformer_rank"].notna(),
        ],
        [
            "mapped_geneformer_gene",
            "raw_token_not_in_geneformer_csv",
            "geneformer_gene_without_raw_token_mapping",
        ],
        default="unmapped",
    )

    ordered_cols = [
        "gene_symbol", "ensembl_id", "token_id", "raw_mapping_status",
        "geneformer_rank", "cosine_shift", "abs_cosine_shift", "cos_mean", "cos_stdev",
        "mean_expr_log2", "is_validated_target",
        "raw_embedding_norm", "raw_embedding_cells_observed",
        density_col, density_col.replace("density", "sparsity"),
        "raw_attention_cells_observed",
        "raw_attention_wwox_to_gene_all_layers",
        "raw_attention_gene_to_wwox_all_layers",
        "raw_attention_wwox_to_gene_last_layer",
        "raw_attention_gene_to_wwox_last_layer",
    ]
    for col in ordered_cols:
        if col not in results.columns:
            results[col] = np.nan
    remaining = [col for col in results.columns if col not in ordered_cols]
    return results[ordered_cols + remaining].sort_values(
        by=["geneformer_rank", "gene_symbol"], na_position="last"
    )


def spearman_summary(label: str, df: pd.DataFrame, x_col: str, y_col: str) -> dict:
    valid = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) < 3:
        return {"metric": label, "rho": np.nan, "p_value": np.nan, "n": len(valid)}
    if spearmanr is None:
        return {
            "metric": label,
            "rho": valid[x_col].corr(valid[y_col], method="spearman"),
            "p_value": np.nan,
            "n": len(valid),
        }
    rho, p_value = spearmanr(valid[x_col], valid[y_col])
    return {"metric": label, "rho": float(rho), "p_value": float(p_value), "n": len(valid)}


def build_summary(results: pd.DataFrame, density_k: int, metadata: dict) -> pd.DataFrame:
    density_col = f"raw_embedding_density_k{min(density_k, max(1, results['raw_embedding_norm'].notna().sum() - 1))}"
    sparsity_col = density_col.replace("density", "sparsity")
    metrics = [
        spearman_summary("raw density vs |cosine shift|", results, density_col, "abs_cosine_shift"),
        spearman_summary("raw sparsity vs |cosine shift|", results, sparsity_col, "abs_cosine_shift"),
        spearman_summary("raw density vs mean expression", results, density_col, "mean_expr_log2"),
        spearman_summary("raw embedding norm vs |cosine shift|", results, "raw_embedding_norm", "abs_cosine_shift"),
        spearman_summary("WWOX->gene attention (all layers) vs |cosine shift|", results, "raw_attention_wwox_to_gene_all_layers", "abs_cosine_shift"),
        spearman_summary("gene->WWOX attention (all layers) vs |cosine shift|", results, "raw_attention_gene_to_wwox_all_layers", "abs_cosine_shift"),
        spearman_summary("WWOX->gene attention (last layer) vs |cosine shift|", results, "raw_attention_wwox_to_gene_last_layer", "abs_cosine_shift"),
        spearman_summary("gene->WWOX attention (last layer) vs |cosine shift|", results, "raw_attention_gene_to_wwox_last_layer", "abs_cosine_shift"),
        spearman_summary("WWOX->gene attention (all layers) vs mean expression", results, "raw_attention_wwox_to_gene_all_layers", "mean_expr_log2"),
        spearman_summary("gene->WWOX attention (all layers) vs mean expression", results, "raw_attention_gene_to_wwox_all_layers", "mean_expr_log2"),
    ]
    summary = pd.DataFrame(metrics)
    for key, value in metadata.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            summary[key] = value
    return summary


def main() -> int:
    args = parse_args()
    ok, errors = check_inputs(args)
    if args.preflight:
        print(json.dumps({
            "ok": ok,
            "errors": errors,
            "project_dir": str(args.project_dir),
            "model_dir": str(args.model_dir) if args.model_dir else None,
            "dataset_path": str(args.dataset_path) if args.dataset_path else None,
            "token_dictionary": str(args.token_dictionary) if args.token_dictionary else None,
            "gene_name_dictionary": str(args.gene_name_dictionary) if args.gene_name_dictionary else None,
        }, indent=2))
        return 0 if ok else 2
    if not ok:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2

    project_dir = args.project_dir.resolve()
    output_dir = project_dir / "benchmark_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    token_dictionary = load_pickle(args.token_dictionary)
    gene_name_dictionary = load_pickle(args.gene_name_dictionary)
    wwox_token_id = int(token_dictionary[args.wwox_ensembl])
    token_to_ensembl, ensembl_to_symbol, symbol_to_token = build_gene_maps(
        project_dir, token_dictionary, gene_name_dictionary
    )

    dataset = load_from_disk(str(args.dataset_path))
    device = resolve_device(args.device)
    embedding_indices = select_cell_indices(len(dataset), args.max_cells)
    attention_indices = select_cell_indices(len(dataset), min(args.attention_cells, args.max_cells))

    print("Raw Geneformer diagnostic configuration:")
    print(f"  model_dir: {args.model_dir}")
    print(f"  dataset_path: {args.dataset_path}")
    print(f"  rows: {len(dataset)}")
    print(f"  WWOX token: {wwox_token_id}")
    print(f"  device: {device}")
    print(f"  embedding cells: {embedding_indices}")
    print(f"  attention cells: {attention_indices}")

    model = load_model(args.model_dir, device, eager_attention=False)
    token_means, token_counts = run_embedding_pass(model, dataset, embedding_indices, device)
    density_df = compute_density(token_means, token_counts, args.density_k)
    del model
    gc.collect()
    if device.type == "mps" and hasattr(torch.mps, "empty_cache"):
        torch.mps.empty_cache()
    if device.type == "cuda":
        torch.cuda.empty_cache()

    attention_df = pd.DataFrame()
    attention_windows: list[dict] = []
    if len(attention_indices) > 0:
        attention_model = load_model(args.model_dir, device, eager_attention=True)
        attention_df, attention_windows = run_attention_pass(
            attention_model,
            dataset,
            attention_indices,
            device,
            wwox_token_id,
            args.attention_window,
        )
        del attention_model
        gc.collect()

    metadata = {
        "model_dir": str(args.model_dir),
        "dataset_path": str(args.dataset_path),
        "token_dictionary": str(args.token_dictionary),
        "gene_name_dictionary": str(args.gene_name_dictionary),
        "wwox_ensembl": args.wwox_ensembl,
        "wwox_token_id": wwox_token_id,
        "dataset_rows": len(dataset),
        "embedding_cells_requested": args.max_cells,
        "embedding_cells_used": len(embedding_indices),
        "attention_cells_requested": args.attention_cells,
        "attention_cells_used": len(attention_indices),
        "attention_window": args.attention_window,
        "density_k": args.density_k,
        "device": str(device),
        "n_raw_embedding_tokens": int(len(density_df)),
        "n_raw_attention_tokens": int(len(attention_df)),
    }

    results = enrich_results(
        project_dir,
        density_df,
        attention_df,
        token_to_ensembl,
        ensembl_to_symbol,
        symbol_to_token,
        args.density_k,
    )
    summary = build_summary(results, args.density_k, metadata)

    result_path = output_dir / "benchmark_raw_embedding_attention.csv"
    summary_path = output_dir / "benchmark_raw_embedding_attention_summary.csv"
    metadata_path = output_dir / "benchmark_raw_embedding_attention_metadata.json"

    results.to_csv(result_path, index=False)
    summary.to_csv(summary_path, index=False)
    metadata_with_windows = dict(metadata)
    metadata_with_windows["attention_windows"] = attention_windows
    metadata_path.write_text(json.dumps(metadata_with_windows, indent=2), encoding="utf-8")

    print(f"Saved gene-level raw diagnostics: {result_path}")
    print(f"Saved raw diagnostic summary: {summary_path}")
    print(f"Saved raw diagnostic metadata: {metadata_path}")
    print(summary[["metric", "rho", "p_value", "n"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
