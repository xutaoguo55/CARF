#!/usr/bin/env python3
"""Smoke tests for CARF model adapters."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ADAPTER = ROOT / "adapters" / "run_adapter.py"
FIXTURES = ROOT / "tests" / "fixtures" / "native_outputs"
TMP = ROOT / "tests" / "tmp"
REQUIRED_COLUMNS = {
    "dataset_id",
    "perturbation_id",
    "model_id",
    "model_name",
    "gene_symbol",
    "score",
    "score_abs",
    "rank",
    "scope_included",
    "score_type",
    "source_file",
}


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def assert_method_scores(path: Path, model_id: str) -> None:
    df = pd.read_csv(path)
    missing = REQUIRED_COLUMNS - set(df.columns)
    assert not missing, f"{path} missing {missing}"
    assert len(df) == 4
    assert set(df["model_id"]) == {model_id}
    assert list(df["rank"]) == [1, 2, 3, 4]
    assert df["score_abs"].is_monotonic_decreasing


def main() -> int:
    TMP.mkdir(parents=True, exist_ok=True)
    gene_score_models = ["scgpt", "scfoundation", "uce", "scbert"]
    expression_models = ["gears", "cpa"]

    for model_id in gene_score_models:
        out = TMP / f"{model_id}_method_scores.csv"
        run(
            [
                sys.executable,
                str(ADAPTER),
                "gene-score-csv",
                "--model-id",
                model_id,
                "--dataset-id",
                "toy_dataset",
                "--perturbation-id",
                "WWOX",
                "--native-output",
                str(FIXTURES / "gene_score_native.csv"),
                "--output",
                str(out),
            ]
        )
        assert_method_scores(out, model_id)

    for model_id in expression_models:
        out = TMP / f"{model_id}_method_scores.csv"
        run(
            [
                sys.executable,
                str(ADAPTER),
                "expression-delta",
                "--model-id",
                model_id,
                "--dataset-id",
                "toy_dataset",
                "--perturbation-id",
                "WWOX",
                "--predicted-expression",
                str(FIXTURES / "predicted_expression.csv"),
                "--baseline-expression",
                str(FIXTURES / "baseline_expression.csv"),
                "--output",
                str(out),
            ]
        )
        assert_method_scores(out, model_id)

    print("Adapter smoke tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
