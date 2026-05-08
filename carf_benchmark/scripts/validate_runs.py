#!/usr/bin/env python3
"""Validate all complete CARF-Benchmark runs listed in benchmark_v1.json."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
CARF_ROOT = SCRIPT_DIR.parents[0]
PROJECT_ROOT = CARF_ROOT.parent
sys.path.insert(0, str(SCRIPT_DIR))

from validate_schema import load_schema, validate_csv  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate complete CARF run inputs.")
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "carf_benchmark" / "configs" / "benchmark_v1.json",
    )
    parser.add_argument("--max-errors", type=int, default=20)
    parser.add_argument(
        "--require-all-inputs",
        action="store_true",
        help="Fail if any configured dataset lacks complete inputs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    schema_paths = {
        "method_scores": PROJECT_ROOT / config["schema"]["method_scores"],
        "ground_truth": PROJECT_ROOT / config["schema"]["ground_truth"],
        "covariates": PROJECT_ROOT / config["schema"]["covariates"],
    }
    schemas = {name: load_schema(path) for name, path in schema_paths.items()}
    any_errors = False
    complete = 0
    skipped = 0

    for dataset in config["datasets"]:
        run_dir = PROJECT_ROOT / dataset["run_dir"]
        input_dir = run_dir / "inputs"
        csv_paths = {
            "method_scores": input_dir / "method_scores.csv",
            "ground_truth": input_dir / "ground_truth.csv",
            "covariates": input_dir / "covariates.csv",
        }
        missing = [name for name, path in csv_paths.items() if not path.exists()]
        if missing:
            skipped += 1
            print(f"SKIP {dataset['dataset_id']}: missing {', '.join(missing)}")
            if args.require_all_inputs:
                any_errors = True
            continue

        complete += 1
        for name, path in csv_paths.items():
            errors = validate_csv(schemas[name], path, args.max_errors)
            if errors:
                any_errors = True
                print(f"FAIL {dataset['dataset_id']} {name}:")
                for error in errors:
                    print(f"  - {error}")
            else:
                print(f"PASS {dataset['dataset_id']} {name}: {path}")

    print(f"Validated complete runs: {complete}; skipped pending runs: {skipped}")
    return 1 if any_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
