#!/usr/bin/env python3
"""Validate CARF-Benchmark CSV files against lightweight schema metadata."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


TRUE_VALUES = {"true", "t", "1", "yes", "y"}
FALSE_VALUES = {"false", "f", "0", "no", "n"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate a CARF-Benchmark CSV file.")
    parser.add_argument("--schema", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--max-errors", type=int, default=20)
    return parser.parse_args()


def load_schema(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def coerce(value: str, expected_type: str) -> bool:
    if value is None or value == "" or value.lower() in {"na", "nan", "null"}:
        return True
    try:
        if expected_type == "string":
            return True
        if expected_type == "number":
            number = float(value)
            return math.isfinite(number)
        if expected_type == "integer":
            return float(value).is_integer()
        if expected_type == "boolean":
            return value.lower() in TRUE_VALUES | FALSE_VALUES
    except ValueError:
        return False
    return True


def validate_csv(schema: dict, csv_path: Path, max_errors: int) -> list[str]:
    errors: list[str] = []
    required = schema.get("required_columns", [])
    column_types = schema.get("column_types", {})
    primary_key = schema.get("primary_key", [])

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames or []
        missing = [col for col in required if col not in columns]
        if missing:
            errors.append(f"Missing required columns: {', '.join(missing)}")
            return errors

        seen_keys: set[tuple[str, ...]] = set()
        for row_number, row in enumerate(reader, start=2):
            if primary_key:
                key = tuple(row.get(col, "") for col in primary_key)
                if key in seen_keys:
                    errors.append(f"Row {row_number}: duplicate primary key {key}")
                seen_keys.add(key)

            for col, expected_type in column_types.items():
                if col in row and not coerce(row[col], expected_type):
                    errors.append(
                        f"Row {row_number}: column {col} value {row[col]!r} is not {expected_type}"
                    )
            if len(errors) >= max_errors:
                errors.append(f"Stopped after {max_errors} errors")
                break

    return errors


def main() -> int:
    args = parse_args()
    schema = load_schema(args.schema)
    errors = validate_csv(schema, args.csv, args.max_errors)
    if errors:
        print(f"Schema validation failed for {args.csv}:")
        for error in errors:
            print(f"  - {error}")
        return 1
    print(f"Schema validation passed: {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
