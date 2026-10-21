"""Verify author-deposited RPE1 cell counts against the raw pseudobulk companion.

The command never opens the single-cell HDF5 file before its full-file source
receipt, fixed byte count, MD5, and SHA256 have passed. It does not download data,
impute missing label counts, normalize expression, or change either source file.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
SINGLECELL = HERE / "raw/rpe1_raw_singlecell_01.h5ad"
BULK = HERE / "raw/rpe1_raw_bulk_01.h5ad"
RECEIPTS = HERE / "provenance/download_receipts"
OUTPUT = HERE / "provenance/rpe1_contents_v1"
EXPECTED_BYTES = 8_700_873_216
EXPECTED_MD5 = "6a2a9d0d2bf4ec147f4d1104043b268c"
BULK_BYTES = 95_350_546
BULK_MD5 = "74765fa87635467a869ea972356ae0e7"
ATOL = 1e-5
RTOL = 1e-5


def file_hashes(path: Path) -> dict:
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(chunk)
            sha.update(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha.hexdigest()}


def source_gate(path: Path, receipt_path: Path, size: int, md5: str, file_id: int) -> dict:
    """Check receipt and raw file bytes without opening an HDF5 object."""
    if not path.exists() or path.stat().st_size != size:
        raise ValueError(f"Source incomplete or wrong size: {path}")
    if not receipt_path.exists():
        raise ValueError(f"Full-file verified download receipt is absent: {receipt_path}")
    receipt = json.loads(receipt_path.read_text())
    record = receipt.get("source_record", {})
    if (receipt.get("status") != "verified_author_deposit"
            or receipt.get("actual_bytes") != size
            or receipt.get("actual_md5") != md5
            or record.get("expected_md5") != md5
            or record.get("file_id") != file_id
            or record.get("file_name") != path.name
            or record.get("bytes", record.get("expected_bytes")) != size):
        raise ValueError(f"Download receipt does not satisfy the fixed source gate: {receipt_path}")
    expected_sha = receipt.get("actual_sha256", "")
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha):
        raise ValueError("Download receipt lacks a valid full-file SHA256")
    before = path.stat()
    print(f"Verifying complete file hashes: {path.name}", flush=True)
    hashes = file_hashes(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Source changed during full-file hashing")
    if hashes["md5"] != md5 or hashes["sha256"] != expected_sha:
        raise ValueError(f"Current source hashes do not match the full-file receipt: {path}")
    return {"path": str(path), "bytes": size, "mtime_ns": after.st_mtime_ns,
            "file_id": file_id, "receipt_sha256": file_hashes(receipt_path)["sha256"],
            **hashes}


def decode_text(values):
    a = np.asarray(values)
    if a.dtype.kind not in "OSU":
        return a
    return np.asarray([v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else v
                       for v in a], dtype=object)


def decode_categories(codes, categories):
    codes, categories = np.asarray(codes), decode_text(categories)
    if codes.ndim != 1 or categories.ndim != 1 or codes.dtype.kind not in "iu":
        raise ValueError("Categorical codes must be a one-dimensional integer array")
    if np.any(codes < -1) or np.any(codes >= len(categories)):
        raise ValueError("Out-of-range categorical code")
    result = np.empty(len(codes), dtype=object)
    result[:] = None
    valid = codes >= 0
    result[valid] = categories[codes[valid]]
    return result


def read_column(group, name):
    """Read plain, modern categorical, or legacy reference-attribute columns.

    Code -1 is preserved as None. In particular it must never index the last
    category. The legacy reader is local to this audit; old scripts are unchanged.
    """
    node = group[name]
    if isinstance(node, h5py.Group):
        if "codes" not in node or "categories" not in node:
            raise ValueError(f"Unsupported categorical group: {node.name}")
        return decode_categories(node["codes"][:], node["categories"][:])
    if node.ndim != 1:
        raise ValueError(f"Metadata column is not one-dimensional: {node.name}")
    if "categories" in node.attrs:
        reference = node.attrs["categories"]
        if isinstance(reference, h5py.Reference):
            if not reference:
                raise ValueError("Null HDF5 categories reference")
            categories = node.file[reference]
        elif isinstance(reference, (str, bytes)):
            key = reference.decode() if isinstance(reference, bytes) else reference
            categories = node.file[key] if key.startswith("/") else node.parent[key]
        else:
            raise ValueError(f"Unsupported categories attribute: {node.name}")
        return decode_categories(node[:], categories[:])
    return decode_text(node[:])


def label_or_none(value):
    if value is None or (isinstance(value, (float, np.floating)) and np.isnan(value)):
        return None
    if isinstance(value, (bytes, np.bytes_)):
        value = value.decode("utf-8")
    if not isinstance(value, str):
        raise ValueError(f"Expected a string identity, received {type(value).__name__}")
    return value if value.strip() else None


def strings(group, name):
    return np.asarray([label_or_none(v) for v in read_column(group, name)], dtype=object)


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def matrix_schema(f):
    result = {"roots": list(f), "X": {"shape": list(f["X"].shape),
              "dtype": str(f["X"].dtype), "chunks": f["X"].chunks,
              "compression": f["X"].compression}}
    for name in ("obs", "var"):
        result[name] = {key: {"kind": "group" if isinstance(node, h5py.Group) else "dataset",
                              "shape": list(node.shape) if isinstance(node, h5py.Dataset) else None,
                              "legacy_categories_reference": bool("categories" in node.attrs)}
                        for key, node in f[name].items()}
    return result


def audit_h5ad_pair(singlecell_path, bulk_path, output, block_rows=256, atol=ATOL, rtol=RTOL):
    """Content audit, separated from the public source gate for tiny fixture tests."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if not 1 <= block_rows <= 512:
        raise ValueError("block_rows must be between 1 and 512")
    with h5py.File(singlecell_path, "r") as sc, h5py.File(bulk_path, "r") as bulk:
        if not isinstance(sc["X"], h5py.Dataset) or sc["X"].ndim != 2:
            raise ValueError("This source audit requires the deposited dense two-dimensional X")
        if not isinstance(bulk["X"], h5py.Dataset) or bulk["X"].ndim != 2:
            raise ValueError("Bulk companion must have a dense two-dimensional X")
        write_json(output / "schema.json", {"singlecell": matrix_schema(sc), "bulk": matrix_schema(bulk)})
        sc_ids, bulk_ids = strings(sc["var"], "gene_id"), strings(bulk["var"], "gene_id")
        sc_names, bulk_names = strings(sc["var"], "gene_name"), strings(bulk["var"], "gene_name")
        if (len(sc_ids) != sc["X"].shape[1] or len(bulk_ids) != bulk["X"].shape[1]
                or not np.array_equal(sc_ids, bulk_ids)
                or not np.array_equal(sc_names, bulk_names)
                or len(set(sc_ids)) != len(sc_ids) or any(g is None for g in sc_ids)
                or any(re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?", g) is None for g in sc_ids)):
            raise ValueError("Gene identity, ordered Ensembl IDs, or decoded symbols disagree")
        labels = strings(sc["obs"], "gene_transcript")
        bulk_labels = strings(bulk["obs"], "gene_transcript")
        expected_counts = np.asarray(read_column(bulk["obs"], "num_cells_filtered"), float)
        if len(labels) != sc["X"].shape[0] or len(bulk_labels) != bulk["X"].shape[0]:
            raise ValueError("Expression rows do not match perturbation-label rows")
        if len(expected_counts) != len(bulk_labels):
            raise ValueError("Bulk count metadata dimension mismatch")
        known_labels = [v for v in bulk_labels if v is not None]
        if len(set(known_labels)) != len(known_labels):
            raise ValueError("Duplicated nonmissing bulk perturbation labels")
        lookup = {v: i for i, v in enumerate(bulk_labels) if v is not None}
        codes = np.asarray([lookup.get(v, -1) for v in labels], dtype=np.int64)
        unknown = Counter(v for v, code in zip(labels, codes) if code < 0)
        observed_counts = np.bincount(codes[codes >= 0], minlength=len(bulk_labels))
        sums = np.zeros(bulk["X"].shape, dtype=np.float64)
        if sums.nbytes > 200_000_000:
            raise ValueError("Aggregation buffer exceeds the fixed 200 MB limit")
        invalid = {"nonfinite": 0, "negative": 0, "noninteger": 0}
        total_umi = 0.0
        minimum, maximum = float("inf"), float("-inf")
        for start in range(0, len(labels), block_rows):
            stop = min(start + block_rows, len(labels))
            x = sc["X"][start:stop]
            finite = np.isfinite(x)
            invalid["nonfinite"] += int((~finite).sum())
            invalid["negative"] += int((finite & (x < 0)).sum())
            invalid["noninteger"] += int((finite & (x != np.floor(x))).sum())
            if finite.any():
                minimum = min(minimum, float(x[finite].min()))
                maximum = max(maximum, float(x[finite].max()))
            total_umi += float(np.sum(x, dtype=np.float64))
            chunk_codes = codes[start:stop]
            valid = chunk_codes >= 0
            np.add.at(sums, chunk_codes[valid], x[valid])
            if start == 0 or stop == len(labels) or stop % (block_rows * 128) == 0:
                print(f"Audited raw X: {stop:,}/{len(labels):,} cells", flush=True)
        label_rows, missing_count_rows, mean_mismatches = [], [], []
        count_mismatches = 0
        invalid_bulk_counts = 0
        mismatch_entries = 0
        compared_entries = 0
        compared_labels = 0
        max_abs, max_rel, max_scaled = 0.0, 0.0, 0.0
        bulk_nonfinite = 0
        bulk_negative = 0
        for i, label in enumerate(bulk_labels):
            n, expected = int(observed_counts[i]), float(expected_counts[i])
            if np.isfinite(expected):
                valid_count = expected >= 0 and expected.is_integer()
                invalid_bulk_counts += int(not valid_count)
                count_ok = valid_count and n == int(expected)
                count_mismatches += int(not count_ok)
                count_status = "match" if count_ok else "mismatch"
            else:
                count_status = "missing_author_count_not_imputed" if np.isnan(expected) else "invalid_infinite_author_count"
                invalid_bulk_counts += int(np.isinf(expected))
                missing_count_rows.append({"bulk_row": i, "gene_transcript": label,
                                           "observed_singlecell_count": n, "status": count_status})
            row = {"bulk_row": i, "gene_transcript": label, "singlecell_count": n,
                   "author_num_cells_filtered": expected if np.isfinite(expected) else None,
                   "count_status": count_status, "mean_status": "no_singlecell_cells_not_imputed",
                   "mean_mismatch_genes": None, "max_absolute_error": None,
                   "max_relative_error_nonzero_reference": None, "max_tolerance_scaled_error": None}
            reference = np.asarray(bulk["X"][i], dtype=np.float64)
            bulk_nonfinite += int((~np.isfinite(reference)).sum())
            bulk_negative += int((np.isfinite(reference) & (reference < 0)).sum())
            if n:
                mean = sums[i] / n
                error = np.abs(mean - reference)
                allowed = atol + rtol * np.abs(reference)
                ok = np.isfinite(mean) & np.isfinite(reference) & (error <= allowed)
                finite_error = np.isfinite(error)
                absolute = float(error[finite_error].max()) if finite_error.any() else None
                nz = finite_error & (reference != 0)
                relative = float((error[nz] / np.abs(reference[nz])).max()) if nz.any() else None
                scaled = float((error[finite_error] / allowed[finite_error]).max()) if finite_error.any() else None
                number = int((~ok).sum())
                compared_labels += 1
                compared_entries += len(reference)
                mismatch_entries += number
                max_abs = max(max_abs, absolute or 0.0)
                max_rel = max(max_rel, relative or 0.0)
                max_scaled = max(max_scaled, scaled or 0.0)
                row.update(mean_status="match" if number == 0 else "mismatch", mean_mismatch_genes=number,
                           max_absolute_error=absolute, max_relative_error_nonzero_reference=relative,
                           max_tolerance_scaled_error=scaled)
                if number:
                    worst = np.flatnonzero(~ok)[np.argsort(np.nan_to_num(error[~ok], nan=np.inf))[-5:]]
                    mean_mismatches.extend({"gene_transcript": label, "gene_id": sc_ids[j],
                        "singlecell_mean": float(mean[j]) if np.isfinite(mean[j]) else None,
                        "author_bulk_mean": float(reference[j]) if np.isfinite(reference[j]) else None,
                        "absolute_error": float(error[j]) if np.isfinite(error[j]) else None}
                        for j in worst)
            label_rows.append(row)
        with (output / "label_comparison.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(label_rows[0]))
            writer.writeheader()
            writer.writerows(label_rows)
        write_json(output / "missing_author_counts.json", missing_count_rows)
        write_json(output / "unmapped_singlecell_labels.json", [
            {"gene_transcript": label, "singlecell_cells": n} for label, n in unknown.items()])
        write_json(output / "mean_mismatch_examples.json", mean_mismatches)
        passed = (not any(invalid.values()) and not unknown and count_mismatches == 0
                  and invalid_bulk_counts == 0 and mismatch_entries == 0 and bulk_nonfinite == 0
                  and bulk_negative == 0
                  and not any(label is None for label in bulk_labels))
        return {"status": "contents_concordant_with_missing_author_count_metadata" if passed and missing_count_rows
                         else "contents_concordant" if passed else "contents_discrepancy",
                "passed": passed, "singlecell_shape": list(sc["X"].shape),
                "singlecell_dtype": str(sc["X"].dtype), "bulk_shape": list(bulk["X"].shape),
                "gene_ids_unique_and_identical_order": True, "decoded_gene_symbols_identical_order": True,
                "raw_X_invalid_entries": invalid, "raw_X_min": minimum if np.isfinite(minimum) else None,
                "raw_X_max": maximum if np.isfinite(maximum) else None,
                "raw_X_total_UMI": int(total_umi) if np.isfinite(total_umi) and total_umi.is_integer() else None,
                "bulk_nonfinite_X_entries": bulk_nonfinite, "bulk_negative_X_entries": bulk_negative,
                "aggregate_buffer_bytes": sums.nbytes,
                "block_rows": block_rows, "cells_with_known_bulk_labels": int(observed_counts.sum()),
                "cells_with_missing_or_unknown_labels": sum(unknown.values()),
                "finite_author_count_rows": int(np.isfinite(expected_counts).sum()),
                "missing_author_count_rows": len(missing_count_rows),
                "invalid_author_count_rows": invalid_bulk_counts, "count_mismatch_rows": count_mismatches,
                "compared_mean_labels": compared_labels, "compared_mean_entries": compared_entries,
                "mean_mismatch_entries": mismatch_entries, "max_absolute_mean_error": max_abs,
                "max_relative_mean_error_nonzero_reference": max_rel, "max_tolerance_scaled_error": max_scaled,
                "atol": atol, "rtol": rtol,
                "tolerance_rule": "abs(singlecell_mean-author_bulk_mean) <= atol + rtol*abs(author_bulk_mean)",
                "missing_data_policy": "Missing author counts are not imputed. Labels with zero observed single cells have no mean comparison. Missing/unknown single-cell labels fail the gate.",
                "scope": "Author-deposit consistency, not FASTQ reprocessing, biological validation, or proof of raw-read provenance."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-rows", type=int, default=256)
    args = parser.parse_args()
    if (OUTPUT / "manifest.json").exists():
        raise SystemExit("Completed content audit exists; refusing to overwrite it")
    start = datetime.now(timezone.utc).isoformat()
    source = source_gate(SINGLECELL, RECEIPTS / "Replogle_RPE1_raw_singlecell.json", EXPECTED_BYTES, EXPECTED_MD5, 35775606)
    companion = source_gate(BULK, RECEIPTS / "Replogle_RPE1_raw_bulk.json", BULK_BYTES, BULK_MD5, 35775581)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    script_hash = file_hashes(Path(__file__))["sha256"]
    try:
        result = audit_h5ad_pair(SINGLECELL, BULK, OUTPUT, block_rows=args.block_rows)
        for item in (source, companion):
            stat = Path(item["path"]).stat()
            if (stat.st_size, stat.st_mtime_ns) != (item["bytes"], item["mtime_ns"]):
                raise ValueError("Verified source changed during content audit")
        if file_hashes(Path(__file__))["sha256"] != script_hash:
            raise ValueError("Audit script changed during execution")
        write_json(OUTPUT / "summary.json", result)
        write_json(OUTPUT / "manifest.json", {"status": result["status"], "started_utc": start,
            "completed_utc": datetime.now(timezone.utc).isoformat(), "sources": [source, companion],
            "runtime": {"python": sys.version, "numpy": np.__version__, "h5py": h5py.__version__},
            "script_sha256": script_hash, "tolerances_fixed_before_opening_singlecell": {"atol": ATOL, "rtol": RTOL},
            "output_sha256": {p.name: file_hashes(p)["sha256"] for p in OUTPUT.iterdir() if p.is_file()}})
        print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False), flush=True)
        if not result["passed"]:
            raise SystemExit(2)
    except Exception as exc:
        write_json(OUTPUT / "failure.json", {"status": "audit_failed_without_source_modification",
            "error": repr(exc), "started_utc": start, "script_sha256": script_hash})
        raise


if __name__ == "__main__":
    main()
