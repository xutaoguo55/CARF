"""Derive pickle-free RPE1 score exports without changing a numeric prediction.

The original score writer used pandas object arrays for string metadata. We
verify their trusted, hash-bound content with a restricted NumPy-only unpickler,
then re-export Unicode IDs obtained from the frozen scope CSVs. Numeric arrays
are read with allow_pickle=False and checked exactly after writing.
"""
from __future__ import annotations
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import pickle
import zipfile

import numpy as np
import numpy._core.multiarray as multiarray
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE / "results/rpe1_measurement_validation_v1"
SOURCE = BASE / "models"
PREP = BASE / "preparation"
OUT = BASE / "models_safe_export"


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


class StringMetadataUnpickler(pickle.Unpickler):
    """Permit only NumPy's minimal object-array reconstruction constructors."""
    def find_class(self, module, name):
        allowed = {("numpy", "ndarray"): np.ndarray, ("numpy", "dtype"): np.dtype,
                   ("numpy._core.multiarray", "_reconstruct"): multiarray._reconstruct,
                   ("numpy.core.multiarray", "_reconstruct"): multiarray._reconstruct}
        if (module, name) not in allowed:
            raise pickle.UnpicklingError(f"Forbidden object metadata constructor: {module}.{name}")
        return allowed[(module, name)]


def read_checked_string_metadata(path, key, expected):
    """Never call np.load(..., allow_pickle=True). Only exact frozen string IDs pass."""
    expected = np.asarray(expected, dtype=str)
    with zipfile.ZipFile(path) as archive:
        info = archive.getinfo(key + ".npy")
        if info.file_size > 2_000_000:
            raise ValueError("Unexpectedly large string metadata member")
        stream = io.BytesIO(archive.read(info))
    version = np.lib.format.read_magic(stream)
    if version == (1, 0):
        shape, _, dtype = np.lib.format.read_array_header_1_0(stream)
    elif version == (2, 0):
        shape, _, dtype = np.lib.format.read_array_header_2_0(stream)
    else:
        raise ValueError(f"Unsupported metadata NPY format version {version}")
    if shape != expected.shape or len(shape) != 1:
        raise ValueError("String metadata shape differs from frozen scope")
    if dtype.hasobject:
        observed = StringMetadataUnpickler(stream).load()
        if type(observed) is not np.ndarray or observed.dtype != np.dtype(object) or observed.shape != shape:
            raise ValueError("Unexpected reconstructed metadata array")
        if not all(type(v) is str for v in observed):
            raise ValueError("Object metadata contains non-string elements")
        if stream.read():
            raise ValueError("Trailing object metadata payload")
    else:
        with np.load(path, allow_pickle=False) as z:
            observed = z[key]
        if observed.dtype.kind not in "US":
            raise ValueError("Metadata is not a string array")
        observed = observed.astype(str)
    if not np.array_equal(observed, expected):
        raise ValueError(f"{key} metadata differs from frozen scope")
    return dict(key=key, original_dtype=str(dtype), n_elements=len(expected),
                exact_frozen_identity_match=True, restricted_unpickler_used=bool(dtype.hasobject))


def safe_copy_score(source, destination, targets, genes, folds):
    target_ids = np.asarray(targets, dtype=str)
    gene_ids = np.asarray(genes, dtype=str)
    folds = np.asarray(folds, dtype=np.int64)
    checks = [read_checked_string_metadata(source, "targets", target_ids),
              read_checked_string_metadata(source, "genes", gene_ids)]
    arrays = {"targets": target_ids, "genes": gene_ids}
    numeric = {}
    with np.load(source, allow_pickle=False) as z:
        for key in z.files:
            if key in ["targets", "genes"]:
                continue
            value = z[key]
            if value.dtype.hasobject or value.dtype.kind not in "biufc":
                raise ValueError(f"Unexpected nonnumeric score member {key}")
            numeric[key] = value
    if not np.array_equal(numeric["folds"], folds):
        raise ValueError("Numeric fold IDs differ from frozen scope")
    if numeric["raw_scores"].shape != (len(target_ids), len(gene_ids)):
        raise ValueError("Raw score dimensions differ from frozen scope")
    arrays.update(numeric)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, **arrays)
    with np.load(destination, allow_pickle=False) as out:
        if set(out.files) != set(arrays):
            raise ValueError("Export member keys changed")
        for key in out.files:
            value = out[key]  # Integration check: every member loads without pickle.
            if value.dtype.hasobject:
                raise ValueError("Object dtype persisted in safe export")
            expected = arrays[key]
            same = np.array_equal(value, expected, equal_nan=True) if value.dtype.kind in "fc" else np.array_equal(value, expected)
            if value.dtype != expected.dtype or not same:
                raise ValueError(f"Export changed array contents or numeric dtype: {key}")
    return checks, sorted(numeric)


def verify_manifest(folder):
    path = folder / "manifest.json"
    record = json.loads(path.read_text())
    for name, expected in record["output_sha256"].items():
        if file_hash(folder/name) != expected:
            raise ValueError(f"Source manifest hash mismatch: {folder/name}")
    return record, file_hash(path)


def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise SystemExit("Safe export directory is nonempty; refusing overwrite")
    start = datetime.now(timezone.utc).isoformat()
    code_hash = file_hash(__file__)
    source_manifest, source_hash = verify_manifest(SOURCE)
    preparation_manifest, preparation_hash = verify_manifest(PREP)
    if source_manifest["preparation_manifest_sha256"] != preparation_hash:
        raise ValueError("Original models do not link to exact preparation")
    if source_manifest["protocol_sha256"] != preparation_manifest["protocol_sha256"]:
        raise ValueError("Original models and preparation protocols differ")
    targets = pd.read_csv(PREP/"evaluation_targets.csv")
    genes = pd.read_csv(PREP/"evaluation_genes.csv")
    target_ids = targets.target_id.to_numpy(dtype=str)
    gene_ids = genes.gene_id.to_numpy(dtype=str)
    folds = targets.fold.to_numpy(dtype=np.int64)
    scope_hashes = {name: file_hash(PREP/name) for name in ["evaluation_targets.csv", "evaluation_genes.csv"]}
    source_index = pd.read_csv(SOURCE/"model_index.csv")
    if len(source_index) != 35:
        raise ValueError("Expected 35 original score exports")
    OUT.mkdir(parents=True)
    exported = []; audit = []
    for record in source_index.to_dict("records"):
        original = SOURCE/record["path"]
        if file_hash(original) != record["sha256"]:
            raise ValueError("Original score index hash mismatch")
        destination = OUT/record["path"]
        metadata_checks, numeric_keys = safe_copy_score(original, destination, target_ids, gene_ids, folds)
        audit.append(dict(path=record["path"], source_sha256=record["sha256"],
                          safe_export_sha256=file_hash(destination), metadata_checks=metadata_checks,
                          numeric_keys_preserved=numeric_keys, all_numeric_arrays_exactly_equal=True,
                          every_export_member_loads_with_allow_pickle_false=True))
        exported.append(dict(record, sha256=file_hash(destination), source_sha256=record["sha256"]))
    pd.DataFrame(exported).to_csv(OUT/"model_index.csv", index=False)
    (OUT/"export_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    if file_hash(__file__) != code_hash or file_hash(SOURCE/"manifest.json") != source_hash or file_hash(PREP/"manifest.json") != preparation_hash:
        raise ValueError("Source or export implementation changed")
    for name, expected in scope_hashes.items():
        if file_hash(PREP/name) != expected:
            raise ValueError("Frozen scope changed during export")
    outputs = {str(p.relative_to(OUT)): file_hash(p) for p in sorted(OUT.rglob("*")) if p.is_file()}
    manifest = dict(status="completed_safe_metadata_export_without_numeric_changes", started_utc=start,
        finished_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=source_manifest["protocol_sha256"],
        preparation_manifest_sha256=preparation_hash, source_model_manifest_sha256=source_hash,
        export_code_sha256=code_hash, frozen_scope_sha256=scope_hashes, output_sha256=outputs,
        n_model_exports=35, numerical_changes=False,
        repair="Original pandas object string IDs were independently checked against frozen scope with a restricted NumPy-only unpickler and exported as Unicode. Every numeric member and dtype is unchanged; no model was refit.",
        source_immutability="Original model files, manifest, runner and protocol remain unchanged.")
    (OUT/"manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Completed 35 pickle-free metadata exports; every numeric array unchanged", flush=True)


if __name__ == "__main__":
    main()
