"""Tiny artificial fixtures test the RPE1 provenance reader, not research outcomes."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from verify_rpe1_contents import audit_h5ad_pair, read_column, source_gate


def put_strings(group, name, values):
    return group.create_dataset(name, data=np.asarray(values, dtype=object), dtype=h5py.string_dtype())


class RPE1MetadataChecks(unittest.TestCase):
    def test_legacy_reference_decodes_missing_without_negative_indexing(self):
        with tempfile.TemporaryDirectory() as td, h5py.File(Path(td)/"legacy.h5", "w") as f:
            g = f.create_group("obs")
            categories = put_strings(g, "category_values", ["alpha", "beta"])
            d = g.create_dataset("labels", data=np.array([1, -1, 0], dtype=np.int16))
            d.attrs["categories"] = categories.ref
            self.assertEqual(read_column(g, "labels").tolist(), ["beta", None, "alpha"])

    def test_modern_categorical_and_plain_numeric(self):
        with tempfile.TemporaryDirectory() as td, h5py.File(Path(td)/"modern.h5", "w") as f:
            g = f.create_group("obs"); cat = g.create_group("labels")
            cat.create_dataset("codes", data=np.array([0, 1, -1]))
            put_strings(cat, "categories", ["a", "b"])
            g.create_dataset("numeric", data=[1.5, 2.5, np.nan])
            self.assertEqual(read_column(g, "labels").tolist(), ["a", "b", None])
            np.testing.assert_equal(read_column(g, "numeric"), [1.5, 2.5, np.nan])

    def test_bad_categorical_codes_raise(self):
        with tempfile.TemporaryDirectory() as td, h5py.File(Path(td)/"bad.h5", "w") as f:
            cat = f.create_group("labels")
            cat.create_dataset("codes", data=np.array([2]))
            put_strings(cat, "categories", ["a", "b"])
            with self.assertRaises(ValueError): read_column(f, "labels")

    def test_source_gate_rejects_incomplete_file_before_hdf5(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/"not_an_hdf5.h5ad"; p.write_bytes(b"unfinished")
            with self.assertRaisesRegex(ValueError, "incomplete or wrong size"):
                source_gate(p, Path(td)/"absent.json", 8_700_873_216, "fixed", 1)
            with self.assertRaisesRegex(ValueError, "receipt is absent"):
                source_gate(p, Path(td)/"absent.json", p.stat().st_size, "fixed", 1)

    def make_pair(self, directory, bad_count=False, bad_mean=False, bad_raw=False):
        directory = Path(directory)
        sc, bulk = directory/"sc.h5ad", directory/"bulk.h5ad"
        x = np.array([[1, 3], [3, 1], [0, 2]], dtype=np.float32)
        if bad_raw: x[0, 0] = 1.5
        for path, matrix, labels in [(sc, x, ["p", "p", "nt"]),
                                    (bulk, np.array([[2, 2], [0, 2], [7, 8]], dtype=np.float32), ["p", "nt", "absent"] )]:
            with h5py.File(path, "w") as f:
                f.create_dataset("X", data=matrix)
                obs = f.create_group("obs"); put_strings(obs, "gene_transcript", labels)
                var = f.create_group("var")
                put_strings(var, "gene_id", ["ENSG00000000001", "ENSG00000000002"])
                cats = put_strings(var, "names", ["A", "B"])
                codes = var.create_dataset("gene_name", data=np.array([0, 1], dtype=np.int16))
                codes.attrs["categories"] = cats.ref
                if path == bulk:
                    obs.create_dataset("num_cells_filtered", data=[3 if bad_count else 2, 1, np.nan])
                    if bad_mean: f["X"][0, 0] = 2.1
        return sc, bulk

    def test_streamed_counts_means_and_missing_count_not_imputed(self):
        with tempfile.TemporaryDirectory() as td:
            sc, bulk = self.make_pair(td)
            result = audit_h5ad_pair(sc, bulk, Path(td)/"out", block_rows=1)
            self.assertTrue(result["passed"])
            self.assertEqual(result["raw_X_total_UMI"], 10)
            self.assertEqual(result["compared_mean_entries"], 4)
            self.assertEqual(result["missing_author_count_rows"], 1)
            missing = json.loads((Path(td)/"out/missing_author_counts.json").read_text())
            self.assertEqual(missing[0]["status"], "missing_author_count_not_imputed")
            self.assertEqual(missing[0]["observed_singlecell_count"], 0)

    def test_count_or_mean_or_integer_discrepancy_fails_without_correction(self):
        for flag, key in [("bad_count", "count_mismatch_rows"), ("bad_mean", "mean_mismatch_entries"),
                          ("bad_raw", "raw_X_invalid_entries")]:
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as td:
                sc, bulk = self.make_pair(td, **{flag: True})
                result = audit_h5ad_pair(sc, bulk, Path(td)/"out", block_rows=2)
                self.assertFalse(result["passed"])
                value = result[key]["noninteger"] if flag == "bad_raw" else result[key]
                self.assertGreater(value, 0)

    def test_missing_author_count_still_compares_mean_when_cells_exist(self):
        with tempfile.TemporaryDirectory() as td:
            sc, bulk = self.make_pair(td)
            with h5py.File(bulk, "r+") as f:
                f["obs/num_cells_filtered"][1] = np.nan
            result = audit_h5ad_pair(sc, bulk, Path(td)/"out", block_rows=2)
            self.assertTrue(result["passed"])
            self.assertEqual(result["compared_mean_entries"], 4)
            missing = json.loads((Path(td)/"out/missing_author_counts.json").read_text())
            self.assertEqual(missing[0]["gene_transcript"], "nt")
            self.assertEqual(missing[0]["observed_singlecell_count"], 1)


if __name__ == "__main__":
    unittest.main()
