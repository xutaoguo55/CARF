import io
from pathlib import Path
import pickle
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from export_rpe1_safe_models import safe_copy_score, read_checked_string_metadata, StringMetadataUnpickler


class UnexpectedPayload:
    def __reduce__(self):
        return (eval, ("1 + 1",))


class RPE1SafeExportChecks(unittest.TestCase):
    def test_actual_generated_object_metadata_becomes_pickle_free(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/"source.npz"; destination = Path(directory)/"safe.npz"
            targets = np.asarray(["ENSG00000000001", "ENSG00000000002"], dtype=object)
            genes = np.asarray(["ENSG00000000003", "ENSG00000000004", "ENSG00000000005"], dtype=object)
            raw = np.arange(6, dtype=float).reshape(2, 3); signed = raw - 2; folds = np.array([0, 1])
            np.savez_compressed(source, targets=targets, genes=genes, raw_scores=raw,
                                signed_predictions=signed, folds=folds, target_zero_variance=np.array([False, True]))
            checks, numeric = safe_copy_score(source, destination, targets, genes, folds)
            self.assertTrue(all(c["restricted_unpickler_used"] for c in checks))
            self.assertIn("signed_predictions", numeric)
            with np.load(destination, allow_pickle=False) as z:
                for key in z.files: self.assertFalse(z[key].dtype.hasobject)
                self.assertEqual(z["targets"].dtype.kind, "U")
                np.testing.assert_array_equal(z["raw_scores"], raw)
                np.testing.assert_array_equal(z["signed_predictions"], signed)

    def test_restricted_unpickler_rejects_unexpected_global(self):
        payload = pickle.dumps(UnexpectedPayload())
        with self.assertRaises(pickle.UnpicklingError): StringMetadataUnpickler(io.BytesIO(payload)).load()

    def test_mismatched_frozen_identity_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"source.npz"
            np.savez(path, targets=np.array(["actual"], dtype=object))
            with self.assertRaises(ValueError): read_checked_string_metadata(path, "targets", ["different"])


if __name__ == "__main__": unittest.main()
