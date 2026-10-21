import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_rpe1_measurement_models import (rng_stream, thin_full_counts,
    normalize_full_counts, frozen_target_folds, fit_selected_targets)


class RPE1ModelChecks(unittest.TestCase):
    def test_normalize_before_gene_subset(self):
        raw = np.array([[2, 3, 5], [4, 1, 15]])
        normalized, total = normalize_full_counts(raw)
        np.testing.assert_array_equal(total, [10, 20])
        np.testing.assert_allclose(normalized[:, :2], np.log1p(raw[:, :2] / total[:, None] * 10000))
        self.assertFalse(np.allclose(normalized[:, :2], np.log1p(raw[:, :2] / raw[:, :2].sum(1)[:, None] * 10000)))

    def test_thinning_streams_reproducible_independent_and_bounded(self):
        raw = np.arange(160).reshape(10, 16)
        a = thin_full_counts(raw, .25, rng_stream("p", "A|target=t1", 5))
        b = thin_full_counts(raw, .25, rng_stream("p", "A|target=t1", 5))
        c = thin_full_counts(raw, .25, rng_stream("p", "A|target=t2", 5))
        np.testing.assert_array_equal(a, b)
        self.assertFalse(np.array_equal(a, c)); self.assertTrue(np.all((a >= 0) & (a <= raw)))
        np.testing.assert_array_equal(thin_full_counts(raw, 1, rng_stream("p", "A", 5)), raw)

    def test_zero_library_fails_without_replacement(self):
        with self.assertRaises(ValueError): normalize_full_counts(np.array([[0, 0], [1, 2]]))

    def test_folds_stable_under_reordering(self):
        names = [f"ENS{i}" for i in range(20)]
        self.assertEqual(dict(zip(names, frozen_target_folds(names, "p"))),
                         dict(zip(names[::-1], frozen_target_folds(names[::-1], "p"))))

    def test_fold_specific_prediction_unaffected_by_heldout_A(self):
        rng = np.random.default_rng(2); post = rng.normal(size=(18, 12))
        folds = np.arange(12) % 3; targets = np.arange(12); evaluated = np.array([0, 1, 2, 3, 4, 5])
        a, _ = fit_selected_targets(post, folds, targets, evaluated, np.zeros(18), 2, .1, [11], "p")
        changed = post.copy(); changed[:, folds == 0] += rng.normal(100, 10, (18, 4))
        b, _ = fit_selected_targets(changed, folds, targets, evaluated, np.zeros(18), 2, .1, [11], "p")
        heldout = folds[evaluated] == 0
        for model in a:
            np.testing.assert_array_equal(a[model][heldout], b[model][heldout])


if __name__ == "__main__": unittest.main()
