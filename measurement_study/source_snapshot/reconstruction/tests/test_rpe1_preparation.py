"""Design invariants for metadata-first RPE1 preparation; no real X read."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_rpe1_measurement import (largest_remainder, proportional_sample, matched_controls,
                                      target_folds, normalized_counts, rng_for)


class RPE1PreparationChecks(unittest.TestCase):
    def test_integer_quotas_and_tie_order(self):
        np.testing.assert_array_equal(largest_remainder([3, 3, 3], 2), [1, 1, 0])
        np.testing.assert_array_equal(largest_remainder([1, 4, 5], 7), [1, 3, 3])
        for n in range(11):
            q = largest_remainder([1, 4, 5], n)
            self.assertEqual(q.sum(), n)
            self.assertTrue(np.all(q <= [1, 4, 5]))
        with self.assertRaises(ValueError): largest_remainder([1, 2], 4)

    def test_sampling_reproducible_and_order_invariant(self):
        cells = np.arange(10)
        strata = np.array([2]*4+[1]*6)
        a = proportional_sample(cells, strata, 5, "test", 2)
        b = proportional_sample(cells[::-1], strata[::-1], 5, "test", 2)
        np.testing.assert_array_equal(a, b)
        self.assertEqual(np.sum(strata[a] == 1), 3)
        self.assertEqual(np.sum(strata[a] == 2), 2)

    def test_exact_GEM_control_matching_and_capacity_failure(self):
        gems = np.array([1, 1, 2, 1, 1, 2, 2, 3])
        selected = np.array([0, 1, 2]); pool = np.array([3, 4, 5, 6, 7])
        control = matched_controls(selected, pool, gems, "matched", 3)
        np.testing.assert_array_equal(np.sort(gems[control]), np.sort(gems[selected]))
        self.assertEqual(len(set(control)), 3)
        with self.assertRaises(ValueError): matched_controls(selected, np.array([3, 5]), gems, "matched", 3)

    def test_folds_are_ID_defined_and_balanced(self):
        ids = [f"ENSG{i:011d}" for i in range(505)]
        a, b = target_folds(ids), target_folds(ids[::-1])
        self.assertEqual(a, b)
        self.assertEqual(sorted([sum(v == k for v in a.values()) for k in range(3)]), [168, 168, 169])

    def test_normalization_full_observed_gene_sum_and_zero_failure(self):
        x = np.array([[1, 9, 90], [2, 18, 180]], dtype=np.uint32)
        norm, totals = normalized_counts(x)
        np.testing.assert_allclose(norm, np.tile(np.log1p([100, 900, 9000]), (2, 1)))
        np.testing.assert_array_equal(totals, [100, 200])
        with self.assertRaises(ValueError): normalized_counts([[0, 0], [1, 2]])
        with self.assertRaises(ValueError): normalized_counts([[0.5, 1]])

    def test_rng_streams_differ_and_are_reproducible(self):
        np.testing.assert_array_equal(rng_for("x", 3).integers(100, size=10), rng_for("x", 3).integers(100, size=10))
        self.assertNotEqual(rng_for("x", 3).integers(2**63), rng_for("y", 3).integers(2**63))


if __name__ == "__main__":
    unittest.main()
