"""Mathematical checks against enumeration, rather than old implementation."""
import itertools
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from carf_core import (topk_weights, rank_metrics, conditional_moments,
                       permutation_pvalue, residualize, bootstrap_mean, bh_adjust)


class TestCore(unittest.TestCase):
    def test_ties_equal_uniform_boundary_randomization(self):
        score = [5, 3, 3, 3, 1]
        y = np.array([0, 1, 0, 1, 0])
        hits = [y[[0, i]].sum() for i in [1, 2, 3]]
        self.assertAlmostEqual(topk_weights(score, 2) @ y, np.mean(hits))

    def test_constant_scores_have_no_arbitrary_gene_preference(self):
        result = rank_metrics(np.zeros(20), [1] * 3 + [0] * 17, 5)
        self.assertEqual(result["scope_enrichment"], 1.0)
        self.assertFalse(result["rankable"])
        self.assertAlmostEqual(result["precision"], 3 / 20)

    def test_scope_null_expectation_and_coverage_are_separate(self):
        y = [1, 0, 0, 0]
        scores = [4, 3, 2, 1]
        metrics = [rank_metrics(scores, p, 2, 3) for p in set(itertools.permutations(y))]
        self.assertAlmostEqual(np.mean([m["scope_enrichment"] for m in metrics]), 1)
        self.assertTrue(all(m["reference_coverage"] == 1 / 3 for m in metrics))

    def test_conditional_moments_match_exhaustive_enumeration(self):
        w = np.array([1, .5, 0, 1, 0, 0])
        y = np.array([1, 0, 0, 1, 1, 0])
        strata = [0, 0, 0, 1, 1, 1]
        vals = [w @ np.array(a + b)
                for a in set(itertools.permutations(y[:3]))
                for b in set(itertools.permutations(y[3:]))]
        m = conditional_moments(w, y, strata)
        self.assertAlmostEqual(m["reference_mean_hits"], np.mean(vals))
        self.assertAlmostEqual(m["reference_variance_hits"], np.var(vals))

    def test_conditional_nonidentifiability(self):
        m = conditional_moments([1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 1, 1])
        self.assertFalse(m["within_stratum_information"])
        p = permutation_pvalue([1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 1, 1], np.random.default_rng(2))
        self.assertEqual(p, 1)

    def test_invalid_scores_fail_instead_of_becoming_zero(self):
        for x in ([1, np.nan], [1, np.inf], []):
            with self.assertRaises(ValueError):
                topk_weights(x, 1)
        with self.assertRaises(ValueError):
            topk_weights([1, 2], 3)

    def test_joint_residuals_are_orthogonal(self):
        rng = np.random.default_rng(2)
        x = rng.normal(size=(100, 3))
        y = x @ [1, 2, 3] + rng.normal(size=100)
        r = residualize(y, x)
        np.testing.assert_allclose(x.T @ r, 0, atol=1e-10)
        np.testing.assert_allclose(residualize(np.ones(100), x), 0)

    def test_bootstrap_retains_multiplicity(self):
        class FixedRNG:
            def integers(self, low, high, size):
                return np.array([[0, 0, 1], [2, 2, 2]])
        ci = bootstrap_mean([1, 2, 9], FixedRNG(), 2)
        np.testing.assert_allclose(ci, np.quantile([4 / 3, 9], [.025, .975]))

    def test_bh(self):
        np.testing.assert_allclose(bh_adjust([.04, .001, .03]), [.04, .003, .04])


if __name__ == "__main__":
    unittest.main()
