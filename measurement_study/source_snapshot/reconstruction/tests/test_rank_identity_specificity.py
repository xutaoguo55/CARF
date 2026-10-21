import unittest

import numpy as np

from rank_identity_specificity import evaluate_cross_task


class IdentitySpecificityTests(unittest.TestCase):
    def test_shared_ranking_has_zero_identity_gain(self):
        shared = np.array([5.0, 4.0, 3.0, 2.0, 1.0])
        scores = np.tile(shared, (3, 1))
        B = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0]], float)
        C = B.copy()
        _, rows = evaluate_cross_task(scores, B, C, k=1)
        np.testing.assert_allclose(rows.same_minus_other_mean, 0.0, atol=0, rtol=0)

    def test_task_aligned_rank_beats_mismatched_ranks(self):
        effects = np.array([
            [5, 1, 0, 0, 0],
            [0, 5, 1, 0, 0],
            [0, 0, 5, 1, 0],
        ], dtype=float)
        # Each task ranking is its own known absolute response profile.
        scores = effects.copy()
        _, rows = evaluate_cross_task(scores, effects, effects, k=1)
        self.assertTrue((rows.same_minus_other_mean > 0).all())
        self.assertTrue(rows.same_exceeds_other_mean.all())


if __name__ == "__main__":
    unittest.main()
