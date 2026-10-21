"""Independent synthetic review; never open real RPE1 outcomes or model files."""
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import evaluate_rpe1_measurement as ev


class RPE1EvaluationReview(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(713)
        self.folds = np.array([0, 0, 1, 1, 2, 2])
        self.B = self.rng.normal(size=(6, 7))
        self.C = self.rng.normal(size=(6, 7))
        self.raw = self.rng.uniform(.05, 2, size=(6, 7))
        self.cov = np.linspace(0, 3, 7)
        self.null = self.rng.normal(size=(4, 7))

    def test_fold_invariant_ranking_has_zero_identity_without_crossfold_mixing(self):
        raw = self.raw.copy()
        for fold in range(3):
            raw[self.folds == fold] = raw[np.flatnonzero(self.folds == fold)[0]]
        for adjustment in ['raw', 'EDR']:
            frame, _, _, utility, _ = ev.basic_evaluation(
                raw, self.cov, self.B, self.C, self.null, self.folds, adjustment, 3)
            np.testing.assert_allclose(frame.identity_gain, 0, atol=1e-14)
            self.assertTrue(np.isnan(utility[self.folds[:, None] != self.folds[None, :]]).all())

    def test_raw_rank_and_response_are_invariant_to_covariate_change(self):
        a = ev.basic_evaluation(self.raw, self.cov, self.B, self.C, self.null, self.folds, 'raw', 3)
        b = ev.basic_evaluation(self.raw, self.cov[::-1] ** 2, self.B, self.C, self.null, self.folds, 'raw', 3)
        for col in ['total_BC', 'other_target_BC', 'identity_gain', 'null_response_mean']:
            np.testing.assert_array_equal(a[0][col], b[0][col])
        for x, y in zip(a[1:], b[1:]):
            np.testing.assert_array_equal(x, y)

    def test_target_permutation_preserves_corresponding_results(self):
        p = np.array([4, 1, 3, 0, 5, 2])
        a = ev.basic_evaluation(self.raw, self.cov, self.B, self.C, self.null, self.folds, 'EDR', 3)
        b = ev.basic_evaluation(self.raw[p], self.cov, self.B[p], self.C[p], self.null, self.folds[p], 'EDR', 3)
        np.testing.assert_allclose(b[0].identity_gain, a[0].identity_gain.to_numpy()[p])
        np.testing.assert_allclose(b[2], a[2][p])
        np.testing.assert_allclose(b[3], a[3][np.ix_(p, p)], equal_nan=True)

    def test_signed_prediction_errors_do_not_use_absolute_or_edr_scores(self):
        signed = self.rng.normal(size=self.raw.shape)
        expected = np.square(signed - self.B).mean(axis=1)
        for adjustment in ['raw', 'EDR']:
            frame, scores, _, _, _ = ev.basic_evaluation(
                np.abs(signed), self.cov, self.B, self.C, self.null,
                self.folds, adjustment, 3, signed)
            np.testing.assert_allclose(frame.unadjusted_prediction_MSE_B, expected)
            if adjustment == 'EDR':
                self.assertTrue((scores < 0).any())
        self.assertFalse(np.allclose(expected, np.square(np.abs(signed) - self.B).mean(axis=1)))

    def test_scale_uses_same_anchor_weights_and_strict_observed_separation(self):
        A = np.ones((6, 7))
        truth = np.ones_like(A)
        null = np.vstack([np.zeros(7), np.full(7, 1.5)])
        scale = ev.scale_reference(A, truth, truth, null, self.folds, 3)
        self.assertAlmostEqual(scale['empirical_A_total'], 1)
        self.assertAlmostEqual(scale['null_partition_mean'], .75)
        self.assertAlmostEqual(scale['total_scale'], .25)
        self.assertFalse(scale['observed_separation'])
        self.assertAlmostEqual(scale['empirical_A_identity'], 0)

    def test_contrast_preserves_direction_and_absolute_changes(self):
        left = pd.DataFrame(dict(total_BC=[2., 3.], identity_gain=[.2, .3],
                                 all_scores_floored=[False, False], score_exact_constant=[False, False]))
        right = left.copy()
        right['total_BC'] = [3., 2.]
        right['identity_gain'] = [.1, .6]
        wl = np.array([[1., 1., 0.], [1., 0., 1.]])
        wr = np.array([[1., 0., 1.], [1., 0., 1.]])
        scale = dict(total_scale=2., observed_separation=True, empirical_A_identity=.5)
        result = ev.covariate_contrast(left, right, wl, wr, scale, 2)
        np.testing.assert_allclose(result.top50_agreement, [.5, 1])
        np.testing.assert_allclose(result.total_BC_change, [1, -1])
        self.assertAlmostEqual(result.total_BC_change.mean(), 0)
        self.assertAlmostEqual(result.total_BC_absolute_change.mean(), 1)
        self.assertAlmostEqual(result.total_absolute_change_over_A_scale.mean(), .5)
        np.testing.assert_allclose(result.identity_gain_change, [-.1, .3])

    def test_all_floor_scores_are_flagged_despite_perfect_weight_agreement(self):
        frame, _, weights, _, _ = ev.basic_evaluation(
            np.zeros_like(self.raw), self.cov, self.B, self.C, self.null, self.folds, 'EDR', 3)
        scale = dict(total_scale=1., observed_separation=True, empirical_A_identity=1.)
        result = ev.covariate_contrast(frame, frame, weights, weights, scale, 3)
        np.testing.assert_allclose(result.top50_agreement, 1)
        self.assertTrue(result.uninformative_at_either_end.all())
        self.assertTrue(frame.all_scores_floored.all())

    def test_unscalable_reference_never_divides_by_its_gap(self):
        frame, _, weights, _, _ = ev.basic_evaluation(
            self.raw, self.cov, self.B, self.C, self.null, self.folds, 'raw', 3)
        scale = dict(total_scale=0., observed_separation=False, empirical_A_identity=0.)
        result = ev.covariate_contrast(frame, frame, weights, weights, scale, 3)
        self.assertTrue(result.total_absolute_change_over_A_scale.isna().all())
        self.assertTrue(result.identity_absolute_change_over_A_identity.isna().all())

    def _manifests(self, model_prep_hash='prep-hash', original_model_hash='synthetic-original-model-hash'):
        cfg_hash = ev.sha256(ev.CONFIG)
        return [({'protocol_sha256': cfg_hash,
                  'status': 'completed_frozen_preparation_no_model_or_thinning'}, 'prep-hash'),
                ({'protocol_sha256': cfg_hash,
                  'status': 'completed_safe_metadata_export_without_numeric_changes',
                  'preparation_manifest_sha256': model_prep_hash,
                  'source_model_manifest_sha256': original_model_hash}, 'model-hash')]

    def test_wrong_preparation_link_rejected_before_loading_any_outcomes(self):
        with tempfile.TemporaryDirectory() as d, \
             patch.object(ev, 'OUT', Path(d)), \
             patch.object(ev, 'verify_manifest', side_effect=self._manifests('other-prep')), \
             patch.object(ev.np, 'load') as loader:
            with self.assertRaisesRegex(ValueError, 'exact preparation'):
                ev.main()
            loader.assert_not_called()

    def test_dc64_with_duplicate_extra_cell_rejected_before_thinning(self):
        z = dict(targets=np.array([f't{i}' for i in range(6)]),
                 genes=np.array([f'g{i}' for i in range(7)]), folds=self.folds,
                 A=self.raw, B=self.B, C=self.C,
                 null_B=np.zeros((30, 7)), null_C=np.zeros((30, 7)))
        candidates = pd.DataFrame(dict(gene_id=z['genes'], source_gene_index=np.arange(7)))
        cells = pd.DataFrame(dict(cell_index=np.arange(400)))
        fake_counts = SimpleNamespace(shape=(400, 8749))
        with tempfile.TemporaryDirectory() as d:
            base = Path(d)
            (base / 'models').mkdir()
            original_manifest = base / 'models/manifest.json'
            original_manifest.write_text('{"synthetic": true}')
            with patch.object(ev, 'OUT', base / 'evaluation'), \
                 patch.object(ev, 'BASE', base), \
                 patch.object(ev, 'verify_manifest', side_effect=self._manifests(
                     original_model_hash=ev.sha256(original_manifest))), \
                 patch.object(ev.np, 'load', side_effect=[z, np.r_[np.arange(64), 0]]), \
                 patch.object(ev.pd, 'read_csv', side_effect=[candidates, cells]), \
                 patch.object(ev, 'load_npz', return_value=SimpleNamespace(toarray=lambda: fake_counts)), \
                 patch.object(ev, 'thin_full_counts') as thinning:
                with self.assertRaisesRegex(ValueError, 'covariate pool'):
                    ev.main()
                thinning.assert_not_called()


if __name__ == '__main__':
    unittest.main()
