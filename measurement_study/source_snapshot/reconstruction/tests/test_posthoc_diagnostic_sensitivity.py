"""Synthetic mathematical checks for the post-hoc diagnostic comparison only."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_posthoc_diagnostic_sensitivity import (percentile_midrank,ranking_scores,
    expression_strata,match_background,evaluate_weights,strict_mean)
from carf_core import topk_weights


class PosthocDiagnosticSensitivityTests(unittest.TestCase):
    def test_midranks_preserve_ties_without_id_breaking(self):
        np.testing.assert_array_equal(percentile_midrank([1,1,3,2]),[.25,.25,.875,.625])
        bins,_=expression_strata([1,1,1,2,2,3])
        self.assertEqual(len(set(bins[:3])),1)
        self.assertEqual(bins[3],bins[4])

    def test_rank_OLS_invariant_under_strictly_monotone_transform(self):
        raw=np.array([[1.,2,2,5,4,3],[4.,3,1,2,5,6]])
        cov=np.array([.1,.5,.2,.6,.3,.4])
        a=ranking_scores(raw,cov,'rank_OLS')
        b=ranking_scores(np.log1p(raw),np.exp(cov),'rank_OLS')
        np.testing.assert_array_equal(a,b)

    def test_rank_OLS_constant_scores_and_constant_expression(self):
        cov=np.arange(10.)
        np.testing.assert_array_equal(ranking_scores(np.ones((2,10)),cov,'rank_OLS'),np.zeros((2,10)))
        raw=np.array([[4.,1,1,3,2]])
        x=percentile_midrank(raw[0])
        np.testing.assert_allclose(ranking_scores(raw,np.ones(5),'rank_OLS')[0],x-x.mean(),atol=1e-15)
        np.testing.assert_array_equal(ranking_scores(cov[None,:],cov,'rank_OLS'),np.zeros((1,10)))

    def test_fractional_boundary_ties_have_total_k(self):
        w=topk_weights([5,4,4,4,1],2)
        np.testing.assert_allclose(w,[1,1/3,1/3,1/3,0])
        self.assertAlmostEqual(w.sum(),2)

    def test_matching_preserves_stratum_mass_and_complement_capacity(self):
        cov=np.arange(20.)
        w=np.zeros((1,20));w[0,[0,4,8,12,16]]=1
        v,supported,info,W,H,bins,count=match_background(w,cov,k=5)
        self.assertTrue(supported[0])
        self.assertAlmostEqual(v.sum(),5)
        self.assertTrue(np.all(w+v<=1+1e-12))
        for b in range(5):self.assertAlmostEqual(v[0,bins==b].sum(),W[0,b])
        np.testing.assert_allclose(info['expression_mean_difference'],[(w-v)@cov/5][0])

    def test_constant_fractional_weights_match_themselves(self):
        w=np.full((2,23),5/23)
        v,ok,*_=match_background(w,np.arange(23.),k=5)
        self.assertTrue(ok.all());np.testing.assert_allclose(v,w,atol=1e-15)

    def test_capacity_failure_is_whole_ranking_NA(self):
        w=np.zeros((2,10));w[0,:2]=1;w[1,[0,4]]=1
        v,ok,info,W,H,bins,count=match_background(w,np.arange(10.),k=2)
        self.assertFalse(ok[0]);self.assertTrue(ok[1]);self.assertTrue(np.isnan(v[0]).all())
        self.assertEqual(info['matching_failed_strata'][0],1)
        self.assertEqual(info['matching_max_capacity_deficit'][0],2)

    def test_empty_strata_and_constant_expression_do_not_force_failure(self):
        w=np.zeros((1,10));w[0,:2]=1
        v,ok,_,_,_,bins,count=match_background(w,np.ones(10),k=2)
        self.assertTrue(ok[0]);self.assertEqual(np.count_nonzero(count),1)
        np.testing.assert_allclose(v[0,2:],.25)

    def test_failure_invalidates_fold_identity_without_shrinking_sources(self):
        rng=np.random.default_rng(9)
        w=np.zeros((4,10));w[:,[0,4]]=1;w[0]=0;w[0,:2]=1
        v,ok,*_=match_background(w,np.arange(10.),k=2)
        f=evaluate_weights(w,v,ok,rng.normal(size=(4,10)),rng.normal(size=(4,10)),np.zeros((3,10)),np.array([0,0,1,1]),k=2)
        self.assertTrue(np.isnan(f['matched_identity_gain'][:2]).all())
        self.assertTrue(np.isfinite(f['matched_total_BC'][1]))
        self.assertTrue(np.isfinite(f['matched_identity_gain'][2:]).all())
        np.testing.assert_array_equal(f['matching_fold_supported'],[False,False,True,True])

    def test_same_fold_shared_ranking_zero_identity(self):
        rng=np.random.default_rng(7)
        w=np.zeros((4,20));w[:2,[0,4]]=1;w[2:,[8,12]]=1
        v,ok,*_=match_background(w,np.arange(20.),k=2)
        f=evaluate_weights(w,v,ok,rng.normal(size=(4,20)),rng.normal(size=(4,20)),rng.normal(size=(3,20)),np.array([0,0,1,1]),k=2)
        np.testing.assert_allclose(f['identity_gain'],0,atol=1e-15)
        np.testing.assert_allclose(f['matched_identity_gain'],0,atol=1e-15)

    def test_weighted_utility_matches_explicit_fixed_source_loops(self):
        rng=np.random.default_rng(4);B=rng.normal(size=(6,30));C=rng.normal(size=B.shape)
        w=np.vstack([topk_weights(rng.normal(size=30),3) for _ in range(6)])
        v,ok,*_=match_background(w,np.ones(30),k=3)
        folds=np.array([0,0,0,1,1,1]);f=evaluate_weights(w,v,ok,B,C,np.zeros((2,30)),folds,k=3)
        for t in range(6):
            sources=[r for r in range(6) if folds[r]==folds[t] and r!=t]
            own=np.dot(w[t]-v[t],B[t]*C[t])/3
            other=np.mean([np.dot(w[r]-v[r],B[t]*C[t])/3 for r in sources])
            self.assertAlmostEqual(f['matched_identity_gain'][t],own-other)

    def test_raw_C_only_invariance_does_not_imply_matched_invariance(self):
        raw=np.arange(20.)[None,:]
        np.testing.assert_array_equal(ranking_scores(raw,np.arange(20.),'raw'),ranking_scores(raw,np.ones(20),'raw'))
        w=np.zeros((1,20));w[0,[0,4]]=1
        v1=match_background(w,np.arange(20.),k=2)[0]
        v2=match_background(w,np.ones(20),k=2)[0]
        self.assertFalse(np.allclose(v1,v2))

    def test_strict_summary_does_not_drop_failed_values(self):
        self.assertTrue(np.isnan(strict_mean([1.,np.nan])))
        self.assertEqual(strict_mean([1.,3.]),2)


if __name__=='__main__':unittest.main()
