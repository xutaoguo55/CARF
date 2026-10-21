"""Scientific invariants for the fixed-scope factorial pilot."""
import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from measurement_diagnostic_pilot import (thin_counts, stream_rng, covariate_mean,
                                         evaluate_scores, ranking_scores)
from carf_core import topk_weights


class MeasurementDiagnosticChecks(unittest.TestCase):
    def test_thinning_is_reproducible_and_count_conservative(self):
        counts=csr_matrix(np.array([[1,9,0],[2,0,20]]))
        a=thin_counts(counts,.5,stream_rng("test",.5,1,"Dm"))
        b=thin_counts(counts,.5,stream_rng("test",.5,1,"Dm"))
        self.assertEqual((a!=b).nnz,0)
        self.assertTrue(np.all(a.toarray()<=counts.toarray()))
        self.assertEqual((counts!=thin_counts(counts,1,stream_rng("test",1,1,"Dm"))).nnz,0)
        self.assertNotEqual(stream_rng("test",.5,1,"Dm").integers(2**32),
                            stream_rng("test",.5,1,"Dc").integers(2**32))

    def test_covariate_denominator_includes_all_retained_genes(self):
        means,_=covariate_mean(csr_matrix([[1,9,90],[2,18,180]]))
        np.testing.assert_allclose(means,np.log1p([100,900,9000]))
        with self.assertRaises(ValueError): covariate_mean(csr_matrix([[0,0],[1,1]]))

    def test_fixed_scores_untouched_by_covariate_only_intervention(self):
        distance=np.array([[1.,2,3,4],[4.,2,1,3]])
        x=ranking_scores(distance,np.array([1.,2,3,4]))[0]
        y=ranking_scores(distance,np.array([4.,1,3,2]))[0]
        np.testing.assert_array_equal(x["raw_distance"],y["raw_distance"])
        self.assertFalse(np.allclose(x["EDR"],y["EDR"]))

    def test_target_invariant_ranking_has_zero_identity_gain(self):
        scores=np.tile([1.,2,3,4],(3,1))
        b=np.array([[1,0,2,1],[0,2,1,0],[2,2,0,1]],float)
        w=np.vstack([topk_weights(x,2) for x in scores])
        result,_,_,_=evaluate_scores(scores,np.array([1.,2,3,4]),np.log10(scores),
                                     b,b,np.zeros((2,4)),w,["a","b","c"],2)
        np.testing.assert_allclose(result.identity_gain,0,atol=1e-15)
        np.testing.assert_allclose(result.top50_overlap,1)
        np.testing.assert_allclose(result.null_response_mean,0)

if __name__=="__main__": unittest.main()
