import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from evaluate_heldout_linear import stable_folds,bilinear_components,bilinear_solve,predict_selected_fold,within_fold_utilities


class HeldoutLinearChecks(unittest.TestCase):
    def test_heldout_outcome_change_does_not_change_fit(self):
        rng=np.random.default_rng(15);Y=rng.normal(size=(12,8))
        train=np.arange(5);test=np.arange(5,8);target_gene=np.arange(8)
        p,_=predict_selected_fold(Y,train,test,target_gene,k=2)
        changed=Y.copy();changed[:,test]=rng.normal(1000,500,size=(12,3))
        q,_=predict_selected_fold(changed,train,test,target_gene,k=2)
        np.testing.assert_array_equal(p,q)

    def test_two_sided_solve_matches_kronecker_normal_equations(self):
        rng=np.random.default_rng(3);Y=rng.normal(size=(9,6));G,b,P,W,_=bilinear_components(Y,np.arange(6),k=2,ridge=.1)
        left=G.T@G+.1*np.eye(2);right=P.T@P+.1*np.eye(2);rhs=G.T@(Y-b)@P
        direct=np.linalg.solve(np.kron(right.T,left),rhs.flatten(order="F")).reshape((2,2),order="F")
        np.testing.assert_allclose(W,direct,rtol=1e-12,atol=1e-12)

    def test_pca_column_centering_and_score_scaling(self):
        Y=np.arange(48,dtype=float).reshape(8,6)+np.random.default_rng(2).normal(size=(8,6))
        G,_,_,_,s=bilinear_components(Y,np.arange(6),k=2)
        np.testing.assert_allclose(G.mean(axis=0),0,atol=1e-12)
        np.testing.assert_allclose(G.T@G,np.diag(s*s),rtol=1e-12,atol=1e-12)

    def test_fold_invariant_mean_has_zero_identity_and_no_crossfold_values(self):
        scores=np.array([[1.,2,3],[1.,2,3],[3.,1,2],[3.,1,2]])
        B=np.arange(12,dtype=float).reshape(4,3)/10;folds=np.array([0,0,1,1])
        _,u,same,other=within_fold_utilities(scores,B,B,folds,2)
        np.testing.assert_allclose(same-other,0,atol=1e-15)
        self.assertTrue(np.isnan(u[:2,2:]).all());self.assertTrue(np.isnan(u[2:,:2]).all())

    def test_fold_assignment_stable_under_input_order(self):
        targets=[f"g{i}" for i in range(14)];a=dict(zip(targets,stable_folds(targets,"ds")))
        reverse=targets[::-1];b=dict(zip(reverse,stable_folds(reverse,"ds")))
        self.assertEqual(a,b);self.assertEqual(set(a.values()),{0,1,2})

if __name__=="__main__":unittest.main()
