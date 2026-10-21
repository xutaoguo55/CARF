import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from evaluate_pearson_measurement import normalized_model_expression,absolute_pearson


class PearsonMeasurementChecks(unittest.TestCase):
    def test_absolute_correlation_and_zero_variance_policy(self):
        expression=np.array([[1.,2.,3.,4.],[4.,3.,2.,1.],[2.,2.,2.,2.],[1.,4.,2.,3.]])
        score,zero_targets,zero_genes=absolute_pearson(expression,[0,2])
        np.testing.assert_allclose(score[0,:2],[1.,1.])
        np.testing.assert_array_equal(score[1],0)
        np.testing.assert_array_equal(score[:,2],0)
        np.testing.assert_array_equal(zero_targets,[False,True])
        np.testing.assert_array_equal(zero_genes,[False,False,True,False])
        self.assertAlmostEqual(score[0,3],abs(np.corrcoef(expression[0],expression[3])[0,1]))

    def test_all_model_gene_denominator_and_zero_cell_failure(self):
        expression,totals=normalized_model_expression([[1,2],[9,18],[90,180]])
        np.testing.assert_allclose(expression[:,0],np.log1p([100,900,9000]))
        np.testing.assert_allclose(expression[:,0],expression[:,1])
        np.testing.assert_array_equal(totals,[100,200])
        with self.assertRaises(ValueError):normalized_model_expression([[0,1],[0,2]])
        with self.assertRaises(ValueError):normalized_model_expression([[.5,1],[1,2]])

if __name__=="__main__":unittest.main()
