# CARF-Benchmark v1 Model Summary

Model-level summaries aggregate complete dataset/perturbation runs.
CBS_mean averages task-normalized CBS values; CBS_anchored_mean reports random/perfect anchored scores.
Intervals are non-parametric bootstrap 95% intervals across available runs.

| rank_model | model_name | model_family | n_runs | n_datasets | CBS_mean | CBS_ci_low | CBS_ci_high | CBS_anchored_mean | PCE_mean | PSR_max_mean | PSR_at_10_p_perm_min | PSR_max_p_perm_min | EBS_mean | CII_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.000 | Geneformer V2 50-cell | foundation_model | 1.000 | 1.000 | 0.684 | 0.684 | 0.684 | 0.202 | 3.265 | 70.394 | 0.000 | 0.000 | 0.617 | 0.898 |
| 2.000 | scTenifoldKnk | network_model | 1.000 | 1.000 | 0.673 | 0.673 | 0.673 | 0.314 | 4.874 | 11.765 | 0.033 | 0.007 | 0.148 | 0.916 |
| 3.000 | Control co-expression baseline | statistical_baseline | 9.000 | 3.000 | 0.638 | 0.601 | 0.669 | 0.399 | 2.445 | 18.443 | 0.000 | 0.000 | 0.441 | 0.000 |
| 4.000 | Mean expression baseline | statistical_baseline | 9.000 | 3.000 | 0.627 | 0.550 | 0.667 | 0.282 | 4.889 | 15.233 | 0.000 | 0.000 | 1.000 | 0.848 |
| 5.000 | GEARS | deep_perturbation_model | 3.000 | 1.000 | 0.497 | 0.412 | 0.591 | 0.389 | 3.509 | 8.728 | 0.000 | 0.000 | 0.621 | 0.829 |
| 6.000 | Linear model | statistical_baseline | 1.000 | 1.000 | 0.430 | 0.430 | 0.430 | 0.322 | 1.581 | 0.000 | 1.000 | 1.000 | 0.036 | 0.129 |
| 7.000 | Pearson correlation | statistical_baseline | 1.000 | 1.000 | 0.337 | 0.337 | 0.337 | 0.229 | 1.581 | 0.000 | 1.000 | 1.000 | 0.314 | 0.000 |
| 8.000 | Mean prediction baseline | statistical_baseline | 9.000 | 3.000 | 0.333 | 0.333 | 0.333 | 0.333 | 0.000 | 0.000 | 1.000 | 1.000 | 0.000 | 1.000 |
