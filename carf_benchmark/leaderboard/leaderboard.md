# CARF-Benchmark v1 Leaderboard

This leaderboard is generated from standardized CARF-Benchmark v1 inputs.
CBS is normalized within each dataset/perturbation task; CBS_anchored uses random/perfect anchors and should be interpreted alongside audit diagnostics.

| rank_overall | dataset_id | perturbation_id | model_name | model_family | coverage | PCE | PSR_at_10 | PSR_at_10_p_perm | PSR_max | EBS | CII | CBS | CBS_anchored |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.000 | norman_2019_combo | KLF1 | Control co-expression baseline | statistical_baseline | 1.000 | 2.619 | 6.579 | 0.000 | 6.579 | 0.497 | 0.000 | 0.694 | 0.487 |
| 2.000 | replogle_2022_genome_scale | FDPS | Control co-expression baseline | statistical_baseline | 1.000 | 2.445 | 4.377 | 0.001 | 4.377 | 0.377 | 0.000 | 0.687 | 0.464 |
| 3.000 | wwox_dlbcl_v1 | WWOX_silencing | Geneformer V2 50-cell | foundation_model | 0.184 | 3.265 | 70.394 | 0.000 | 70.394 | 0.617 | 0.898 | 0.684 | 0.202 |
| 4.000 | adamson_2016_perturbseq | DNAJC19 | Control co-expression baseline | statistical_baseline | 1.000 | 2.356 | 2.545 | 0.006 | 2.545 | 0.381 | 0.000 | 0.680 | 0.548 |
| 5.000 | replogle_2022_genome_scale | FDPS+HUS1 | Control co-expression baseline | statistical_baseline | 1.000 | 2.517 | 3.387 | 0.009 | 3.387 | 0.415 | 0.000 | 0.676 | 0.424 |
| 6.000 | wwox_dlbcl_v1 | WWOX_silencing | scTenifoldKnk | network_model | 0.092 | 4.874 | 11.765 | 0.033 | 11.765 | 0.148 | 0.916 | 0.673 | 0.314 |
| 7.000 | adamson_2016_perturbseq | ASCC3 | Mean expression baseline | statistical_baseline | 1.000 | 4.783 | 6.024 | 0.000 | 6.885 | 1.000 | 0.780 | 0.667 | 0.359 |
| 8.000 | adamson_2016_perturbseq | SEC61B | Mean expression baseline | statistical_baseline | 1.000 | 9.310 | 13.630 | 0.002 | 13.630 | 1.000 | 0.843 | 0.667 | 0.134 |
| 9.000 | adamson_2016_perturbseq | DNAJC19 | Mean expression baseline | statistical_baseline | 1.000 | 3.659 | 3.272 | 0.000 | 3.272 | 1.000 | 0.829 | 0.667 | 0.468 |
| 10.000 | norman_2019_combo | CEBPE | Mean expression baseline | statistical_baseline | 1.000 | 3.372 | 2.628 | 0.000 | 2.647 | 1.000 | 0.610 | 0.667 | 0.495 |
| 11.000 | replogle_2022_genome_scale | FDPS+HUS1 | Mean expression baseline | statistical_baseline | 1.000 | 3.367 | 4.741 | 0.000 | 4.877 | 1.000 | 0.863 | 0.667 | 0.323 |
| 12.000 | replogle_2022_genome_scale | HUS1 | Mean expression baseline | statistical_baseline | 1.000 | 6.428 | 83.650 | 0.012 | 83.650 | 1.000 | 0.898 | 0.667 | 0.057 |
| 13.000 | replogle_2022_genome_scale | FDPS | Mean expression baseline | statistical_baseline | 1.000 | 3.559 | 5.836 | 0.000 | 5.836 | 1.000 | 0.919 | 0.667 | 0.348 |
| 14.000 | norman_2019_combo | CEBPE | Control co-expression baseline | statistical_baseline | 1.000 | 2.957 | 1.460 | 0.236 | 1.985 | 0.663 | 0.000 | 0.655 | 0.508 |
| 15.000 | norman_2019_combo | KLF1 | Mean expression baseline | statistical_baseline | 1.000 | 4.527 | 4.935 | 0.001 | 6.250 | 1.000 | 0.910 | 0.650 | 0.345 |
| 16.000 | replogle_2022_genome_scale | HUS1 | Control co-expression baseline | statistical_baseline | 1.000 | 1.943 | 83.650 | 0.013 | 83.650 | 0.391 | 0.000 | 0.637 | 0.259 |
| 17.000 | norman_2019_combo | BAK1 | Control co-expression baseline | statistical_baseline | 1.000 | 2.144 | 50.250 | 0.019 | 50.250 | 0.459 | 0.000 | 0.626 | 0.215 |
| 18.000 | norman_2019_combo | BAK1 | GEARS | deep_perturbation_model | 0.371 | 6.378 | 0.000 | 1.000 | 22.356 | 0.671 | 0.945 | 0.591 | 0.157 |
| 19.000 | adamson_2016_perturbseq | ASCC3 | Control co-expression baseline | statistical_baseline | 1.000 | 2.484 | 3.442 | 0.019 | 4.131 | 0.418 | 0.000 | 0.567 | 0.402 |
| 20.000 | adamson_2016_perturbseq | SEC61B | Control co-expression baseline | statistical_baseline | 1.000 | 2.546 | 9.086 | 0.019 | 9.086 | 0.371 | 0.000 | 0.523 | 0.286 |
| 21.000 | norman_2019_combo | CEBPE | GEARS | deep_perturbation_model | 0.371 | 1.931 | 1.083 | 0.001 | 1.083 | 0.521 | 0.603 | 0.487 | 0.587 |
| 22.000 | wwox_dlbcl_v1 | WWOX_silencing | Linear model | statistical_baseline | 1.000 | 1.581 | 0.000 | 1.000 | 0.000 | 0.036 | 0.129 | 0.430 | 0.322 |
| 23.000 | norman_2019_combo | KLF1 | GEARS | deep_perturbation_model | 0.371 | 2.216 | 2.744 | 0.000 | 2.744 | 0.671 | 0.938 | 0.412 | 0.423 |
| 24.000 | wwox_dlbcl_v1 | WWOX_silencing | Pearson correlation | statistical_baseline | 1.000 | 1.581 | 0.000 | 1.000 | 0.000 | 0.314 | 0.000 | 0.337 | 0.229 |
| 25.000 | adamson_2016_perturbseq | ASCC3 | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 26.000 | adamson_2016_perturbseq | SEC61B | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 27.000 | adamson_2016_perturbseq | DNAJC19 | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 28.000 | norman_2019_combo | KLF1 | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 29.000 | norman_2019_combo | BAK1 | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 30.000 | norman_2019_combo | CEBPE | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 31.000 | replogle_2022_genome_scale | FDPS+HUS1 | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 32.000 | replogle_2022_genome_scale | HUS1 | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 33.000 | replogle_2022_genome_scale | FDPS | Mean prediction baseline | statistical_baseline | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.333 | 0.333 |
| 34.000 | norman_2019_combo | BAK1 | Mean expression baseline | statistical_baseline | 1.000 | 4.991 | 0.000 | 1.000 | 10.050 | 1.000 | 0.980 | 0.328 | 0.008 |

Active adapter converters for scGPT, scFoundation, UCE, scBERT, GEARS, CPA, and mean baselines are listed in `carf_benchmark/registry/models.csv`.
