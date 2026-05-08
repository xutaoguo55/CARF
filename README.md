# CARF: Confounder-Adjusted Ranking Framework for Perturbation Inference Benchmarks

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

A general-purpose methodology for revealing the confounding architecture underlying apparent benchmark performance in perturbation inference. CARF operationalizes a falsifiable null hypothesis: **a method's perturbation signal should survive adjustment for measured confounders** (expression level, gene family, co-expression).

## Citation

If you use CARF in your research, please cite:

> Guo X. CARF: A Confounder-Adjusted Ranking Framework for Perturbation Inference Benchmarks — Demonstrated on WWOX Tumor Suppressor Silencing in Diffuse Large B-Cell Lymphoma. Zenodo. DOI: [10.5281/zenodo.20088473](https://doi.org/10.5281/zenodo.20088473)

## Overview

CARF provides seven quantitative metrics and six confounder-adjusted diagnostics:

**Metrics**: PCE (Precision-Coverage Efficiency), PSR (Perturbation Specificity Ratio), EBS (Expression Bias Score), CII (Co-expression Independence Index), CBS (Composite Benchmark Score), Rank Variance Decomposition, Platform Transfer Bound

**Diagnostics**: PSD, EDR, CPS, Statistical Mediation Decomposition, E-value Sensitivity, Conformal Sensitivity

## Quick Start

### Requirements

- R >= 4.2
- Python >= 3.9 (for Geneformer inference only)
- Required R packages: `ggplot2`, `patchwork`, `reshape2`, `scales`, `ggrepel`, `dplyr`
- Geneformer dependencies: see `requirements-carf.txt`

### Reproduce results

```bash
cd code
Rscript 14_generate_figure1.R       # CARF overview
Rscript 09_generate_figures.R       # Main results
Rscript 11_generate_causal_figures.R # Causal framework
Rscript 16_carf_report_card.R       # Report card
Rscript 15_expression_distance_baseline.R  # Expression distance baseline
Rscript 17_embedding_geometry_validation.R # Embedding geometry validation
```

All outputs are generated at 170mm width, 300 DPI, with Okabe-Ito colorblind-friendly palette.

### Reproduce full benchmark

Script execution order is documented in `run_all_figures.R`. Key pipeline steps:

1. `01_prepare_dlbcl_data.R` — Download and preprocess GEO data
2. `02_run_geneformer_perturbation.py` — Geneformer zero-shot perturbation (requires GPU)
3. `03_run_baselines.R` — Pearson, linear model, scTenifoldKnk
4. `04_evaluate_benchmark.R` — Compute PSR, EBS, CII, CBS
5. `05_cross_platform_validation.R` — Cross-platform transfer analysis
6. `06_analyze_failure_modes.R` — Expression bias and embedding analysis
7. `07_compute_benchmark_metrics.R` — Mathematical metrics and CBS_anchor

## Data Availability

All expression data are from public GEO repositories:

| Dataset | Accession | Platform | Use |
|---------|-----------|----------|-----|
| DLBCL discovery | [GSE10846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846) | Affymetrix HG-U133 Plus 2.0 | Discovery (n=420) |
| DLBCL validation | [GSE32918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32918) | Affymetrix HG-U133 Plus 2.0 | Validation (n=172) |
| DLBCL cross-platform | [GSE87371](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87371) | Illumina Human DASL | Cross-platform (n=223) |
| DLBCL scRNA-seq | [GSE182434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182434) | 10x Genomics | scRNA-seq reference |
| DLBCL survival | [GSE31312](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31312) | Affymetrix HG-U133 Plus 2.0 | Survival (n=414) |

Benchmark results (CSV files) are in `benchmark_results/`.

## Repository Structure

```
.
├── code/                          # All analysis scripts
│   ├── 01_prepare_dlbcl_data.R
│   ├── 02_run_geneformer_perturbation.py
│   ├── 03_run_baselines.R
│   ├── 04_evaluate_benchmark.R
│   ├── 05_cross_platform_validation.R
│   ├── 06_analyze_failure_modes.R
│   ├── 07_compute_benchmark_metrics.R
│   ├── 09_generate_figures.R
│   ├── 10_causal_framework.R
│   ├── 11_generate_causal_figures.R
│   ├── 14_generate_figure1.R
│   ├── 15_expression_distance_baseline.R
│   ├── 16_carf_report_card.R
│   ├── 17_embedding_geometry_validation.R
│   ├── 18_compute_cbs_anchor.R
│   ├── 19_geneformer_stability.R
│   └── common_config.R
├── benchmark_results/              # All CSV output files
├── carf_benchmark/                 # Benchmark framework and leaderboard
├── tests/                          # Unit tests
├── Makefile                        # Reproducible pipeline
├── Dockerfile                      # Containerized environment
├── requirements-carf.txt           # Python dependencies
└── README.md
```

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/) (CC BY-NC 4.0).

**In short:**
- ✅ You may share, copy, and redistribute the material in any medium or format
- ✅ You may adapt, remix, transform, and build upon the material
- ✅ You must give appropriate credit, provide a link to the license, and indicate if changes were made
- ❌ You may not use the material for commercial purposes without permission

For commercial use inquiries, please contact the corresponding author.

## Important Caveats

1. **Single-case calibration**: CARF's metrics and diagnostic thresholds were developed on the WWOX/DLBCL case study. Provisional thresholds are heuristics, not binary pass/fail criteria. Independent meta-benchmark validation is needed before prescriptive use.

2. **Geneformer domain mismatch**: Geneformer V2 was pretrained on single-cell RNA-seq. The WWOX/DLBCL evaluation uses bulk microarray data — this tests generalization, not intended use.

3. **Ground truth circularity**: The 17-gene validation set was derived from a pipeline in which Pearson correlation with WWOX was a primary screening step. Performance metrics for correlation-based methods are structurally advantaged.

4. **CBS reference-set dependence**: CBS uses min-max normalization across methods. Adding a new method changes all CBS values. Use `CBS_anchor` (in `benchmark_results/benchmark_cbs_anchor.csv`) for cross-benchmark comparisons.

5. **Geneformer 50-cell subset**: Primary Geneformer results use 50 of 420 samples. Cross-subset stability analysis confirms rank ordering is stable (ρ=0.906 between 10-cell and 50-cell subsets), but PSR shows substantial bootstrap variance (95% CI: [27.74, 127.65] at k=10).

## Contributing

This is research code accompanying a journal submission. Bug reports and suggestions are welcome via GitHub Issues. For substantive methodological contributions, please contact the authors.

## Contact

Xutao Guo — gxt827@126.com

Department of Hematology, Nanfang Hospital, Southern Medical University, Guangzhou, China
