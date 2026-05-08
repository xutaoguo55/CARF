# Reproducibility Guide: Foundation Model Benchmark for Perturbation Inference

## Overview

This repository reproduces all analyses in the manuscript "Benchmarking computational methods for inferring transcriptional consequences of tumor suppressor gene silencing." Four methods are compared for inferring transcriptional consequences of WWOX silencing in DLBCL against 17 experimentally validated targets.

## Software Requirements

### R (v4.4+)
Required packages with tested versions:
```
 install.packages(c("ggplot2", "patchwork", "reshape2", "scales", "mediation", "testthat"))
# Tested versions: ggplot2_3.5.1, patchwork_1.2.0, reshape2_1.4.4,
# scales_1.3.0, mediation_4.5.0
```
For exact reproducibility, a `renv.lock` file is provided in the repository root.

### Python (v3.10+)
Required for Geneformer perturbation inference:
```
pip install transformer-geneformer==0.1.2
pip install torch==2.4.0 scanpy==1.10.1 anndata==0.10.6
pip install transformers datasets pandas scipy
```
Note: Geneformer's original code hardcodes CUDA device placement. Apple Silicon (MPS) compatibility requires patching `emb_extractor.py`, `perturber_utils.py`, and `in_silico_perturber.py` in the Geneformer package. Patch files are provided in `patches/` with instructions in `patches/README_MPS.md`.

## Data Sources

| Data | Source | Accession |
|------|--------|-----------|
| DLBCL expression (discovery) | GEO | GSE10846 |
| DLBCL expression (validation, Affymetrix) | GEO | GSE32918 |
| DLBCL expression (validation, Illumina) | GEO | GSE87371 |
| Geneformer V2 model | HuggingFace | ctheodoris/Geneformer |
| Validated WWOX targets | Literature curation | 29 genes (17 validated) |

Pre-computed benchmark results are included in `benchmark_results/`. To regenerate from raw data, follow the pipeline below.

## Pipeline Execution Order

Scripts in the root directory (`01_` – `06_`) are the original analysis pipeline. Scripts in `code/` contain benchmark computation, figure generation, sensitivity analysis, raw-model diagnostics, and integrity checks.

### Step 1: Data Preparation
```
Rscript 01_prepare_dlbcl_data.R
```
**Inputs**: Raw GSE10846 expression data, ground truth gene list
**Outputs**: `benchmark_results/GSE10846_gene_expression_log2.csv`, `benchmark_results/pseudo_counts.csv`

### Step 2: Geneformer Perturbation Inference
```
python 02_run_geneformer_perturbation.py
```
**Inputs**: Preprocessed expression data, Geneformer model checkpoint
**Outputs**: `benchmark_results/benchmark_geneformer_10cell.csv`, `benchmark_results/benchmark_geneformer_50cell.csv`, `benchmark_results/geneformer_10cell_ranking.csv`

Requires GPU with ≥12GB VRAM (or CPU with patience). The 50-cell perturbation embedding uses 50 random cells per donor × 420 donors = 21,000 embeddings.

### Step 3: Baseline Methods
```
Rscript 03_run_baselines.R
```
**Inputs**: Expression data from Step 1
**Outputs**: `benchmark_results/baseline_pearson.csv`, `benchmark_results/baseline_lm.csv`

### Step 4: Benchmark Evaluation
```
Rscript 04_evaluate_benchmark.R
```
**Inputs**: All method outputs (Steps 2–3)
**Outputs**: `benchmark_results/benchmark_all_methods.csv`, `benchmark_results/benchmark_29genes.csv`, `benchmark_results/ground_truth_29.csv`

### Step 5: Cross-Platform Validation
```
Rscript 05_cross_platform_validation.R
```
**Inputs**: GSE32918, GSE87371 expression matrices
**Outputs**: `benchmark_results/cross_platform_validation.csv`

### Step 6: Failure Mode Analysis
```
Rscript 06_analyze_failure_modes.R
```
**Inputs**: All benchmark results from Steps 1–4
**Outputs**: `benchmark_results/failure_mode_analysis.csv`

### Step 7: Compute Benchmark Metrics
```
Rscript code/07_compute_benchmark_metrics.R
```
**Inputs**: All method outputs
**Outputs**: `benchmark_results/benchmark_mathematical_metrics.csv`, `benchmark_results/benchmark_metrics_summary.csv`, `benchmark_results/benchmark_psr_curves.csv`

### Step 8: Methodological Innovations
```
Rscript code/08_methodological_innovations.R
```
**Inputs**: All benchmark results
**Outputs**: `benchmark_results/benchmark_psr_decomposition.csv`, `benchmark_results/benchmark_expression_deconfounded.csv`, `benchmark_results/benchmark_causal_perturbation_signal.csv`, `benchmark_results/benchmark_bootstrap_ci.csv`

### Step 9: Causal Framework
```
Rscript code/10_causal_framework.R
```
**Inputs**: All benchmark results
**Outputs**: `benchmark_results/benchmark_causal_mediation.csv`, `benchmark_results/benchmark_evalue_sensitivity.csv`, `benchmark_results/benchmark_conformal_psr.csv`

### Step 10: Generate Figures
```
Rscript code/09_generate_figures.R        # Figures 2-7
Rscript code/11_generate_causal_figures.R # Figure 8
```
**Outputs**: `figures_gigascience/figure2_*.pdf` – `figures_gigascience/figure8_*.pdf`

### Step 11: Generate Supplementary Tables
```
Rscript code/12_generate_supplementary_tables.R
```
**Outputs**: core supplementary tables. Re-run this script after Step 15 to copy raw embedding/attention outputs into Tables S19-S20.

### Step 12: Editor Revision Analysis (Symmetric Diagnostics)
```
Rscript code/13_editor_revision_analysis.R
```
**Inputs**: All benchmark results
**Outputs**: `benchmark_results/benchmark_rank_stability.csv`, `benchmark_results/benchmark_scTenifoldKnk_edr.csv`, `benchmark_results/benchmark_scTenifoldKnk_cps_decomposition.csv`, `benchmark_results/benchmark_scTenifoldKnk_evalue.csv`

### Step 13: Generate Figure 1 (Study Design Schematic)
```
Rscript code/14_generate_figure1.R
```
**Outputs**: `figures_gigascience/figure1_carf_overview.pdf`

### Step 14: Additional Sensitivity Figures and Report Card
```
Rscript code/15_expression_distance_baseline.R
Rscript code/16_carf_report_card.R
Rscript code/17_embedding_geometry_validation.R
Rscript code/run_all_figures.R
```
**Outputs**: Supplementary expression-distance/geometry diagnostics and `figures_gigascience/figure9_carf_report_card.pdf`

### Step 15: Raw Geneformer Hidden-State Density and Attention Diagnostics
```
python3 code/19_raw_embedding_attention_analysis.py --preflight
python3 code/19_raw_embedding_attention_analysis.py --max-cells 8 --attention-cells 4 --attention-window 512 --density-k 20 --device cpu
```
**Inputs**: Local Geneformer model checkpoint, token dictionary, gene-name dictionary, and tokenized `.dataset` directory. These can be supplied through command-line flags or `GENEFORMER_MODEL_DIR`, `GENEFORMER_DATASET_DIR`, `GENEFORMER_TOKEN_DICTIONARY`, and `GENEFORMER_GENE_NAME_DICTIONARY`.
**Outputs**: `benchmark_results/benchmark_raw_embedding_attention.csv`, `benchmark_results/benchmark_raw_embedding_attention_summary.csv`, `benchmark_results/benchmark_raw_embedding_attention_metadata.json`

The script computes full-sequence final-layer hidden-state density for sampled cells. Attention is summarized in a WWOX-centered rank-order window because retaining full 3,992-token attention tensors across all layers/cells is quadratic in sequence length.

### Step 16: Final Supplementary Tables, Unit Tests, and Integrity Verification
```
Rscript --vanilla code/12_generate_supplementary_tables.R
Rscript --vanilla tests/testthat.R
Rscript --vanilla code/18_generate_manifest.R
Rscript --vanilla code/18_generate_manifest.R --verify
```
**Outputs**: finalized `supplementary_tables/Table_S1_*.csv` – `Table_S20_*.csv`, `MANIFEST.sha256.csv`, `MANIFEST.sha256.txt`, and `MANIFEST.sha256.verify.csv`

### Step 17: CARF-Benchmark v1 Standardized Leaderboard
```
make carf-v1
make carf-v1-sources
make carf-v1-adapter-smoke
```
Equivalent explicit commands:
```
Rscript --vanilla code/20_prepare_carf_benchmark_v1.R
python3 carf_benchmark/scripts/validate_schema.py --schema carf_benchmark/schema/method_scores.schema.json --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/method_scores.csv
python3 carf_benchmark/scripts/validate_schema.py --schema carf_benchmark/schema/ground_truth.schema.json --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/ground_truth.csv
python3 carf_benchmark/scripts/validate_schema.py --schema carf_benchmark/schema/covariates.schema.json --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/covariates.csv
Rscript --vanilla code/21_run_carf_benchmark_v1.R
python3 carf_benchmark/datasets/prepare_public_perturbseq.py materialize-sources --dataset-id all
python3 carf_benchmark/tests/test_adapters.py
```
**Outputs**: standardized seed inputs in `carf_benchmark/runs/wwox_dlbcl_v1/inputs/` and leaderboard outputs in `carf_benchmark/leaderboard/`, including run-level leaderboard, model-level summary with bootstrap intervals, audit matrix, PSR curves, and dataset readiness.

CARF-Benchmark v1 recomputes metrics from full standardized method-score files. This makes it suitable as a reusable benchmark standard; values can differ from legacy manuscript tables that used curated 29-gene comparison tables for selected case-study diagnostics.

### Step 18: Public Perturb-seq Dataset Activation and Model Adapters
```
python3 carf_benchmark/datasets/prepare_public_perturbseq.py list-sources
CARF_DATASETS=adamson_2016_perturbseq make carf-v1-download-datasets
python3 carf_benchmark/datasets/prepare_public_perturbseq.py prepare \
  --dataset-id adamson_2016_perturbseq \
  --h5ad carf_benchmark/raw/adamson_2016_perturbseq/Adamson.h5ad \
  --max-perturbations 3
```
**Outputs**: public source manifests under `carf_benchmark/runs/<dataset_id>/source_manifest.json`; after h5ad download/preparation, full CARF inputs under `carf_benchmark/runs/<dataset_id>/inputs/`.

This repository includes one prepared public Perturb-seq smoke-scale run:
`replogle_2022_genome_scale`, using three perturbations from
`Replogle_exp6.h5ad.gz`. It includes zero-delta, mean-expression, and target
co-expression baselines. The raw h5ad is excluded from MANIFEST and Docker
context; the standardized CARF inputs are included.

Native scGPT, scFoundation, UCE, scBERT, GEARS, and CPA outputs are standardized with:
```
python3 carf_benchmark/adapters/run_adapter.py gene-score-csv --help
python3 carf_benchmark/adapters/run_adapter.py expression-delta --help
```
The adapter layer converts real native prediction files into `method_scores.csv`; it does not generate synthetic model predictions when external model outputs are absent.

## Random Seed Documentation

All analyses with stochastic components use fixed random seeds for reproducibility:
- Bootstrap resampling (PSR/CBS 95% CI): `set.seed(42)`
- Cross-platform ρ bootstrap CI: `set.seed(42)`
- Split-conformal prediction (training/calibration split): `set.seed(42)`
- scTenifoldKnk tensor decomposition: deterministic (no random component)
- Geneformer perturbation inference: deterministic given the model checkpoint and input data

## Geneformer MPS Compatibility

The 50-cell Geneformer perturbation run was performed on Apple Silicon (MPS). Geneformer's original code hardcodes CUDA device placement. To reproduce on MPS, apply the patches provided in `patches/`:
1. `emb_extractor.py`: auto-detect device (CUDA > MPS > CPU)
2. `perturber_utils.py`: recognize MPS as valid accelerator
3. `in_silico_perturber.py`: set multiprocess start method to "fork" (Python 3.13+)

On MPS, each 3,992-token sample processes in ~12 minutes. On CUDA GPU (T4), approximately 1-2 minutes per sample.

## Expected Outputs Summary

| Category | Files | Description |
|----------|-------|-------------|
| Figures | 12 PDFs (Figures 1-9 + Supplementary Figures S1-S2) | Main and supplementary figures |
| Main results | 15 CSVs | Benchmark metrics, PSR curves, EDR, CPS |
| Causal results | 3 CSVs | Mediation, E-values, Conformal prediction |
| Supplementary | 22 CSVs | Tables S1–S20, including split S6 and S17 files |

## Key Results to Verify

After running the full pipeline, these canonical values should match:

| Metric | Value | File |
|--------|-------|------|
| Geneformer CBS | 0.794 | `benchmark_mathematical_metrics.csv` |
| Geneformer PSR k=10 | 70.39 | `benchmark_psr_curves.csv` |
| scTenifoldKnk PSR k=10 (full HVG scope) | 23.53 | `benchmark_psr_curves.csv` |
| scTenifoldKnk PSR k=10 (EDR/CPS subset) | 10.54 | Table 10 in manuscript |
| Geneformer EBS | 0.617 | `benchmark_mathematical_metrics.csv` |
| Raw density vs \|cosine shift\| | ρ=0.150 | `benchmark_raw_embedding_attention_summary.csv` |
| Raw embedding norm vs \|cosine shift\| | ρ=-0.516 | `benchmark_raw_embedding_attention_summary.csv` |
| gene→WWOX all-layer attention vs \|cosine shift\| | ρ=0.326 | `benchmark_raw_embedding_attention_summary.csv` |
| Expression variance (R²) | 45.7% | `benchmark_causal_perturbation_signal.csv` |
| EDR PSR k=10 (deconfounded) | 0.00 | Section 3.8 of manuscript |
| CPS Mann-Whitney p | 0.39 | `benchmark_causal_perturbation_signal.csv` |
| Causal mediation: expression | 63.1% | `benchmark_causal_mediation.csv` |
| E-value (GF PSR original) | 140.3 | `benchmark_evalue_sensitivity.csv` |

## Runtime Estimates

| Step | Hardware | Approximate Time |
|------|----------|-----------------|
| Step 1 (data prep) | Any | < 5 min |
| Step 2 (Geneformer) | GPU 12GB | ~8 hours |
| Step 3 (baselines) | Any | < 10 min |
| Steps 4-14 | Any | < 1 hour total |
| Raw hidden-state/attention diagnostics | CPU, sampled cells | ~1 min for 8 embedding cells + 4 attention windows |
| CARF-Benchmark v1 seed leaderboard | Any | < 1 min |
| Public Perturb-seq source manifests | Any | < 1 min |
| Public h5ad download/preparation | Network + RAM dependent | minutes to hours depending on selected datasets |

## Directory Structure After Full Run

```
foundation_model_benchmark/
├── 01_prepare_dlbcl_data.R
├── 02_run_geneformer_perturbation.py
├── 03_run_baselines.R
├── 04_evaluate_benchmark.R
├── 05_cross_platform_validation.R
├── 06_analyze_failure_modes.R
├── code/
│   ├── 07_compute_benchmark_metrics.R
│   ├── 08_methodological_innovations.R
│   ├── 09_generate_figures.R
│   ├── 10_causal_framework.R
│   ├── 11_generate_causal_figures.R
│   ├── 12_generate_supplementary_tables.R
│   ├── 18_generate_manifest.R
│   ├── 19_raw_embedding_attention_analysis.py
│   ├── 20_prepare_carf_benchmark_v1.R
│   └── 21_run_carf_benchmark_v1.R
├── carf_benchmark/
│   ├── schema/
│   ├── registry/
│   ├── datasets/
│   ├── adapters/
│   ├── runs/
│   └── leaderboard/
├── R/
│   └── carf_metrics.R
├── tests/
│   └── testthat/
├── benchmark_results/
│   ├── *.csv (benchmark result files)
│   └── GSE10846_filtered_v2.loom
├── figures_gigascience/
│   └── *.pdf (main and supplementary figures)
├── supplementary_tables/
│   └── Table_S*.csv (S1-S20)
├── MANIFEST.sha256.csv
├── MANIFEST.sha256.txt
├── manuscript/
│   └── manuscript_draft.md
├── data/
│   └── tmp/
└── README_reproducibility.md (this file)
```

## Contact

For questions about reproducibility, contact the corresponding author.
