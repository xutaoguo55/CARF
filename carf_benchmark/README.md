# CARF-Benchmark v1

CARF-Benchmark v1 is a standardized confounder-audit layer for perturbation
prediction benchmarks. It is designed to sit downstream of any model that
produces gene-level perturbation scores or ranks.

The seed run in this repository uses the WWOX/DLBCL benchmark and four active
methods: Pearson correlation, linear model, scTenifoldKnk, and Geneformer.
The registry now also activates public Perturb-seq data sources and executable
adapter conversion for scGPT, scFoundation, UCE, scBERT, GEARS, CPA, and a
mean-prediction baseline.

## What v1 Standardizes

- `method_scores.csv`: one row per model, perturbation, and scored gene
- `ground_truth.csv`: positive and negative validation labels
- `covariates.csv`: expression, co-expression, and gene-family covariates
- `run_manifest.json`: run metadata and metric settings
- `leaderboard.csv`: recomputed CARF metrics from standardized inputs
- `audit_matrix.csv`: diagnostic flags for expression bias, co-expression dependence, and top-k enrichment

## Seed Benchmark

Build the standardized WWOX/DLBCL seed inputs:

```bash
Rscript --vanilla code/20_prepare_carf_benchmark_v1.R
```

Validate schema conformance:

```bash
python3 carf_benchmark/scripts/validate_schema.py \
  --schema carf_benchmark/schema/method_scores.schema.json \
  --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/method_scores.csv

python3 carf_benchmark/scripts/validate_schema.py \
  --schema carf_benchmark/schema/ground_truth.schema.json \
  --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/ground_truth.csv

python3 carf_benchmark/scripts/validate_schema.py \
  --schema carf_benchmark/schema/covariates.schema.json \
  --csv carf_benchmark/runs/wwox_dlbcl_v1/inputs/covariates.csv
```

Generate the leaderboard:

```bash
Rscript --vanilla code/21_run_carf_benchmark_v1.R
```

Outputs:

- `carf_benchmark/leaderboard/leaderboard.csv`
- `carf_benchmark/leaderboard/leaderboard.md`
- `carf_benchmark/leaderboard/leaderboard_summary.csv`
- `carf_benchmark/leaderboard/leaderboard_summary.md`
- `carf_benchmark/leaderboard/audit_matrix.csv`
- `carf_benchmark/leaderboard/psr_curves.csv`
- `carf_benchmark/leaderboard/dataset_readiness.csv`

## Active Public Perturb-seq Run

The repository also includes a prepared Replogle public Perturb-seq smoke-scale
run generated from `Replogle_exp6.h5ad.gz`:

- dataset: `replogle_2022_genome_scale`
- perturbations: `FDPS+HUS1`, `HUS1`, `FDPS`
- rows per standardized input: 15,057
- baseline models: `mean_baseline`, `mean_expression_baseline`,
  `control_coexpression_baseline`

The raw h5ad lives under `carf_benchmark/raw/` when downloaded and is excluded
from Docker context and MANIFEST checksums; the standardized inputs are tracked.

The expression and co-expression baselines are deliberate negative/sanity
controls. They quantify whether a benchmark can be solved by high abundance or
target co-expression alone, which is essential for interpreting foundation-model
claims.

## Important Metric Note

The v1 seed leaderboard recomputes metrics from complete standardized method
score files. This is the correct behavior for a reusable benchmark standard.
Some values can therefore differ from legacy manuscript tables that used a
curated 29-gene comparison table for selected diagnostics. Treat the v1
leaderboard as the standard-run output and the manuscript tables as the
case-study narrative output.

## Adding a New Dataset

1. Add the dataset to `carf_benchmark/registry/datasets.csv` and, for public
   h5ad sources, `carf_benchmark/registry/public_dataset_sources.csv`.
2. Materialize source manifests:

```bash
make carf-v1-sources
```

3. Download a public h5ad when storage/network budget allows:

```bash
CARF_DATASETS=adamson_2016_perturbseq make carf-v1-download-datasets
```

4. Prepare CARF inputs from the downloaded h5ad:

```bash
python3 carf_benchmark/datasets/prepare_public_perturbseq.py prepare \
  --dataset-id adamson_2016_perturbseq \
  --h5ad carf_benchmark/raw/adamson_2016_perturbseq/Adamson.h5ad \
  --max-perturbations 3
```

5. The preparation step writes:
   - `method_scores.csv`
   - `ground_truth.csv`
   - `covariates.csv`
   - `run_manifest.json`
6. Validate all CSVs with `validate_schema.py`.
7. Add the run to `carf_benchmark/configs/benchmark_v1.json`.
8. Run `Rscript --vanilla code/21_run_carf_benchmark_v1.R`.

## Adding a New Model

1. Add the model to `carf_benchmark/registry/models.csv`.
2. Check `carf_benchmark/adapters/model_output_specs.csv` for accepted native
   output modes.
3. Convert native outputs with `carf_benchmark/adapters/run_adapter.py`.
4. Preserve native outputs in the run directory and point `source_file` to them.
5. Ensure ranks are computed within the model's actual evaluable scope.
6. Append adapter output into the run `inputs/method_scores.csv`.
7. Re-run schema validation and the leaderboard.

Examples:

```bash
python3 carf_benchmark/adapters/run_adapter.py gene-score-csv \
  --model-id scgpt \
  --dataset-id adamson_2016_perturbseq \
  --perturbation-id DDIT3 \
  --native-output scgpt_DDIT3_gene_scores.csv \
  --output carf_benchmark/runs/adamson_2016_perturbseq/adapters/scgpt_DDIT3.csv

python3 carf_benchmark/adapters/run_adapter.py expression-delta \
  --model-id gears \
  --dataset-id norman_2019_combo \
  --perturbation-id AHR_FEV \
  --predicted-expression gears_AHR_FEV_predicted_expression.csv \
  --baseline-expression control_expression.csv \
  --output carf_benchmark/runs/norman_2019_combo/adapters/gears_AHR_FEV.csv
```

## Minimal Adapter Contract

Every adapter must output:

- `dataset_id`
- `perturbation_id`
- `model_id`
- `model_name`
- `gene_symbol`
- `score`
- `score_abs`
- `rank`
- `scope_included`
- `score_type`
- `source_file`

Optional but strongly recommended:

- `gene_id`
- `direction`
- `model_family`
- `model_version`
- `adapter_version`

## CI and Docker

The GitHub Actions workflow `.github/workflows/carf-benchmark.yml` runs:

1. R unit tests
2. Python adapter smoke tests
3. public source-manifest materialization
4. seed input preparation
5. schema validation
6. leaderboard generation
7. artifact upload for run-level leaderboard, model summary, audit matrix, PSR
   curves, and dataset readiness

The repository `Dockerfile` provides a reproducible environment for running
the same commands locally or in CI.
