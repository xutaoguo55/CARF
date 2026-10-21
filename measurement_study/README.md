# CARF measurement study — version 2.0.0

Code and reproducibility data for the CARF measurement study of expression adjustment, ranking stability and perturbation-identity specificity. The study separates changes to model-input measurements from changes to the adjustment covariate, and evaluates raw ranking, expression-dependence residualization (EDR), and posthoc rank-OLS and expression-stratum matching diagnostics.

- Version archive: https://doi.org/10.5281/zenodo.23537089
- Code: https://github.com/xutaoguo55/CARF (new study under `measurement_study/`; earlier CARF files and versions remain available)
- Original project contributions: **CC BY-NC 4.0**; upstream data retain **CC BY 4.0**, as detailed in `THIRD_PARTY_NOTICES.md`.

The full archive contains the data needed for the portable workflow below. A source-only GitHub checkout does not contain all score exports and result tables: download and extract `CARF-measurement-study-2.0.0.zip` from the version archive and run the commands inside the extracted package. The manuscript is distributed separately.

## Portable score-to-table reproduction

Use Python 3.12 and an empty virtual environment. From the extracted package directory:

```sh
python3.12 -m venv ../carf-measurement-2.0.0-env
../carf-measurement-2.0.0-env/bin/python -m pip install -r requirements-review.txt
../carf-measurement-2.0.0-env/bin/python verify_release.py --output ../carf-measurement-2.0.0-check
```

Use a new output directory for each execution. The workflow verifies every included payload digest, validates safe Unicode metadata with `allow_pickle=False`, reconstructs ranking and utility calculations, runs the mathematical tests and renders the original and journal figure versions. It does not require the originating workstation or its Python environment.

The numerical scope is:

- 35 safe model-score exports; 128 targets, 665 evaluation genes and fixed three-fold membership.
- All original raw-ranking metrics, six independently calculated EDR conditions, all 14 nested Dc covariate vectors, and the 65,536-row original target table and its summaries.
- All 256 contexts × three rankings = **768 posthoc ranking/matching evaluations**, **98,304 target rows** and **491,520 stratum rows**, including strict summaries and paired/factorial contrasts. These are recalculated from frozen scores, responses and covariates; no model is refitted.
- The original mathematical test suite (66 tests), plus three tests for the public provenance adapter. Tests cover fractional ties, OLS residualization, target-held-out PCA/bilinear ridge, safe metadata, and matching-failure propagation.
- Both figure display versions: five PNG/PDF figures each, with **11 original and 10 journal source CSVs** compared byte-for-byte. Rendering can depend on fonts and PDF timestamps; source-table equality is the numeric check.

For integrity only, without third-party dependencies:

```sh
python3 verify_release.py --checksums-only
```

For the original figure version only:

```sh
../carf-measurement-2.0.0-env/bin/python source_snapshot/reconstruction/plot_rpe1_measurement.py \
  --evaluation data/rpe1/evaluation \
  --review-summary provenance/rpe1_independent_review/summary.json \
  --out ../carf-measurement-2.0.0-figures
```

Measured verification results and software versions are retained under `validation/`. The earlier isolated execution used a fresh virtual environment and a relocated package. Public-package verification separately checks the metadata adapter and the full portable numerical workflow. Neither procedure redownloads the 8.7 GB raw source or refits the native R model.

## Package contents and provenance

`source_snapshot/reconstruction/` contains the unchanged executed analysis scripts, import dependencies, configurations and original tests. `data/rpe1/` contains the 35 safe score exports, selected Dc count measurements, frozen evaluation references and all top-level evaluation tables. `data/posthoc_diagnostic_sensitivity_v1/` contains the complete supplementary result tables and independent audit. `figures/rpe1/` and `figures/journal/` contain figures and source CSVs. `data/development/` retains the Adamson/Norman development summaries used to describe the project history, including q=0.25 covariate-only/joint comparisons and held-out-linear baseline identity and overlap. These are development results, not independent replication.

`release_manifest.json` and `SHA256SUMS` describe the files actually included here. Original study manifests also list omitted large raw caches and the repetitive grids of 512 original evaluation NPZs and 256 posthoc context NPZs; these grids are **not** included or compared by the compact verifier. Their omissions and original hashes are explicit in `provenance/omitted_repeated_grids.json`. The original full-grid audits are retained separately.

Scientific data tables, numeric score files, figures, executed analysis sources and frozen configurations are byte-identical to the sealed study package. Public metadata views replace originating absolute paths with source-relative paths or placeholders. `provenance/public_artifact_mapping.json` records each changed metadata file's original digest, public digest and reason. In particular, two provenance manifests have path-only redactions. The non-code posthoc design note also omits internal workflow wording; its mathematical rules are unchanged and its freeze receipt retains the original design digest. The portable adapter first validates their public bytes and then resolves their declared original identities for frozen-manifest links. It does **not** claim that the public and original bytes are identical or recreate the omitted private metadata bytes. All scientific output hashes remain unchanged. Internal execution transcripts and prior packaging history are excluded; published verification summaries retain measured results.

## Interpretation and limitations

Native all-floor targets and the intentional constant-ranking control are preserved. Perfect ranking agreement for a constant score is not evidence of identifiable Top50 genes. Operational ranking p10 is calculated after averaging thinning seeds within each target; operational utility p90 is across target-by-seed absolute changes.

The rank-OLS and five-stratum matching rules were fixed after inspection of the original RPE1 results. They are **posthoc sensitivity analyses**, not independent confirmation. All original models, targets, genes, folds, response arms and null products are retained. Matching uses five expression-rank strata and complementary fractional weights; it does not enforce exact expression balance or establish exchangeability. All observed matching conditions were supported. Artificial capacity-deficit cases test propagation of missing values to the full fold. Small native-score differences remain flagged even when rank-OLS ranks them.

Technical GEM partitions are not biological replicates. Reproducible calculation does not establish model superiority, universal benefit from adjustment, statistical calibration or mechanistic validity. The bilinear baseline is an architecture-based adaptation; it does not reproduce the original published benchmark. Historical K562 curated feature/preprocessing limits remain applicable.

The original posthoc runner is retained unchanged, but its source-manifest identities bind the original frozen execution and require omitted evaluation grids. It cannot be pointed at a newly timestamped raw-to-model run without an explicitly new protocol/source record. `verify_posthoc_review.py` is the supported compact score-to-table route; do not weaken frozen hash gates or silently reuse the earlier freeze receipt for a different execution.

## External raw data and optional full RPE1 rerun

`provenance/external_raw_requirements.json` gives exact URLs, byte counts, MD5 and SHA256 values:

- RPE1 single cells: Figshare file 35775606, 8,700,873,216 bytes.
- RPE1 author pseudobulk companion: Figshare file 35775581, 95,350,546 bytes.
- Adamson/Norman development H5ADs: Zenodo record 14638780, respectively 2,275,857,650 and 3,419,847,558 bytes after decompression.

The following source-level sequence describes the executed RPE1 pipeline. **A full raw-to-model rerun and fresh R installation have not been tested from this distribution.** Use a separate empty directory, sufficient disk/RAM, and an independently provisioned R environment. The native runner requires scTenifoldKnk 2.0.0, Matrix, jsonlite and digest and uses two cores. `provenance/original_R_sessionInfo.txt` records the observed runtime; this archive does not provide or validate a portable R installation recipe.

Capture the Python interpreter before changing directories. Replace the illustrative `rerun_dir` with a new empty directory. `Rscript` and required R packages must already be available. The Python requirements do not install them.

```sh
review_python="$(cd ../carf-measurement-2.0.0-env/bin && pwd)/python"
rerun_dir="/absolute/path/to/new-empty-study"
mkdir -p "$rerun_dir"
cp -R source_snapshot/. "$rerun_dir/"
cd "$rerun_dir"
"$review_python" reconstruction/acquire_replogle_rpe1.py
curl -fL --retry 3 https://ndownloader.figshare.com/files/35775581 \
  -o reconstruction/raw/rpe1_raw_bulk_01.h5ad
"$review_python" reconstruction/verify_rpe1_contents.py
"$review_python" reconstruction/prepare_rpe1_measurement.py
"$review_python" reconstruction/run_rpe1_measurement_models.py
"$review_python" reconstruction/export_rpe1_safe_models.py
"$review_python" reconstruction/evaluate_rpe1_measurement.py
```

These commands have been inspected, not run as a fresh full pipeline for this distribution. Run from that study root so the R runner can resolve `reconstruction/.Rlib`. Scripts refuse to overwrite completed outputs. The safe-export stage repairs the original writer's object-typed string metadata, checks string identity against the frozen scope CSVs and preserves every numeric array and dtype. It performs no scientific transformation or refit.

Full K562 reruns additionally need the legacy external-input mapping and omitted fit/preparation archives documented in `provenance/path_and_dependency_audit.json`. No one-command historical K562 preprocessing or full-model rerun is claimed.
