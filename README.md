# CARF 2.0: measurement sensitivity in expression adjusted perturbation rankings

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23537089.svg)](https://doi.org/10.5281/zenodo.23537089)

This version accompanies **Separating model input and covariate measurement effects in expression adjusted perturbation rankings**. It is a bounded empirical evaluation of how measurement changes affect ranking agreement, response utility and target identity. It does not establish causal confounder removal, universal superiority of expression adjustment, a new prediction architecture or a validated practical decision threshold.

## Start here

- [Current study source and verification guide](measurement_study/README.md)
- [Version 2.0.0 code and derived-data archive](https://doi.org/10.5281/zenodo.23537089)
- [GitHub version 2.0.0](https://github.com/xutaoguo55/CARF/releases/tag/v2.0.0)

Download the full archive from Zenodo to run the portable verification workflow. The `measurement_study/` directory exposes the source, tests and verification documentation; the complete numeric payload and source receipts are in the archive. The archive README gives exact commands, dependencies and validation limits.

## Research scope

The experiment separates model-input and expression-covariate measurement changes while holding evaluation outcomes fixed. Two K562 studies provide development analyses. The RPE1 analysis uses 128 evaluation targets, 665 common readout genes, 35 model conditions and 512 ranking conditions with target-disjoint prediction folds. An explicitly post hoc comparison evaluates rank-OLS and five-expression-stratum backgrounds across all 256 existing measurement contexts and 768 ranking conditions.

High ranking agreement can coexist with changes in target-identity information. Adjustment effects depend on thinning level, fold, diagnostic and background; adverse and non-monotonic findings are retained. Shared response is not automatically false biology, and technical cell or GEM-group splits are not independent biological replicates.

## Reproducibility and source data

The archive retains source checksums, frozen scopes, predictions, reference responses, complete evaluation tables and figure inputs. Portable verification reconstructs results from those frozen predictions and references. It is not a fresh raw-to-model rerun, a new R model fit or independent biological validation. Large original count matrices are obtained from their public deposits.

Source datasets: [GSE90546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546), [GSE133344](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344), [Zenodo 14638780](https://doi.org/10.5281/zenodo.14638780), and [Figshare+ 20029387 v1](https://doi.org/10.25452/figshare.plus.20029387.v1). The study uses deposited processed count matrices and does not claim reproduction from FASTQ files.

## Historical version

Version 1 is preserved at [Zenodo](https://doi.org/10.5281/zenodo.20088473), [Git tag v1.0.0](https://github.com/xutaoguo55/CARF/tree/v1.0.0), and the [previous README](https://github.com/xutaoguo55/CARF/blob/196965458d1693598adb061d87a56b049326d5a1/README.md). Legacy code and results remain in their original directories. They belong to the earlier WWOX/DLBCL research direction and do not define or validate the claims of this rebuilt study. Use `measurement_study/` and the version 2 archive for the current analysis; version 1 workflow commands do not reproduce it.

## License and citation

Original project material is licensed under [CC BY-NC 4.0](LICENSE). Third-party source data and dependencies retain their own licenses and attribution requirements; see the archive notices. The repository license does not replace those terms.

> Guo, Xutao. CARF 2.0: Separating model input and covariate measurement effects in expression adjusted perturbation rankings (2.0.0). Zenodo. https://doi.org/10.5281/zenodo.23537089

This DOI identifies the software and derived-data release, not a published journal article.

Contact: Xutao Guo, gxt827@126.com.
