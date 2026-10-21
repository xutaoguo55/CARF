# License scope and original data attribution

Original CARF project contributions in this package (analysis code, documentation and analysis outputs, to the extent the project holds rights) are licensed under Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0). The complete terms are in `LICENSE`: https://creativecommons.org/licenses/by-nc/4.0/.

This grant does not replace upstream rights. Original source data, extracted source subsets and third-party metadata retain their original terms; the project adds no noncommercial restriction to their independently licensed content. No external Python or R package binaries or environments are bundled.

## RPE1 Perturb-seq

Replogle JM et al. *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq*. Cell (2022). https://doi.org/10.1016/j.cell.2022.05.013.

The processed source dataset by Joseph Replogle and Jonathan Weissman is https://doi.org/10.25452/figshare.plus.20029387.v1. The deposited record identifies CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Its public metadata snapshot is `provenance/figshare_article_20029387.json`; raw file URLs and digests are in `provenance/external_raw_requirements.json`.

Source-derived subsets in `data/rpe1/preparation/`, including `Dc_counts.npz`, cell/gene/target identities and selection/split tables, retain CC BY 4.0 for the underlying source content. They have been selected, reindexed and partitioned for this analysis. Normalized means, response contrasts and model outputs are derived as documented in the frozen protocol and source code. These transformations do not imply endorsement by the data creators. The full upstream raw matrices are not redistributed here.

## K562 development inputs

The curated Adamson/Norman H5AD input deposit is **Wei, Zhiting, Perturbation generalization H5ad datasets**, https://doi.org/10.5281/zenodo.14638780, whose metadata identifies CC BY 4.0. Source and verification records are in `provenance/development/`. Original studies are associated with GEO GSE90546 (Adamson) and GSE133344 (Norman). The license stated for the curated deposit is not represented as a new license grant for all independently hosted GEO material. The package contains development summaries and provenance, not the full H5AD inputs. Historical curated preprocessing is not recreated by the portable workflow.

## External software and published model architecture

Python dependencies (`numpy`, `pandas`, `scipy`, `h5py`, `matplotlib`) and R dependencies (`scTenifoldKnk`, `scTenifoldNet`, `Matrix`, `jsonlite`, `digest` and their dependencies) are obtained separately and remain under their own licenses. The project license does not apply to these dependencies. The original observed runtime and package versions are retained in `provenance/`.

The target-held-out bilinear baseline is an architecture-based adaptation of the published linear perturbation-prediction baseline, not a reproduction of its original benchmark. The design audit, original paper and official repository references are in `provenance/linear_baseline_design_audit_v1.json`. No original study result is claimed to come from a second validated deep prediction model.

The CC BY 4.0 legal text is supplied as `LICENSE-THIRD-PARTY-CC-BY-4.0.txt`. Public source metadata and license text retain their attribution; the package contains no claim of ownership over the original datasets or external software.
