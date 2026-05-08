#!/usr/bin/env python3
"""Fix drug perturbation covariates: compute co-expression with known drug targets."""
import pandas as pd
import numpy as np
import os as _os

_PROJ_DIR = _os.environ.get("CARF_PROJ_DIR",
    _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
RUNS_DIR = _os.path.join(_PROJ_DIR, "carf_benchmark", "runs", "drug_perturbation_hsts")

# Known drug targets (Entrez IDs)
DRUG_TARGETS = {
    "Vorinostat": [3065, 3066, 8841, 10013],     # HDAC1, HDAC2, HDAC3, HDAC6
    "Imatinib": [613, 25, 3815, 5156, 5159],      # BCR, ABL1, KIT, PDGFRA, PDGFRB
    "Topotecan": [7150, 116447],                    # TOP1, TOP1MT
}

# Load expression matrix
expr = pd.read_csv("/tmp/drug_perturbation/expr_matrix.csv", index_col=0)
conds = pd.read_csv("/tmp/drug_perturbation/sample_conditions.csv", index_col=0)

# DMSO control samples
dmso_samples = conds[conds["drug"] == "DMSO"].index.tolist()
dmso_expr = expr[dmso_samples]

print(f"Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")
print(f"DMSO controls: {len(dmso_samples)} samples")

# Load existing covariates
cov = pd.read_csv(os.path.join(RUNS_DIR, "inputs", "covariates.csv"))
print(f"Covariates: {len(cov)} rows")

# For each drug, compute co-expression with targets
for drug, target_entrez in DRUG_TARGETS.items():
    # Find which targets are in expression matrix (index is int64)
    found_targets = [t for t in target_entrez if t in expr.index]
    print(f"\n{drug}: {len(found_targets)}/{len(target_entrez)} targets found in expr matrix")
    for t in found_targets:
        print(f"  Entrez {t}: present")

    if not found_targets:
        print(f"  WARNING: No targets found for {drug}, skipping")
        continue

    # Get target expression in DMSO controls
    target_expr = dmso_expr.loc[found_targets]

    # Compute Pearson correlation of each gene with each target in DMSO
    all_cors = []
    for gene_idx in expr.index:
        gene_expr = dmso_expr.loc[gene_idx]
        gene_cors = []
        for target_idx in target_expr.index:
            if gene_expr.std() > 0 and target_expr.loc[target_idx].std() > 0:
                cor = np.corrcoef(gene_expr, target_expr.loc[target_idx])[0, 1]
                if not np.isnan(cor):
                    gene_cors.append(cor)
        if gene_cors:
            all_cors.append(max(abs(c) for c in gene_cors))
        else:
            all_cors.append(0.0)

    coexpr_map = dict(zip(expr.index, all_cors))

    # Update covariates for this perturbation
    pert_mask = cov["perturbation_id"] == drug
    n_updated = 0
    for idx in cov[pert_mask].index:
        gene = int(cov.loc[idx, "gene_symbol"])
        if gene in coexpr_map:
            cov.loc[idx, "abs_coexpression_with_perturbed_gene"] = coexpr_map[gene]
            n_updated += 1

    print(f"  Updated {n_updated} genes with co-expression values")
    print(f"  Mean coexpr: {np.mean(all_cors):.4f}, Max: {np.max(all_cors):.4f}")

# Save updated covariates
cov.to_csv(os.path.join(RUNS_DIR, "inputs", "covariates.csv"), index=False)
print(f"\nSaved updated covariates to {RUNS_DIR}/inputs/covariates.csv")
