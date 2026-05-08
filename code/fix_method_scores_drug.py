#!/usr/bin/env python3
"""Fix control_coexpression_baseline method_scores to use drug target co-expression."""
import pandas as pd
import numpy as np
import os as _os

_PROJ_DIR = _os.environ.get("CARF_PROJ_DIR",
    _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
RUNS_DIR = _os.path.join(_PROJ_DIR, "carf_benchmark", "runs", "drug_perturbation_hsts")

ms = pd.read_csv(f"{RUNS_DIR}/inputs/method_scores.csv")
cov = pd.read_csv(f"{RUNS_DIR}/inputs/covariates.csv")

# Build lookup: (perturbation_id, gene_symbol) -> abs_coexpression
coexpr_lookup = {}
for _, row in cov.iterrows():
    key = (row["perturbation_id"], int(row["gene_symbol"]))
    coexpr_lookup[key] = row["abs_coexpression_with_perturbed_gene"]

# Update control_coexpression_baseline scores
mask = ms["model_id"] == "control_coexpression_baseline"
n_updated = 0
for idx in ms[mask].index:
    pert = ms.loc[idx, "perturbation_id"]
    gene = int(ms.loc[idx, "gene_symbol"])
    key = (pert, gene)
    if key in coexpr_lookup:
        val = coexpr_lookup[key]
        ms.loc[idx, "score"] = val
        ms.loc[idx, "score_abs"] = abs(val)
        ms.loc[idx, "direction"] = "up" if val >= 0 else "down"
        n_updated += 1

# Re-rank within each perturbation+model group
for (pert, mod), grp in ms.groupby(["perturbation_id", "model_id"]):
    indices = grp.index
    # Rank by score_abs descending
    new_ranks = ms.loc[indices, "score_abs"].rank(ascending=False, method="min").astype(int)
    ms.loc[indices, "rank"] = new_ranks

ms.loc[mask, "score_type"] = "correlation"

ms.to_csv(f"{RUNS_DIR}/inputs/method_scores.csv", index=False)
print(f"Updated {n_updated} co-expression baseline scores")
print(f"Saved updated method_scores.csv")

# Show summary
coexpr_ms = ms[ms["model_id"] == "control_coexpression_baseline"]
for pert in coexpr_ms["perturbation_id"].unique():
    sub = coexpr_ms[coexpr_ms["perturbation_id"] == pert]
    print(f"  {pert}: score_abs range [{sub['score_abs'].min():.4f}, {sub['score_abs'].max():.4f}], mean={sub['score_abs'].mean():.4f}")
