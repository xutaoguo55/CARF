#!/usr/bin/env python3
"""
Run GEARS inference on Norman 2019 and convert to CARF method_scores.csv format.

GEARS outputs predicted expression delta for each gene after perturbation.
Higher absolute delta = stronger predicted perturbation effect.
"""

import sys, os, warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'gears_repo'))
from gears import PertData, GEARS

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GEARS_DIR = os.path.join(BASE_DIR, '..', 'gears_repo')
RUNS_DIR = os.path.join(BASE_DIR, '..', 'carf_benchmark', 'runs')
os.chdir(GEARS_DIR)

# Dataset config
DATASET_ID = 'norman_2019_combo'
PERTURBATIONS = ['KLF1', 'BAK1', 'CEBPE']

MODEL_ID = 'gears'
MODEL_NAME = 'GEARS'
MODEL_FAMILY = 'deep_perturbation_model'
MODEL_VERSION = 'pretrained_norman_v1'
SCORE_TYPE = 'predicted_delta'

print("=== GEARS → CARF: Norman 2019 Inference ===\n")

# ---- Load GEARS data and model ----
print("Loading PertData...")
pert_data = PertData(GEARS_DIR)
pert_data.load(data_path=os.path.join(GEARS_DIR, 'norman_umi_go'))
pert_data.prepare_split(split='no_test', seed=1)
pert_data.get_dataloader(batch_size=32, test_batch_size=128)

print("Loading pre-trained GEARS model...")
gears_model = GEARS(pert_data, device='cpu',
                    weight_bias_track=False,
                    proj_name='gears',
                    exp_name='gears_norman')
gears_model.load_pretrained(os.path.join(GEARS_DIR, 'model_ckpt', 'model_ckpt'))
print("Model loaded.\n")

# Get gene mapping: GEARS gene_names is a pandas Series with ENSEMBL ID index and gene symbol values
gene_names_list = list(pert_data.gene_names.values)  # list of gene symbols
n_genes = len(gene_names_list)
print(f"GEARS gene universe: {n_genes} genes")

# ---- Read existing method_scores to match gene format ----
existing_ms_file = os.path.join(RUNS_DIR, DATASET_ID, 'inputs', 'method_scores.csv')
existing_ms = pd.read_csv(existing_ms_file)

# Build gene universe from existing method_scores (for one perturbation)
existing_genes = existing_ms[existing_ms['perturbation_id'] == PERTURBATIONS[0]][['gene_symbol']].drop_duplicates()
existing_gene_set = set(existing_genes['gene_symbol'].values)
print(f"Existing CARF gene universe: {len(existing_gene_set)} genes")

# Check overlap
gears_gene_set = set(str(g) for g in gene_names_list)
overlap = gears_gene_set & existing_gene_set
print(f"Gene overlap: {len(overlap)} ({100*len(overlap)/len(existing_gene_set):.1f}% of CARF genes)")

# ---- Run inference for each perturbation ----
# GEARS predict() returns post-perturbation expression levels, not deltas.
# Compute delta = predicted_post - control_mean
ctrl_expr = gears_model.ctrl_expression.numpy()  # mean expression in control cells

all_rows = []

for pert in PERTURBATIONS:
    print(f"\nRunning GEARS inference for {pert}...")
    result = gears_model.predict([[pert]])
    pred_post = result[pert]  # post-perturbation expression, shape (n_genes,)

    # Compute perturbation delta
    pred_delta = pred_post - ctrl_expr

    # Build gene-level predictions: map GEARS gene name → predicted delta
    gene_preds = {}
    for i in range(n_genes):
        gene_name = str(gene_names_list[i])
        gene_preds[gene_name] = float(pred_delta[i])

    # Match with CARF gene universe
    matched_genes = []
    for gene_symbol in existing_gene_set:
        if gene_symbol in gene_preds:
            matched_genes.append({
                'gene_symbol': gene_symbol,
                'score': gene_preds[gene_symbol],
            })

    df_pred = pd.DataFrame(matched_genes)
    df_pred['score_abs'] = df_pred['score'].abs()
    df_pred = df_pred.sort_values('score_abs', ascending=False).reset_index(drop=True)
    df_pred['rank'] = range(1, len(df_pred) + 1)
    df_pred['direction'] = df_pred['score'].apply(lambda x: 'up' if x > 0 else ('down' if x < 0 else 'zero'))

    # Add CARF metadata
    df_pred['dataset_id'] = DATASET_ID
    df_pred['perturbation_id'] = pert
    df_pred['model_id'] = MODEL_ID
    df_pred['model_name'] = MODEL_NAME
    df_pred['model_family'] = MODEL_FAMILY
    df_pred['model_version'] = MODEL_VERSION
    df_pred['gene_id'] = ''
    df_pred['scope_included'] = True
    df_pred['score_type'] = SCORE_TYPE
    df_pred['source_file'] = f'gears_inference_{pert}.npy'
    df_pred['adapter_version'] = 'gears_to_carf_v1'

    # Match column order of existing method_scores
    col_order = ['dataset_id', 'perturbation_id', 'model_id', 'model_name',
                 'model_family', 'model_version', 'gene_id', 'gene_symbol',
                 'score', 'score_abs', 'rank', 'direction', 'scope_included',
                 'score_type', 'source_file', 'adapter_version']
    df_pred = df_pred[col_order]

    n_matched = len(df_pred)
    n_lost = len(existing_gene_set) - n_matched
    print(f"  Matched: {n_matched} genes, Lost: {n_lost} genes")
    print(f"  Score range: [{df_pred['score'].min():.4f}, {df_pred['score'].max():.4f}]")
    print(f"  Top 5 by abs score: {df_pred['gene_symbol'].head(5).tolist()}")

    all_rows.append(df_pred)

# ---- Merge with existing and save ----
print("\n--- Merging with existing method_scores ---")

# Read existing method_scores and filter out any previous GEARS entries
existing_ms = existing_ms[existing_ms['model_id'] != MODEL_ID]

# Append GEARS rows
gears_ms = pd.concat(all_rows, ignore_index=True)
merged_ms = pd.concat([existing_ms, gears_ms], ignore_index=True)

# Save
output_file = os.path.join(RUNS_DIR, DATASET_ID, 'inputs', 'method_scores.csv')
# Backup first
import shutil
backup_file = output_file.replace('.csv', '_backup_before_gears.csv')
shutil.copy(output_file, backup_file)
print(f"Backup saved: {backup_file}")

merged_ms.to_csv(output_file, index=False)
print(f"Saved: {output_file}")
print(f"Total rows: {len(merged_ms)}")
print(f"Models in file: {merged_ms['model_id'].unique().tolist()}")

# Summary
print("\n=== Summary ===")
for pert in PERTURBATIONS:
    sub = merged_ms[(merged_ms['perturbation_id'] == pert) & (merged_ms['model_id'] == MODEL_ID)]
    print(f"  {pert}: {len(sub)} genes, score range [{sub['score'].min():.3f}, {sub['score'].max():.3f}]")

print("\nDone!")
