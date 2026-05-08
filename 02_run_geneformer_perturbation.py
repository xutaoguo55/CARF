#!/usr/bin/env python3
"""
Google Colab Notebook — Geneformer Zero-Shot WWOX Deletion Perturbation
========================================================================

PIPELINE STATUS (2026-04-29):
✅ Tokenization validated locally (420 cells, 3,992 tokens each, WWOX captured)
✅ MPS (Apple Silicon) patches applied — model runs on non-CUDA GPUs
✅ 3-cell validation run completed successfully with emb_mode="cls_and_gene"
   - 3 cells in 37.3 min on MPS (MacBook Pro)
   - Output: gene_embs_dict with 3,989 genes × 3 cosine values each
   - Cosine similarities ~0.9999 (deleting 1/3,992 genes = tiny effect)
⏳ Full 420-cell run pending (on Colab T4 GPU, est. 1-2 hours)

CRITICAL FINDINGS:
1. Geneformer's default model_input_size=4096 truncates WWOX (rank ~7,000
   in bulk data) from 99.5% of samples. Solution: pre-filter to top 4,000
   expressed genes + 22 candidate genes. WWOX token #15444 verified in all
   420 cells after pre-filtering.

2. Geneformer requires CUDA GPU by default. We patched 3 source files:
   emb_extractor.py, perturber_utils.py, in_silico_perturber.py to support
   MPS (Apple Silicon) and auto-detect the best available device.

3. Python 3.13 multiprocessing uses "spawn" by default, which fails with
   datasets.map(). Solution: set multiprocess start_method to "fork" before
   importing geneformer.

4. MPS performance: ~12-14 min/cell with 4,096-token sequences. Full 420
   cells would take ~87 hours on MPS. Colab T4 GPU recommended.

Files to upload to Google Drive (MyDrive/geneformer_benchmark/):
  GSE10846_filtered_v2.loom (4.7 MB) — pre-filtered loom file
  Optional: tokenized .dataset as a zip (faster than re-tokenizing)

Simplified Colab workflow:
  1. Upload GSE10846_filtered_v2.loom to Google Drive
  2. Run Cells 1 (install), 3-4 (dictionary files, tokenize)
  3. Run Cell 5 (load model), Cell 6 (perturbation, cls_and_gene mode)
  4. Run Cell 7 (stats), Cells 8-11 (inspect, validate, save)

KEY PARAMETERS:
  - emb_mode="cls_and_gene" (NOT "cell" or "cls") for gene-level ranking
  - input_data_file (NOT input_data_directory) for perturb_data()
  - forward_batch_size=32 for T4 GPU, =1 for MPS

# ===========================================================================
# CELL 1: Install Dependencies (~5 min)
# ===========================================================================
"""
!pip install -q scanpy anndata loompy mygene datasets

# Git LFS (required for Geneformer model weights)
!apt-get update -qq && apt-get install -y -qq git-lfs
!git lfs install

# Clone Geneformer from HuggingFace
!if [ ! -d "Geneformer" ]; then
    git clone https://huggingface.co/ctheodoris/Geneformer
fi

%cd Geneformer
!pip install -q .

import os, sys, glob, pickle
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
print("Dependencies installed")
print(f"Geneformer package at: {os.path.dirname(__file__) if '__file__' in dir() else os.getcwd()}")
"""

# ===========================================================================
# CELL 2: Mount Google Drive & Setup
# ===========================================================================
"""
from google.colab import drive
drive.mount('/content/drive')

# === CONFIGURE THESE PATHS ===
DRIVE_DATA_DIR = "/content/drive/MyDrive/geneformer_benchmark"

# Create local working directories (faster I/O)
!mkdir -p /content/data/raw_data
!mkdir -p /content/data/tokenized
!mkdir -p /content/data/perturbation_output

# Copy input files from Drive to local Colab VM
# Upload these TWO files to your Drive BEFORE running:
#   1. GSE10846_pseudo_counts_transposed.csv
#   2. gene_symbol_to_ensembl.csv
!cp "$DRIVE_DATA_DIR"/GSE10846_pseudo_counts_transposed.csv /content/data/
!cp "$DRIVE_DATA_DIR"/gene_symbol_to_ensembl.csv /content/data/

print("Data loaded")
!ls -lh /content/data/
"""

# ===========================================================================
# CELL 3: Map Gene Symbols → ENSEMBL & Create .loom File
# ===========================================================================
"""
import scanpy as sc
import loompy

# Load pseudo-count matrix (samples × genes, gene symbols as columns)
counts_df = pd.read_csv("/content/data/GSE10846_pseudo_counts_transposed.csv", index_col=0)
print(f"Pseudo-count matrix: {counts_df.shape[0]} samples x {counts_df.shape[1]} genes")

# Load gene symbol → ENSEMBL mapping (from R org.Hs.eg.db)
symbol_to_ensembl = pd.read_csv("/content/data/gene_symbol_to_ensembl.csv")
symbol_to_ensembl = symbol_to_ensembl.drop_duplicates(subset="SYMBOL", keep="first")

# Build mapping dict (one ENSEMBL ID per symbol)
symbol_map = {}
for _, row in symbol_to_ensembl.iterrows():
    if row['SYMBOL'] not in symbol_map:
        symbol_map[row['SYMBOL']] = row['ENSEMBL']

# Also try to fill missing genes using mygene
missing = [s for s in counts_df.columns if s not in symbol_map]
if missing:
    print(f"Mapping {len(missing)} missing genes via mygene...")
    import mygene
    mg = mygene.MyGeneInfo()
    results = mg.querymany(missing, scopes='symbol', fields='ensembl.gene', species='human')
    for r in results:
        if 'ensembl' in r and 'symbol' in r:
            ensembl_info = r['ensembl']
            if isinstance(ensembl_info, dict):
                eid = ensembl_info.get('gene')
            elif isinstance(ensembl_info, list):
                eid = ensembl_info[0].get('gene') if ensembl_info else None
            else:
                eid = None
            if eid:
                symbol_map[r['symbol']] = eid

print(f"Total gene symbol → ENSEMBL mappings: {len(symbol_map)}")

# Map columns to ENSEMBL
valid_symbols = [s for s in counts_df.columns if s in symbol_map]
ensembl_ids = [symbol_map[s] for s in valid_symbols]
counts_mapped = counts_df[valid_symbols].copy()
counts_mapped.columns = ensembl_ids

# Remove duplicate ENSEMBL IDs (keep first occurrence)
counts_mapped = counts_mapped.loc[:, ~counts_mapped.columns.duplicated()]

print(f"Mapped expression matrix: {counts_mapped.shape[0]} samples x {counts_mapped.shape[1]} ENSEMBL genes")

# Create AnnData → .loom (Geneformer tokenizer reads .loom files)
adata = sc.AnnData(
    X=counts_mapped.values.astype(np.float32),
    obs=pd.DataFrame(index=counts_mapped.index),
    var=pd.DataFrame(index=counts_mapped.columns)
)

loom_path = "/content/data/raw_data/GSE10846_bulk_pseudo_counts.loom"
loompy.create(
    loom_path,
    adata.X.T,  # loom convention: genes × cells
    row_attrs={"Gene": adata.var_names.to_numpy()},
    col_attrs={"CellID": adata.obs_names.to_numpy()}
)
print(f".loom file created: {loom_path}")
print(f"  Shape (genes × samples): {adata.X.T.shape}")
"""

# ===========================================================================
# CELL 4: Find Geneformer Dictionary Files
# ===========================================================================
"""
# The dictionary files are in the geneformer/ package directory
import geneformer

GENERFORMER_PKG = os.path.dirname(geneformer.__path__[0])
print(f"Geneformer package dir: {GENERFORMER_PKG}")

# Find .pkl files
pkl_files = []
for root, dirs, files in os.walk(GENERFORMER_PKG):
    for f in files:
        if f.endswith('.pkl'):
            pkl_files.append(os.path.join(root, f))

print(f"\nFound {len(pkl_files)} .pkl files:")
for f in sorted(pkl_files):
    print(f"  {f}")

# Locate required files
token_dict_file = None
gene_median_file = None

for f in pkl_files:
    basename = os.path.basename(f)
    if 'token_dictionary' in basename:
        token_dict_file = f
    if 'gene_median_dictionary' in basename:
        gene_median_file = f

print(f"\ntoken_dictionary: {token_dict_file}")
print(f"gene_median_dictionary: {gene_median_file}")

# Check WWOX is in vocabulary
if token_dict_file:
    with open(token_dict_file, 'rb') as f:
        token_dict = pickle.load(f)
    wwox_ensembl = "ENSG00000186153"
    if wwox_ensembl in token_dict:
        print(f"\nWWOX ({wwox_ensembl}) in vocabulary: token #{token_dict[wwox_ensembl]}")
    else:
        print(f"\nWARNING: WWOX not in vocabulary! Checking alternatives...")
        wwox_matches = [(k, v) for k, v in token_dict.items()
                        if isinstance(k, str) and 'ENSG00000186153' in k]
        print(f"  Matches: {wwox_matches}")
    print(f"Vocabulary size: {len(token_dict)}")
"""

# ===========================================================================
# CELL 5: Tokenize Bulk Data with TranscriptomeTokenizer
# ===========================================================================
"""
from geneformer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer(
    custom_attr_name_dict={},  # No cell metadata for bulk data
    nproc=4,
    gene_median_file=gene_median_file,
    token_dictionary_file=token_dict_file
)

print("Tokenizing bulk pseudo-count data...")
print(f"  Input: /content/data/raw_data/")
print(f"  Output: /content/data/tokenized/")

tk.tokenize_data(
    data_directory="/content/data/raw_data",
    output_directory="/content/data/tokenized",
    output_prefix="GSE10846_bulk"
)

print("Tokenization complete")
!ls -lh /content/data/tokenized/
"""

# ===========================================================================
# CELL 6: Load Pretrained Geneformer V2 Model
# ===========================================================================
"""
# For zero-shot perturbation, we need the pretrained model
# The model weights are in the cloned HuggingFace repo root (/content/Geneformer/)
# Files needed: config.json, model.safetensors

MODEL_DIR = "/content/Geneformer"
print(f"Model directory contents:")
for f in ['config.json', 'model.safetensors', 'generation_config.json']:
    path = os.path.join(MODEL_DIR, f)
    if os.path.exists(path):
        size_mb = os.path.getsize(path) / (1024*1024)
        print(f"  {f} ({size_mb:.1f} MB)")
    else:
        print(f"  {f} — MISSING")

# Load model via Geneformer API
from geneformer import GeneformerModel

model = GeneformerModel(
    model_type="Pretrained",
    model_version="V2",
    num_classes=0,
    emb_layer=-1,
    forward_batch_size=64
)
print("Geneformer V2 Pretrained model loaded")
"""

# ===========================================================================
# CELL 7: Run Zero-Shot WWOX Deletion Perturbation
# ===========================================================================
"""
from geneformer import InSilicoPerturber

WWWOX_ENSEMBL = "ENSG00000186153"

# CRITICAL: Use cls_and_gene for gene-level ranking
# cell or cls mode only gives cell-level cosine shift (needs "gene" in emb_mode)
# Validated on MPS: 3 cells, 37.3 min, 3,989 genes with cosine shifts
isp = InSilicoPerturber(
    perturb_type="delete",          # Simulate WWOX silencing
    genes_to_perturb=[WWWOX_ENSEMBL],
    model_type="Pretrained",
    num_classes=0,
    emb_mode="cls_and_gene",       # REQUIRED for gene-level ranking
    max_ncells=None,               # Use all 420 samples
    emb_layer=-1,
    forward_batch_size=32,         # T4 GPU batch size
    nproc=4
)

print("Running zero-shot WWOX deletion perturbation...")
print(f"  Model type: Pretrained (zero-shot)")
print(f"  Perturbation: delete {WWWOX_ENSEMBL} (WWOX)")
print(f"  emb_mode: cls_and_gene")
print(f"  Estimated: 1-2 hours on Colab T4 GPU")

# IMPORTANT: input_data_file expects path to .dataset directory
isp.perturb_data(
    model_directory=MODEL_DIR,
    input_data_file="/content/data/tokenized/GSE10846_filtered.dataset",
    output_directory="/content/data/perturbation_output",
    output_prefix="wwox_deletion"
)

print("Perturbation complete")
!ls -lh /content/data/perturbation_output/
"""

# ===========================================================================
# CELL 8: Extract Perturbation Statistics
# ===========================================================================
"""
from geneformer import InSilicoPerturberStats

# Correct API: no null_dist_data or nproc parameters
ispstats = InSilicoPerturberStats(
    mode="aggregate_gene_shifts",
    genes_perturbed=["ENSG00000186153"],
    combos=0,
)

print("Computing perturbation statistics (aggregate_gene_shifts)...")
ispstats.get_stats(
    input_data_directory="/content/data/perturbation_output",
    null_dist_data_directory=None,
    output_directory="/content/data/perturbation_output",
    output_prefix="wwox_deletion"
)

print("Statistics complete")
!ls -lh /content/data/perturbation_output/
"""

# ===========================================================================
# CELL 9: Load & Inspect Results
# ===========================================================================
"""
import glob

result_files = glob.glob("/content/data/perturbation_output/*.csv")
print("Result files:")
for f in sorted(result_files):
    print(f"  {os.path.basename(f)}  ({os.path.getsize(f)/1024:.1f} KB)")

# Load and display each result file
for f in sorted(result_files):
    df = pd.read_csv(f)
    basename = os.path.basename(f)
    print(f"\n{'='*60}")
    print(f"File: {basename}")
    print(f"Shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")
    print(f"\nFirst 10 rows:")
    print(df.head(10))
"""

# ===========================================================================
# CELL 10: Map ENSEMBL → Gene Symbols & Build Ranked List
# ===========================================================================
"""
# Invert mapping: ENSEMBL ID → gene symbol
ensembl_to_symbol = {}
for sym, eid in symbol_map.items():
    if eid not in ensembl_to_symbol:
        ensembl_to_symbol[eid] = sym

print(f"ENSEMBL → symbol mappings: {len(ensembl_to_symbol)}")

# Load the perturbation stats output
# Geneformer output typically has columns like:
#   gene, cosine_shift, p_value, fdr, etc.
# We need to identify the right file and columns

result_df = None
shift_col = None
pval_col = None

for f in sorted(result_files):
    df = pd.read_csv(f)
    cols = df.columns.tolist()

    # Identify key columns
    gene_col_candidates = ['gene', 'ensembl_id', 'Gene', 'gene_name', 'feature']
    shift_col_candidates = ['cosine_shift', 'cosine_sim_shift', 'shift',
                            'delta_cosine', 'effect_size', 'mean_shift']
    pval_col_candidates = ['p_value', 'pvalue', 'p_val', 'pval',
                           'p_adj', 'fdr', 'padj', 'q_value']

    gene_col = next((c for c in gene_col_candidates if c in cols), None)
    shift_col = next((c for c in shift_col_candidates if c in cols), None)
    pval_col = next((c for c in pval_col_candidates if c in cols), None)

    if gene_col and shift_col:
        result_df = df.copy()
        result_df['gene_symbol'] = result_df[gene_col].map(ensembl_to_symbol)
        print(f"Using file: {os.path.basename(f)}")
        print(f"  Gene column: {gene_col}")
        print(f"  Shift column: {shift_col}")
        print(f"  P-value column: {pval_col}")
        break

if result_df is None:
    print("ERROR: Could not identify result columns. Dumping all files:")
    for f in sorted(result_files):
        df = pd.read_csv(f)
        print(f"\n{os.path.basename(f)} columns: {df.columns.tolist()}")
        print(df.head(3))
else:
    # Sort by absolute shift (genes most affected by WWOX deletion)
    result_df['abs_shift'] = result_df[shift_col].abs()
    result_df = result_df.sort_values('abs_shift', ascending=False)
    result_df['rank'] = range(1, len(result_df) + 1)

    n_total = len(result_df)
    n_mapped = result_df['gene_symbol'].notna().sum()
    print(f"\nTotal genes in output: {n_total}")
    print(f"Genes with symbol mapping: {n_mapped}")
    print(f"\nTop 30 genes by |{shift_col}|:")
    print(result_df[['rank', 'gene_symbol', gene_col, shift_col]].head(30))
"""

# ===========================================================================
# CELL 11: Benchmark Against 17 Validated WWOX Targets
# ===========================================================================
"""
VALIDATED = [
    "PTGS2","MMP1","CXCL6","GJB2","NFKB1","RELA","RELB","NFKB2",
    "TMEM176A","TMEM176B","CD40","PCDHB5","PCDHB7","PCDHB10",
    "PCDHB14","PCDHB16","RAB34"
]

NON_VALIDATED = [
    "PCDHB4","PCDHB8","PCDHB11","PCDHB13","PCDHB15",
    "DSC2","DSG2","DSC3","GJB5","GJB6","GJB3","GJA4"
]

print("=" * 70)
print("GENEFORMER ZERO-SHOT WWOX DELETION — GROUND TRUTH VALIDATION")
print("=" * 70)

validated_ranks = []
validated_shifts = []

print(f"\n{'Gene':<15} {'Rank':>8} {'Percentile':>10} {shift_col:>12}")
print("-" * 50)

for gene in VALIDATED:
    rows = result_df[result_df['gene_symbol'] == gene]
    if len(rows) > 0:
        rank = rows['rank'].values[0]
        shift = rows[shift_col].values[0]
        pct = rank / n_total * 100
        validated_ranks.append(rank)
        validated_shifts.append(shift)
        print(f"{gene:<15} {rank:>8} {pct:>9.1f}% {shift:>12.6f}")
    else:
        print(f"{gene:<15} {'NOT FOUND':>8}")

print(f"\nNon-validated genes:")
for gene in NON_VALIDATED:
    rows = result_df[result_df['gene_symbol'] == gene]
    if len(rows) > 0:
        rank = rows['rank'].values[0]
        print(f"  {gene:<15} rank={rank}/{n_total}")

print(f"\n{'='*50}")
print(f"SUMMARY:")
print(f"  Validated genes found: {len(validated_ranks)}/17")
if validated_ranks:
    print(f"  Mean rank: {np.mean(validated_ranks):.1f} / {n_total}")
    print(f"  Median rank: {np.median(validated_ranks):.1f}")
    print(f"  Top-75 hits: {sum(1 for r in validated_ranks if r <= 75)}/17")
    print(f"  Top-100 hits: {sum(1 for r in validated_ranks if r <= 100)}/17")
    print(f"  Top-500 hits: {sum(1 for r in validated_ranks if r <= 500)}/17")
"""

# ===========================================================================
# CELL 12: Save Results for Download
# ===========================================================================
"""
# Save the complete Geneformer ranking
output_path = "/content/data/geneformer_wwox_deletion_ranking.csv"
result_df.to_csv(output_path, index=False)
print(f"Full ranking saved: {output_path}")
print(f"  {len(result_df)} genes ranked by WWOX deletion effect")

# Save benchmark summary
summary = {
    'method': 'Geneformer_V2_zero_shot',
    'perturbation': 'WWOX_deletion',
    'n_samples': counts_df.shape[0],
    'n_genes_tokenized': len(ensembl_ids),
    'n_genes_output': n_total,
    'n_validated_found': len(validated_ranks),
    'total_validated': 17,
}
if validated_ranks:
    summary['mean_rank'] = np.mean(validated_ranks)
    summary['median_rank'] = np.median(validated_ranks)
    summary['top75_hits'] = sum(1 for r in validated_ranks if r <= 75)
    summary['top100_hits'] = sum(1 for r in validated_ranks if r <= 100)

summary_df = pd.DataFrame([summary])
summary_df.to_csv("/content/data/geneformer_benchmark_summary.csv", index=False)

# Copy results back to Google Drive
!cp /content/data/geneformer_wwox_deletion_ranking.csv "$DRIVE_DATA_DIR"/
!cp /content/data/geneformer_benchmark_summary.csv "$DRIVE_DATA_DIR"/

print("Results copied to Google Drive:")
print(f"  {DRIVE_DATA_DIR}/geneformer_wwox_deletion_ranking.csv")
print(f"  {DRIVE_DATA_DIR}/geneformer_benchmark_summary.csv")
print("\nDownload these files for local analysis with evaluate_benchmark.R")
"""

print("\n" + "=" * 70)
print("GENEFORMER NOTEBOOK COMPLETE")
print("=" * 70)
print("""
Next Steps:
  1. Copy the two result CSV files from Google Drive to your local machine
  2. Place them in: foundation_model_benchmark/benchmark_results/
  3. Run: Rscript foundation_model_benchmark/04_evaluate_benchmark.R
     → This will generate unified metrics, figures, and the benchmark table
""")
