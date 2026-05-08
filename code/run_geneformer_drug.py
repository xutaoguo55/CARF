#!/usr/bin/env python3
"""Run Geneformer zero-shot perturbation on drug perturbation data (GSE96760).

Approach: Delete each drug's known protein target genes from the rank-value
encoding, measure cosine shift in embedding space, aggregate across targets.
"""
import os, sys, json, warnings, time, shutil, pickle
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')

# Critical: set fork BEFORE importing geneformer
from multiprocess import set_start_method
set_start_method("fork", force=True)

os.chdir("/tmp/Geneformer")

# ── Config ──────────────────────────────────────────────────────────
MODEL_DIR = "/tmp/Geneformer/Geneformer-V2-104M"  # Model checkpoint
GENEFORMER_DIR = "/tmp/Geneformer"  # Source code + token dicts
MODEL_VERSION = "V2"  # Uses 104M model, 4096 token input size
OUT_DIR = "/tmp/geneformer_drug_output"
TOKEN_DIR = f"{OUT_DIR}/tokenized"
LOOM_DIR = f"{OUT_DIR}/loom"
PERTURB_DIR = f"{OUT_DIR}/perturbation"
for d in [OUT_DIR, TOKEN_DIR, LOOM_DIR, PERTURB_DIR]:
    os.makedirs(d, exist_ok=True)

# Drug target Ensembl IDs (mapped from Entrez via MyGene.info)
DRUG_TARGET_ENSEMBL = {
    "Vorinostat": {
        "HDAC1": "ENSG00000116478",
        "HDAC2": "ENSG00000196591",
        "HDAC3": "ENSG00000171720",
        "HDAC6": "ENSG00000094631",
    },
    "Imatinib": {
        "BCR": "ENSG00000186716",
        "ABL1": "ENSG00000097007",
        "KIT": "ENSG00000157404",
        "PDGFRA": "ENSG00000134853",
        "PDGFRB": "ENSG00000113721",
    },
    "Topotecan": {
        "TOP1": "ENSG00000198900",
        "TOP1MT": "ENSG00000184428",
    },
}

# Drug doses to use (matching CARF ground truth)
DRUG_DOSES = {"Vorinostat": 0.2, "Imatinib": 2.0, "Topotecan": 0.0138}

N_TOP_GENES = 4000  # Pre-filter to top expressed genes (Geneformer token limit ~4096)

# ── Step 1: Load data and mapping ───────────────────────────────────
print("=" * 60)
print("STEP 1: Loading data and creating .loom files")
print("=" * 60)

expr = pd.read_csv("/tmp/drug_perturbation/expr_matrix.csv", index_col=0)
conds = pd.read_csv("/tmp/drug_perturbation/sample_conditions.csv", index_col=0)

with open("/tmp/drug_perturbation/entrez_to_ensembl.json") as f:
    entrez_to_ensembl = json.load(f)
print(f"Loaded {len(entrez_to_ensembl)} Entrez->Ensembl mappings")

# Build reverse map: Ensembl -> Entrez
ensembl_to_entrez = {v["ensembl"]: int(k) for k, v in entrez_to_ensembl.items()}

# DMSO control samples at 6h
dmso_samples = conds[(conds["drug"] == "DMSO") & (conds["time"] == 6)].index.tolist()
dmso_expr = expr[dmso_samples]

def prepare_loom(drug_name, drug_samples, dmso_samples):
    """Create a .loom file with DMSO controls + drug-treated samples."""
    import loompy

    all_samples = dmso_samples + drug_samples
    sub_expr = expr[all_samples].copy()

    # Compute mean expression in DMSO for gene ranking
    dmso_mean = sub_expr[dmso_samples].mean(axis=1)

    # Get drug target Ensembl IDs
    target_ensembl_list = list(DRUG_TARGET_ENSEMBL[drug_name].values())

    # Build Ensembl ID column from mapping
    ensembl_ids = []
    entrez_kept = []
    for eid in sub_expr.index:
        eid_str = str(int(eid))
        if eid_str in entrez_to_ensembl:
            ensembl_ids.append(entrez_to_ensembl[eid_str]["ensembl"])
            entrez_kept.append(eid_str)
        else:
            ensembl_ids.append(f"MISSING_{eid_str}")
            entrez_kept.append(eid_str)

    sub_expr["ensembl_id"] = ensembl_ids
    sub_expr["entrez_id"] = entrez_kept
    sub_expr["dmso_mean"] = dmso_mean

    # Keep only genes with Ensembl mapping AND in token dict
    with open(f"{GENEFORMER_DIR}/geneformer/token_dictionary_gc104M.pkl", "rb") as f:
        token_dict = pickle.load(f)
    has_token = sub_expr["ensembl_id"].isin(token_dict.keys())
    sub_expr = sub_expr[has_token]
    print(f"  {drug_name}: {sub_expr.shape[0]} genes with valid Ensembl tokens")

    # Ensure drug targets are included
    target_rows = sub_expr[sub_expr["ensembl_id"].isin(target_ensembl_list)]
    print(f"  {drug_name}: {len(target_rows)} drug targets found in expression")

    # Select top N by DMSO mean expression + drug targets
    sub_expr = sub_expr.sort_values("dmso_mean", ascending=False)
    top_genes = sub_expr.head(N_TOP_GENES).index.tolist()
    # Ensure drug targets are included
    for trow in target_rows.index:
        if trow not in top_genes:
            top_genes.append(trow)

    sub_expr = sub_expr.loc[top_genes]
    print(f"  {drug_name}: {len(top_genes)} genes selected (top {N_TOP_GENES} + targets)")

    # Extract expression values (without metadata columns) — genes × cells for loompy
    expr_cols = [c for c in sub_expr.columns if c in all_samples]
    expr_data = sub_expr[expr_cols].values.astype(np.int32)  # genes × cells

    # Row (gene) attributes — length = n_genes
    row_attrs = {
        "ensembl_id": np.array(sub_expr["ensembl_id"].values),
        "entrez_id": np.array(sub_expr["entrez_id"].values),
    }

    # Column (cell/sample) attributes — length = n_cells
    n_counts = expr_data.sum(axis=0)  # total counts per cell (sum over genes)
    col_attrs = {
        "n_counts": n_counts,
        "sample_id": np.array(expr_cols),
        "condition": np.array(["DMSO" if s in dmso_samples else drug_name for s in expr_cols]),
    }

    loom_path = f"{LOOM_DIR}/{drug_name}_drug_perturbation.loom"
    loompy.create(loom_path, expr_data, row_attrs, col_attrs, file_attrs={"drug": drug_name})
    print(f"  {drug_name}: Created {loom_path} ({expr_data.shape[0]} genes × {expr_data.shape[1]} cells)")

    return loom_path

for drug in ["Vorinostat", "Imatinib", "Topotecan"]:
    dose = DRUG_DOSES[drug]
    drug_samples = conds[(conds["drug"] == drug) & (conds["dose"] == dose) & (conds["time"] == 6)].index.tolist()
    print(f"\n{drug}: {len(drug_samples)} treated samples, {len(dmso_samples)} DMSO controls")
    prepare_loom(drug, drug_samples, dmso_samples)

print("\nDone creating .loom files. Ready for tokenization.")

# ── Step 2: Tokenize .loom files ───────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Tokenizing .loom files")
print("=" * 60)

from geneformer import TranscriptomeTokenizer

# Custom attribute: pass through condition label
tk = TranscriptomeTokenizer(
    custom_attr_name_dict={"condition": "condition", "sample_id": "sample_id"},
    nproc=1,
    model_version=MODEL_VERSION,
)
tk.tokenize_data(
    data_directory=LOOM_DIR,
    output_directory=TOKEN_DIR,
    output_prefix="drug_perturbation",
    file_format="loom",
)

print("Tokenization complete. Tokenized files:")
for f in sorted(os.listdir(TOKEN_DIR)):
    fp = os.path.join(TOKEN_DIR, f)
    if os.path.isdir(fp):
        print(f"  {f}/")
    else:
        print(f"  {f} ({os.path.getsize(fp)/1024:.0f} KB)")

# ── Step 3: Split combined dataset by condition ─────────────────────
print("\n" + "=" * 60)
print("STEP 3: Splitting dataset by drug condition")
print("=" * 60)

from datasets import Dataset

combined_ds_path = os.path.join(TOKEN_DIR, "drug_perturbation.dataset")
ds = Dataset.load_from_disk(combined_ds_path)
print(f"Combined dataset: {len(ds)} samples")
print(f"Conditions: {dict(zip(*np.unique(ds['condition'], return_counts=True)))}")

# Create per-drug filtered datasets
for drug_name in DRUG_TARGET_ENSEMBL:
    drug_ds = ds.filter(lambda x: x["condition"] in ["DMSO", drug_name])
    drug_ds_path = os.path.join(TOKEN_DIR, f"{drug_name}_filtered.dataset")
    drug_ds.save_to_disk(drug_ds_path)
    print(f"  {drug_name}: {len(drug_ds)} samples saved to {drug_ds_path}")

# ── Step 4: Run perturbation inference ──────────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Running perturbation inference")
print("=" * 60)

from geneformer import InSilicoPerturber

all_cosines = {}  # drug_name -> {entrez_id: max_abs_cosine_shift}

for drug_name, target_dict in DRUG_TARGET_ENSEMBL.items():
    print(f"\n--- {drug_name} ---")
    target_ensembl = list(target_dict.values())
    target_names = list(target_dict.keys())

    # Use per-drug filtered dataset
    dataset_path = os.path.join(TOKEN_DIR, f"{drug_name}_filtered.dataset")
    if not os.path.exists(dataset_path):
        print(f"  WARNING: No filtered dataset found at {dataset_path}")
        continue
    print(f"  Dataset: {dataset_path}")
    print(f"  Targets: {', '.join(f'{n}={e}' for n, e in zip(target_names, target_ensembl))}")

    # Run perturbation for each target gene individually
    drug_cosine_shifts = {}  # target_name -> {ensembl_id: cosine_shift}

    for target_name, target_ensg in zip(target_names, target_ensembl):
        pert_out = os.path.join(PERTURB_DIR, f"{drug_name}_{target_name}")
        if os.path.exists(pert_out):
            shutil.rmtree(pert_out)
        os.makedirs(pert_out, exist_ok=True)

        isp = InSilicoPerturber(
            perturb_type="delete",
            genes_to_perturb=[target_ensg],
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cls_and_gene",
            max_ncells=100,
            emb_layer=-1,
            forward_batch_size=1,
            nproc=1,
            model_version=MODEL_VERSION,
        )

        start = time.time()
        try:
            isp.perturb_data(
                model_directory=MODEL_DIR,
                input_data_file=dataset_path,
                output_directory=pert_out,
                output_prefix=f"{drug_name}_{target_name}",
            )
            elapsed = time.time() - start
            print(f"  {target_name}: completed in {elapsed:.1f}s")

            # Read cosine shift results
            cos_file = os.path.join(pert_out, f"{drug_name}_{target_name}_cosine_shifts.parquet")
            if os.path.exists(cos_file):
                cos_df = pd.read_parquet(cos_file)
                print(f"    Cosine shifts: {len(cos_df)} genes")
                drug_cosine_shifts[target_name] = cos_df
            else:
                print(f"    WARNING: No cosine shift file found at {cos_file}")
                # List what was produced
                for f in sorted(os.listdir(pert_out)):
                    print(f"      {f}")
        except Exception as e:
            elapsed = time.time() - start
            print(f"  {target_name}: FAILED after {elapsed:.1f}s: {e}")
            import traceback
            traceback.print_exc()

    # Aggregate across targets: max absolute cosine shift per gene
    if drug_cosine_shifts:
        # Build gene-level aggregation
        gene_shifts = {}  # ensembl_id -> list of cosine_shifts
        for tname, cos_df in drug_cosine_shifts.items():
            for _, row in cos_df.iterrows():
                ensg = row.get("ensembl_id") or row.get("gene")
                if ensg is None:
                    continue
                shift_val = row.get("cosine_shift") or row.get("shift") or 0.0
                if ensg not in gene_shifts:
                    gene_shifts[ensg] = []
                gene_shifts[ensg].append(shift_val)

        # Max absolute cosine shift
        aggregated = {}
        for ensg, shifts in gene_shifts.items():
            # Find shift with max absolute value
            best = max(shifts, key=lambda x: abs(x))
            entrez = ensembl_to_entrez.get(ensg)
            if entrez is not None:
                aggregated[entrez] = best

        all_cosines[drug_name] = aggregated
        print(f"  Aggregated: {len(aggregated)} genes with cosine shifts")
        print(f"    Range: [{min(aggregated.values()):.6f}, {max(aggregated.values()):.6f}]")
        print(f"    Mean abs: {np.mean([abs(v) for v in aggregated.values()]):.6f}")

# ── Step 5: Integrate into CARF method_scores.csv ──────────────────
print("\n" + "=" * 60)
print("STEP 5: Integrating Geneformer results into CARF")
print("=" * 60)

# Auto-detect project root relative to this script
import os as _os
_PROJ_DIR = _os.environ.get("CARF_PROJ_DIR",
    _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
RUNS_DIR = _os.path.join(_PROJ_DIR, "carf_benchmark", "runs", "drug_perturbation_hsts")
ms = pd.read_csv(f"{RUNS_DIR}/inputs/method_scores.csv")

# Check existing model_ids to avoid conflict
existing_models = ms["model_id"].unique()
print(f"Existing models: {list(existing_models)}")

# Add Geneformer scores for each drug
new_rows = []
for drug_name, cosine_dict in all_cosines.items():
    pert_id = drug_name  # matches perturbation_id in CARF
    for entrez_id, cosine_shift in cosine_dict.items():
        new_rows.append({
            "perturbation_id": pert_id,
            "gene_symbol": entrez_id,
            "model_id": "geneformer_v2_104m",
            "model_name": "Geneformer V2-104M (zero-shot drug target deletion)",
            "model_family": "foundation_model",
            "score": cosine_shift,
            "score_abs": abs(cosine_shift),
            "direction": "up" if cosine_shift >= 0 else "down",
            "score_type": "cosine_shift",
            "rank": 0,  # will be re-ranked below
        })

gf_df = pd.DataFrame(new_rows)
if len(gf_df) > 0:
    print(f"Geneformer results: {len(gf_df)} rows")
    for pert in gf_df["perturbation_id"].unique():
        sub = gf_df[gf_df["perturbation_id"] == pert]
        print(f"  {pert}: {len(sub)} genes, abs range [{sub['score_abs'].min():.6f}, {sub['score_abs'].max():.6f}]")

    # Re-rank within each perturbation+model group
    for (pert, mod), grp in gf_df.groupby(["perturbation_id", "model_id"]):
        indices = grp.index
        new_ranks = gf_df.loc[indices, "score_abs"].rank(ascending=False, method="min").astype(int)
        gf_df.loc[indices, "rank"] = new_ranks

    # Remove old Geneformer entries if any
    ms = ms[ms["model_id"] != "geneformer_v2_104m"]

    # Append new
    ms = pd.concat([ms, gf_df], ignore_index=True)

    # Save
    ms.to_csv(f"{RUNS_DIR}/inputs/method_scores.csv", index=False)
    print(f"Updated method_scores.csv: {len(ms)} rows, {len(ms['model_id'].unique())} models")

    # Save raw cosine shifts for reference
    cosine_out = f"{RUNS_DIR}/geneformer_cosine_shifts.json"
    with open(cosine_out, "w") as f:
        json_out = {drug: {str(k): v for k, v in shifts.items()}
                    for drug, shifts in all_cosines.items()}
        json.dump(json_out, f, indent=2)
    print(f"Saved raw cosine shifts to {cosine_out}")
else:
    print("WARNING: No Geneformer results generated")

print("\nDONE. Ready to re-run EDR/CPS diagnostics.")

