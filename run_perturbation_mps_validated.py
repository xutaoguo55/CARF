#!/usr/bin/env python3
"""Run Geneformer zero-shot WWOX deletion perturbation with MPS support."""
import os, sys, warnings, time, glob, shutil
warnings.filterwarnings('ignore')

# CRITICAL: Set fork BEFORE importing geneformer to avoid spawn issues on Python 3.13
from multiprocess import set_start_method
set_start_method("fork", force=True)

os.chdir("/tmp/Geneformer")

from geneformer import InSilicoPerturber

data_dir = "/tmp/geneformer_test"
out_dir = f"{data_dir}/perturbation_output"
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.makedirs(out_dir, exist_ok=True)

MODEL_DIR = "/tmp/Geneformer"
WWWOX_ENSEMBL = "ENSG00000186153"

import torch
print(f"Device: {'cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu'}")

# Use cls_and_gene to get gene-level embedding shifts (needed for gene ranking)
isp = InSilicoPerturber(
    perturb_type="delete",
    genes_to_perturb=[WWWOX_ENSEMBL],
    model_type="Pretrained",
    num_classes=0,
    emb_mode="cls_and_gene",
    max_ncells=3,
    emb_layer=-1,
    forward_batch_size=1,
    nproc=1
)

print(f"Configuration:")
print(f"  Model: {MODEL_DIR}")
print(f"  Input: {data_dir}/tokenized_final/GSE10846_filtered.dataset")
print(f"  emb_mode: {isp.emb_mode}")
print(f"  ncells: 3")
print(f"  nproc: 1")

start = time.time()
print("\nRunning perturbation...")
try:
    isp.perturb_data(
        model_directory=MODEL_DIR,
        input_data_file=f"{data_dir}/tokenized_final/GSE10846_filtered.dataset",
        output_directory=out_dir,
        output_prefix="wwox_deletion"
    )
    elapsed = time.time() - start
    print(f"\nPerturbation complete in {elapsed:.1f}s")

    for f in sorted(glob.glob(f"{out_dir}/*")):
        print(f"  {os.path.basename(f)} ({os.path.getsize(f)/1024:.1f} KB)")
except Exception as e:
    elapsed = time.time() - start
    print(f"\nFailed after {elapsed:.1f}s: {e}")
    import traceback
    traceback.print_exc()
