#!/usr/bin/env python3
"""Patch Geneformer hardcoded CUDA references to use CPU/MPS."""
import re, os, sys

GENEFORMER_DIR = sys.argv[1] if len(sys.argv) > 1 else "/tmp/Geneformer"

files_to_patch = [
    "geneformer/emb_extractor.py",
    "geneformer/in_silico_perturber.py",
    "geneformer/perturber_utils.py",
]

for fp in files_to_patch:
    path = os.path.join(GENEFORMER_DIR, fp)
    if not os.path.exists(path):
        print(f"SKIP: {path} not found")
        continue
    with open(path) as f:
        content = f.read()
    original = content
    
    if "def _get_device" not in content:
        content = content.replace(
            "import torch",
            'import torch\n\ndef _get_device():\n    if torch.cuda.is_available():\n        return "cuda"\n    else:\n        return "cpu"\n'
        )
    
    content = content.replace('device="cuda"', 'device=_get_device()')
    content = content.replace('.to("cuda")', '.to(_get_device())')
    content = re.sub(
        r'torch\.cuda\.empty_cache\(\)',
        'torch.cuda.empty_cache() if torch.cuda.is_available() else None',
        content
    )
    
    if content != original:
        with open(path, "w") as f:
            f.write(content)
        print(f"PATCHED: {fp}")
    else:
        print(f"OK: {fp}")

