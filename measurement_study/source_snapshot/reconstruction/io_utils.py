from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def input_path(dataset):
    entries = json.loads((ROOT / "evidence/recovery_manifest.json").read_text())["external_inputs"]
    return Path(next(x["source"] for x in entries if Path(x["source"]).name == dataset + ".h5ad"))


def text_array(values):
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def column(group, name):
    x = group[name]
    if isinstance(x, h5py.Dataset):
        return text_array(x[:]) if x.dtype.kind in "OSU" else x[:]
    codes = x["codes"][:]
    cats = text_array(x["categories"][:])
    if np.any(codes < 0):
        raise ValueError(f"Missing values in categorical field {name}")
    return cats[codes]


def gene_symbols(f):
    for key in ["Gene_symbol", "gene_name", "_index"]:
        if key in f["var"]:
            return column(f["var"], key)
    raise ValueError("No supported gene symbol field")


def dense_rows(f, indices):
    indices = np.asarray(indices, int)
    if len(np.unique(indices)) != len(indices):
        raise ValueError("Cell split contains duplicate cells")
    # h5py requires sorted unique index lists. Moment calculations ignore order.
    return f["X"][np.sort(indices), :].astype(np.float64)


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
