"""Equal-cell-budget calibration across disjoint technical library groups."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from carf_core import residualize, topk_weights, rank_metrics, welch_reference
from io_utils import HERE, input_path, column, dense_rows, gene_symbols, sha256, write_json
from run_pilot import rng_for, coexpression


def sample_matched_controls(pert, pool, gemgroup, rng):
    """Match library composition, rather than just total control cell count."""
    sampled = []
    for g in np.unique(gemgroup[pert]):
        n = int((gemgroup[pert] == g).sum())
        candidates = pool[gemgroup[pool] == g]
        if len(candidates) < n:
            raise ValueError(f"Insufficient matched controls in library {g}")
        sampled.extend(rng.choice(candidates, n, replace=False))
    return np.asarray(sampled, int)


def control_pools(controls, gemgroup, seed, dataset, mode, rotation):
    rng = rng_for(seed, dataset + mode + "_control_pools")
    if mode == "disjoint_library_pools":
        return [controls[(gemgroup[controls] - 1) % 4 == (i + rotation) % 4] for i in range(4)]
    groups = [[] for _ in range(4)]
    for g in np.unique(gemgroup[controls]):
        parts = np.array_split(rng.permutation(controls[gemgroup[controls] == g]), 4)
        for i in range(4):
            groups[i].extend(parts[i])
    return [np.asarray(x, int) for x in groups]


def responses(f, parts, n, cfg):
    pa, pb, pc, ca, cb, cc = parts
    all_indices = np.concatenate(parts)
    if len(np.unique(all_indices)) != len(all_indices):
        raise ValueError("Repeated cells between task arms")
    A = dense_rows(f, pa).mean(0) - dense_rows(f, ca).mean(0)
    B = dense_rows(f, pb).mean(0) - dense_rows(f, cb).mean(0)
    C, q, y = welch_reference(dense_rows(f, pc), dense_rows(f, cc),
                               cfg["reference_fdr"], cfg["reference_min_abs_delta"])
    return {"A": A, "B": B, "C": C, "q": q, "supported": y, "n": n}


def score(task, genes, D, cmean, cvar, wrong, cfg, seed):
    keep = genes != task["name"]
    A, B, C, y = [task[key][keep] for key in ["A", "B", "C", "supported"]]
    coexpr, _ = coexpression(D, genes, task["name"])
    coexpr = coexpr[keep]
    logscore = np.log10(np.maximum(np.abs(A), 1e-10))
    edr = residualize(logscore, cmean[keep])
    methods = {
        "repeatability_positive": np.abs(A),
        "legacy_EDR": edr,
        "legacy_two_stage_CPS": residualize(edr, coexpr),
        "joint_CPS": residualize(logscore, np.column_stack([cmean[keep], coexpr])),
        "control_mean_only": cmean[keep],
        "control_variance_only": cvar[keep],
        "control_target_coexpression_only": coexpr,
        "constant_zero_delta": np.zeros(len(A)),
        "random_ranking": rng_for(seed, task["name"] + "_rank").random(len(A)),
    }
    if wrong is not None:
        methods["wrong_perturbation"] = np.abs(wrong[keep])
    rows = []
    for method, values in methods.items():
        for k in [cfg["primary_k"]] + cfg["secondary_k"]:
            w = topk_weights(values, k)
            m = rank_metrics(values, y, k)
            m.update(task=task["name"], task_type=task["type"], method=method,
                     mean_BC_product=float(w @ (B * C) / k),
                     mean_abs_response_B=float(w @ np.abs(B) / k),
                     A_direction_aligned_response_B=float(w @ (np.sign(A) * B) / k))
            rows.append(m)
    return rows


def run(dataset, mode, seed, rotation, cfg, cfg_path):
    source = input_path(dataset)
    receipt_path = HERE / "provenance" / (dataset + "_file_verification.json")
    receipt = json.loads(receipt_path.read_text())
    if sha256(source) != receipt["local_sha256"]:
        raise ValueError("Curated source hash changed")
    out = HERE / "results/library_calibration_v2" / f"{dataset}_{mode}_seed{seed}"
    if (out / "manifest.json").exists():
        raise SystemExit(f"Completed run exists: {out}")
    out.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    code_files = [Path(__file__), HERE / "carf_core.py", HERE / "io_utils.py", HERE / "run_pilot.py"]
    code_hashes = {p.name: sha256(p) for p in code_files}
    metadata_path = HERE / "provenance" / (dataset + "_cells_with_original_guides.csv.gz")
    metadata = pd.read_csv(metadata_path).set_index("cell_barcode", verify_integrity=True)
    tasks, split_records = [], []
    n = cfg["cells_per_arm"]
    with h5py.File(source, "r") as f:
        genes = gene_symbols(f)
        if len(np.unique(genes)) != len(genes):
            raise ValueError("Duplicate gene identity")
        barcodes = column(f["obs"], f["obs"].attrs["_index"])
        labels = column(f["obs"], "perturbation")
        gemgroup = metadata.loc[barcodes, "gemgroup"].to_numpy(dtype=int)
        if not np.array_equal(metadata.loc[barcodes, "perturbation"].to_numpy(), labels):
            raise ValueError("Original metadata no longer agrees with source labels")
        controls = np.flatnonzero(labels == "control")
        pools = control_pools(controls, gemgroup, seed, dataset, mode, rotation)
        all_ctrl = np.concatenate(pools)
        assert len(np.unique(all_ctrl)) == len(all_ctrl)
        D = dense_rows(f, pools[3])
        cmean, cvar = D.mean(0), D.var(0, ddof=1)
        table = pd.crosstab(labels, (gemgroup - 1) % 4).reindex(columns=range(4), fill_value=0)
        table["single_gene"] = [(s != "control") and ("+" not in s) for s in table.index]
        table["included"] = table.single_gene & (table[list(range(4))].min(axis=1) >= n)
        table.to_csv(out / "task_inclusion_by_library.csv")
        for condition in table.index[table.included]:
            rng = rng_for(seed, dataset + mode + condition)
            pert = np.flatnonzero(labels == condition)
            if mode == "disjoint_library_pools":
                par = [rng.choice(pert[(gemgroup[pert] - 1) % 4 == (i+rotation) % 4], n, replace=False) for i in range(3)]
            else:
                selected = rng.choice(pert, 3*n, replace=False)
                par = list(np.array_split(selected, 3))
            car = [sample_matched_controls(par[i], pools[i], gemgroup, rng) for i in range(3)]
            parts = par + car
            if np.intersect1d(np.concatenate(parts), pools[3]).size:
                raise ValueError("Covariate control pool overlaps evaluation")
            task = responses(f, parts, n, cfg)
            task.update(name=condition, type="real_single_gene_perturbation")
            tasks.append(task)
            for arm, idx in zip(["pert_A", "pert_B", "pert_C", "ctrl_A", "ctrl_B", "ctrl_C"], parts):
                split_records.extend({"task": condition, "arm": arm, "cell_index": int(j), "gemgroup": int(gemgroup[j])} for j in idx)
        for b in range(cfg["control_null_partitions_per_run"]):
            name = f"control_vs_control_{b:03d}"
            rng = rng_for(seed, dataset + mode + name)
            par, car = [], []
            for i in range(3):
                p = rng.choice(pools[i], n, replace=False)
                remaining = pools[i][~np.isin(pools[i], p)]
                c = sample_matched_controls(p, remaining, gemgroup, rng)
                par.append(p)
                car.append(c)
            task = responses(f, par+car, n, cfg)
            task.update(name=name, type="control_split_null")
            tasks.append(task)
            for arm, idx in zip(["pert_A", "pert_B", "pert_C", "ctrl_A", "ctrl_B", "ctrl_C"], par+car):
                split_records.extend({"task": name, "arm": arm, "cell_index": int(j), "gemgroup": int(gemgroup[j])} for j in idx)
    split_records.extend({"task": "all_tasks", "arm": "covariates_D", "cell_index": int(j), "gemgroup": int(gemgroup[j])} for j in pools[3])
    pd.DataFrame(split_records).to_csv(out / "cell_splits.csv.gz", index=False)
    np.savez_compressed(out / "responses.npz", genes=genes, tasks=np.array([t["name"] for t in tasks]),
                        **{key: np.stack([t[key] for t in tasks]) for key in ["A", "B", "C", "q", "supported"]})
    pd.DataFrame({"gene_symbol": genes, "mean_D": cmean, "variance_D": cvar}).to_csv(out / "covariates.csv.gz", index=False)
    real = [t for t in tasks if t["type"] == "real_single_gene_perturbation"]
    wrong = {t["name"]: real[(i+1) % len(real)]["A"] for i, t in enumerate(real)}
    rows = []
    for task in tasks:
        rows.extend(score(task, genes, D, cmean, cvar, wrong.get(task["name"]), cfg, seed))
    metrics = pd.DataFrame(rows).assign(dataset=dataset, mode=mode, seed=seed, cells_per_arm=n)
    metrics.to_csv(out / "metrics.csv.gz", index=False)
    write_json(out / "manifest.json", {
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": dataset, "mode": mode, "seed": seed, "rotation": rotation,
        "real_tasks": len(real), "null_tasks": cfg["control_null_partitions_per_run"],
        "source_sha256": receipt["local_sha256"], "metadata_sha256": sha256(metadata_path),
        "protocol": cfg, "protocol_sha256": sha256(cfg_path), "code_sha256": code_hashes,
        "status": "development_technical_library_validation_not_biological_replicates",
        "outputs": {p.name: sha256(p) for p in out.iterdir() if p.is_file()},
    })
    print(f"Completed {dataset}, {mode}, seed {seed}: {len(real)} real, {cfg['control_null_partitions_per_run']} null tasks", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["Adamson", "Norman"], required=True)
    args = parser.parse_args()
    cfg_path = HERE / "config/library_calibration_v2.json"
    cfg = json.loads(cfg_path.read_text())
    for seed, rotation in zip(cfg["seeds"], cfg["library_rotations"]):
        for mode in cfg["modes"]:
            run(args.dataset, mode, seed, rotation, cfg, cfg_path)


if __name__ == "__main__":
    main()
