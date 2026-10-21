"""Fit frozen RPE1 measurement models without loading evaluation responses.

All normalization denominators contain the complete 8,749 observed genes. The
supervised model trains separately within each held-out target fold. This file
exports predictions and numerical diagnostics; it never evaluates B/C outcomes.
"""
from __future__ import annotations
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess

import h5py
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, save_npz, load_npz

from evaluate_heldout_linear import predict_selected_fold, bilinear_solve
from evaluate_pearson_measurement import absolute_pearson
from io_utils import HERE, ROOT, sha256, write_json

CONFIG = HERE / "config/rpe1_measurement_validation_v1.json"
PREP = HERE / "results/rpe1_measurement_validation_v1/preparation"
OUT = HERE / "results/rpe1_measurement_validation_v1/models"
RUNNER = HERE / "run_measurement_model.R"


def checked_hash(path, expected):
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"Source hash mismatch: {path}")
    return actual


def rng_stream(protocol_id, description, seed):
    payload = f"{protocol_id}|{description}|{seed}".encode()
    return np.random.default_rng(int.from_bytes(hashlib.sha256(payload).digest()[:16], "little"))


def validate_counts(counts):
    x = np.asarray(counts)
    if x.ndim != 2 or not np.isfinite(x).all() or (x < 0).any() or not np.equal(x, np.floor(x)).all():
        raise ValueError("Expected finite nonnegative integer cell-by-gene counts")
    return x


def thin_full_counts(counts, q, rng):
    x = validate_counts(counts)
    if not 0 < q <= 1:
        raise ValueError("Invalid molecular-thinning probability")
    return x.astype(np.int64, copy=True) if q == 1 else rng.binomial(x.astype(np.int64), q)


def normalize_full_counts(counts):
    """Normalize before any model/candidate gene restriction."""
    x = validate_counts(counts)
    totals = x.sum(axis=1, dtype=np.float64)
    if np.any(totals <= 0):
        raise ValueError("Zero observed-gene library: retain failure, never replace cell")
    return np.log1p(x.astype(np.float64) * (10000.0 / totals[:, None])), totals


def frozen_target_folds(targets, protocol_id):
    if len(set(targets)) != len(targets):
        raise ValueError("Duplicate target identity")
    ordered = sorted(targets, key=lambda t: (hashlib.sha256(f"{protocol_id}|fold|{t}".encode()).hexdigest(), t))
    mapping = {t: i % 3 for i, t in enumerate(ordered)}
    return np.array([mapping[t] for t in targets], dtype=int)


def fit_selected_targets(post, all_folds, target_gene_indices, eval_indices, reference,
                         components=10, ridge=.1, permutation_seeds=(), protocol_id="rpe1_measurement_validation_v1"):
    """Predict selected held-out targets; no outcome or test-A argument is used in fit."""
    eval_indices = np.asarray(eval_indices, int)
    predictions = {"heldout_linear": np.full((len(eval_indices), post.shape[0]), np.nan),
                   "training_condition_mean": np.full((len(eval_indices), post.shape[0]), np.nan)}
    predictions.update({f"training_label_permutation_{s}": np.full_like(predictions["heldout_linear"], np.nan)
                        for s in permutation_seeds})
    fits = []
    for fold in np.unique(all_folds):
        train = np.flatnonzero(all_folds != fold)
        positions = np.flatnonzero(all_folds[eval_indices] == fold)
        test = eval_indices[positions]
        if not len(test):
            raise ValueError("Missing evaluation target in frozen fold")
        prediction, (G, b, P, W, singular_values) = predict_selected_fold(
            post, train, test, target_gene_indices, components, ridge)
        predictions["heldout_linear"][positions] = (prediction - reference[:, None]).T
        predictions["training_condition_mean"][positions] = b[:, 0] - reference
        orders = []; permutation_W = []
        for seed in permutation_seeds:
            order = rng_stream(protocol_id, f"linear_label_permutation|fold={fold}", seed).permutation(len(train))
            wp = bilinear_solve(post[:, train] - b, G, P[order], ridge)
            pp = G @ wp @ G[np.asarray(target_gene_indices)[test]].T + b - reference[:, None]
            predictions[f"training_label_permutation_{seed}"][positions] = pp.T
            orders.append(order); permutation_W.append(wp)
        fits.append(dict(fold=int(fold), training_indices=train, test_indices=test,
                         evaluation_positions=positions, G=G, b=b, P=P, W=W,
                         singular_values=singular_values,
                         permutation_orders=np.asarray(orders, int).reshape(len(permutation_seeds), len(train)),
                         permutation_W=np.asarray(permutation_W, float).reshape(len(permutation_seeds), components, components),
                         permutation_seeds=np.asarray(permutation_seeds, int)))
    if any(not np.isfinite(p).all() for p in predictions.values()):
        raise ValueError("Unfilled or nonfinite held-out prediction")
    return predictions, fits


def assert_native_configuration(cfg):
    expected = dict(package="scTenifoldKnk", version="2.0.0", nc_nNet=10, nc_nCells=400,
                    nc_nComp=3, nc_q=.9, ma_method="heat", heatT=10, dr_direction=False, cores=2)
    for name, value in expected.items():
        if cfg["native"][name] != value:
            raise ValueError(f"Existing native runner does not implement frozen {name}")


def native_fit(counts, genes, targets, q, seed, algorithm_seed, cfg, source_hashes):
    folder = OUT / f"q{q:g}_thin{seed}" / "native" / f"algorithm{algorithm_seed}"
    folder.mkdir(parents=True, exist_ok=True)
    mmwrite(folder / "control_counts.mtx", csr_matrix(counts.T))
    (folder / "genes.txt").write_text("\n".join(genes) + "\n")
    (folder / "targets.txt").write_text("\n".join(targets) + "\n")
    input_names = ["control_counts.mtx", "genes.txt", "targets.txt"]
    write_json(folder / "input_manifest.json", dict(dataset="RPE1", q=q, thinning_seed=seed,
        model_seed=algorithm_seed, protocol_sha256=sha256(CONFIG), runner_sha256=sha256(RUNNER),
        preparation_input_sha256=source_hashes, input_umi_mass=int(counts.sum()),
        input_shape=list(counts.T.shape), input_sha256={n: sha256(folder/n) for n in input_names}))
    print(f"Fitting native q={q} thinning_seed={seed} algorithm_seed={algorithm_seed}", flush=True)
    with (folder / "fit.log").open("w") as log:
        result = subprocess.run(["Rscript", str(RUNNER), str(folder)], cwd=ROOT,
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Native model failed; retained log: {folder/'fit.log'}")
    receipt = json.loads((folder/"model_manifest.json").read_text())
    checked_hash(folder/"input_manifest.json", receipt["input_manifest_sha256"])
    checked_hash(folder/"distances.csv", receipt["distances_sha256"])
    frame = pd.read_csv(folder/"distances.csv", index_col=0)
    if set(frame.index) != set(targets) or set(frame.columns) != set(genes):
        raise ValueError("Native output target/gene identity mismatch")
    scores = frame.loc[targets, genes].to_numpy(float)
    if not np.isfinite(scores).all() or (scores < 0).any():
        raise ValueError("Nonfinite or negative native score")
    return scores


def save_model(name, scores, data, q, seed, algorithm_seed, index, numerical,
               signed=None, target_zero_variance=None, candidate_zero_variance=None):
    scores = np.asarray(scores, float)
    if scores.shape != (len(data["evaluation_targets"]), len(data["candidate_genes"])):
        raise ValueError("Prediction does not match frozen evaluation scope")
    if not np.isfinite(scores).all() or (scores < 0).any():
        raise ValueError("Invalid raw score")
    folder = OUT / f"q{q:g}_thin{seed}" / "scores"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{name}_algorithm{algorithm_seed}.npz"
    arrays = dict(targets=data["evaluation_targets"], genes=data["candidate_genes"],
                  folds=data["evaluation_folds"], raw_scores=scores)
    if signed is not None:
        if signed.shape != scores.shape or not np.isfinite(signed).all():
            raise ValueError("Invalid signed prediction")
        arrays["signed_predictions"] = signed
    if target_zero_variance is not None:
        arrays["target_zero_variance"] = target_zero_variance
        arrays["candidate_zero_variance"] = candidate_zero_variance
    np.savez_compressed(path, **arrays)
    index.append(dict(model=name, q=q, thinning_seed=seed, algorithm_seed=algorithm_seed,
                      path=str(path.relative_to(OUT)), sha256=sha256(path),
                      n_targets=len(scores), n_candidates=scores.shape[1],
                      signed_predictions_available=signed is not None))
    for i, target in enumerate(data["evaluation_targets"]):
        numerical.append(dict(model=name, q=q, thinning_seed=seed, algorithm_seed=algorithm_seed,
            target=target, fold=int(data["evaluation_folds"][i]),
            score_floor_count=int((scores[i] <= 1e-10).sum()),
            all_scores_floored=bool(np.all(scores[i] <= 1e-10)),
            raw_exact_constant=bool(np.ptp(scores[i]) == 0),
            target_zero_variance=bool(target_zero_variance[i]) if target_zero_variance is not None else None,
            candidate_zero_variance_count=int(candidate_zero_variance.sum()) if candidate_zero_variance is not None else None,
            min_score=float(scores[i].min()), max_score=float(scores[i].max())))


def text_values(values):
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def load_preparation(cfg, protocol_hash):
    manifest = json.loads((PREP/"manifest.json").read_text())
    if manifest["protocol_sha256"] != protocol_hash:
        raise ValueError("Preparation used a different frozen protocol")
    if manifest["status"] != "completed_frozen_preparation_no_model_or_thinning" or manifest["source_sha256"] != cfg["source"]["sha256"]:
        raise ValueError("Preparation is incomplete or identifies a different raw source")
    required = ["genes_all.csv", "targets_all.csv", "evaluation_targets.csv", "model_genes.csv",
                "evaluation_genes.csv", "Dm_counts.npz", "Dm_cells.csv", "A_raw_counts.h5", "A_post_mean.npy",
                "metadata_manifest.json", "gene_freeze_manifest.json"]
    hashes = {}
    for name in required:
        expected = manifest["output_sha256"][name]
        hashes[name] = checked_hash(PREP/name, expected)
    genes = pd.read_csv(PREP/"genes_all.csv")
    targets = pd.read_csv(PREP/"targets_all.csv")
    evaluated = pd.read_csv(PREP/"evaluation_targets.csv")
    model = pd.read_csv(PREP/"model_genes.csv")
    candidates = pd.read_csv(PREP/"evaluation_genes.csv")
    dm_cells = pd.read_csv(PREP/"Dm_cells.csv")
    all_genes = genes.gene_id.astype(str).to_numpy()
    all_targets = targets.target_id.astype(str).to_numpy()
    all_folds = targets.fold.to_numpy(int)
    evaluation_targets = evaluated.target_id.astype(str).to_numpy()
    evaluation_indices = evaluated.all_target_index.to_numpy(int)
    model_indices = model.source_gene_index.to_numpy(int)
    candidate_indices = candidates.source_gene_index.to_numpy(int)
    if len(all_genes) != 8749 or len(all_targets) != 505:
        raise ValueError("Metadata-frozen source dimensions changed")
    if len(set(all_genes)) != len(all_genes) or len(set(evaluation_targets)) != len(evaluation_targets):
        raise ValueError("Duplicate source/evaluation identity")
    if not np.array_equal(genes.source_gene_index.to_numpy(int), np.arange(len(all_genes))):
        raise ValueError("Source gene order changed")
    if not np.array_equal(all_targets[evaluation_indices], evaluation_targets):
        raise ValueError("Evaluation target-to-source mapping changed")
    if not np.array_equal(all_folds, frozen_target_folds(all_targets, cfg["protocol_id"])):
        raise ValueError("Frozen target folds differ from prespecified SHA256 assignment")
    evaluation_folds = all_folds[evaluation_indices]
    if not np.array_equal(evaluated.fold.to_numpy(int), evaluation_folds):
        raise ValueError("Evaluation target fold mismatch")
    if not 64 <= len(evaluation_targets) <= 128 or np.min(np.bincount(evaluation_folds, minlength=3)) < 10:
        raise ValueError("Insufficient frozen evaluation targets/fold coverage")
    if len(set(model_indices)) != len(model_indices) or len(set(candidate_indices)) != len(candidate_indices):
        raise ValueError("Duplicate model/candidate gene")
    if len(candidate_indices) <= 50 or not set(candidate_indices) <= set(model_indices):
        raise ValueError("Invalid frozen candidate scope")
    gene_map = {g: i for i, g in enumerate(all_genes)}
    target_gene_indices = np.asarray([gene_map[t] for t in all_targets], int)
    model_map = {int(g): i for i, g in enumerate(model_indices)}
    evaluation_target_model_indices = np.asarray([model_map[gene_map[t]] for t in evaluation_targets], int)
    candidate_model_indices = np.asarray([model_map[int(i)] for i in candidate_indices], int)
    dm_counts = load_npz(PREP/"Dm_counts.npz").toarray()
    validate_counts(dm_counts)
    if dm_counts.shape != (400, 8749) or len(dm_cells) != 400:
        raise ValueError("Frozen Dm dimensions changed")
    with h5py.File(PREP/"A_raw_counts.h5") as f:
        if f.attrs["source_sha256"] != cfg["source"]["sha256"]:
            raise ValueError("A raw cache identifies a different source")
        if f["counts"].shape != (505, 40, 8749):
            raise ValueError("Frozen A raw count cache shape changed")
        if not np.array_equal(text_values(f["targets"][:]), all_targets) or not np.array_equal(text_values(f["genes"][:]), all_genes):
            raise ValueError("A cache gene/target identity changed")
        a_cells = f["cell_indices"][:]
        if a_cells.shape != (505, 40) or len(np.unique(a_cells)) != a_cells.size:
            raise ValueError("Duplicate/incomplete frozen A cell indices")
        if "cell_index" not in dm_cells or len(set(dm_cells.cell_index)) != 400:
            raise ValueError("Missing/duplicate Dm source cell indices")
        if np.intersect1d(a_cells, dm_cells.cell_index.to_numpy(int)).size:
            raise ValueError("Dm overlaps A training cells")
    for name in ["genes_all.csv", "targets_all.csv", "evaluation_targets.csv", "model_genes.csv", "evaluation_genes.csv", "Dm_cells.csv"]:
        (OUT/name).write_bytes((PREP/name).read_bytes())
    write_json(OUT/"scope_freeze.json", dict(created_before_any_fit_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=protocol_hash, preparation_manifest_sha256=sha256(PREP/"manifest.json"),
        required_preparation_inputs=hashes, full_gene_count=8749, eligible_training_targets=505,
        evaluation_target_count=len(evaluation_targets), model_gene_count=len(model_indices),
        evaluation_gene_count=len(candidate_indices), evaluation_targets_by_fold=np.bincount(evaluation_folds, minlength=3).tolist(),
        read_boundaries="No B/C, control-A, empirical-anchor, null, or adjustment-covariate values are loaded by model fitting."))
    return dict(all_genes=all_genes, all_targets=all_targets, all_folds=all_folds,
        evaluation_targets=evaluation_targets, evaluation_indices=evaluation_indices,
        evaluation_folds=evaluation_folds, model_indices=model_indices,
        model_genes=all_genes[model_indices], candidate_indices=candidate_indices,
        candidate_genes=all_genes[candidate_indices], target_gene_indices=target_gene_indices,
        evaluation_target_model_indices=evaluation_target_model_indices,
        candidate_model_indices=candidate_model_indices, dm_counts=dm_counts,
        preparation_hashes=hashes, preparation_manifest_sha256=sha256(PREP/"manifest.json"))


def condition_inputs(data, q, seed, cfg):
    folder = OUT/f"q{q:g}_thin{seed}"
    folder.mkdir(parents=True, exist_ok=True)
    dm_description = f"thin|Dm|q={q:g}"
    dm = thin_full_counts(data["dm_counts"], q, rng_stream(cfg["protocol_id"], dm_description, seed))
    dm_normalized, dm_totals = normalize_full_counts(dm)
    save_npz(folder/"Dm_thinned_counts.npz", csr_matrix(dm))
    np.savez_compressed(folder/"Dm_normalized_reference.npz", genes=data["all_genes"],
                        mean=dm_normalized.mean(axis=0), cell_library_sums=dm_totals)
    post = np.empty((8749, 505), float)
    audits = [dict(pool="Dm", target="", q=q, thinning_seed=seed, rng_description=dm_description,
        n_cells=400, source_UMI=int(data["dm_counts"].sum()), retained_UMI=int(dm.sum()),
        min_retained_cell_UMI=int(dm_totals.min()), zero_library_cells=0)]
    with h5py.File(PREP/"A_raw_counts.h5") as f:
        for i, target in enumerate(data["all_targets"]):
            raw = f["counts"][i]
            description = f"thin|A|target={target}|q={q:g}"
            thinned = thin_full_counts(raw, q, rng_stream(cfg["protocol_id"], description, seed))
            normalized, totals = normalize_full_counts(thinned)
            post[:, i] = normalized.mean(axis=0)
            audits.append(dict(pool="A", target=target, q=q, thinning_seed=seed, rng_description=description,
                n_cells=40, source_UMI=int(raw.sum()), retained_UMI=int(thinned.sum()),
                min_retained_cell_UMI=int(totals.min()), zero_library_cells=0))
    error = None
    if q == 1:
        baseline = np.load(PREP/"A_post_mean.npy", allow_pickle=False)
        if baseline.shape != (505, 8749):
            raise ValueError("Frozen A mean shape changed")
        error = float(np.max(np.abs(post.T-baseline)))
        if error > 1e-12:
            raise ValueError("q1 count-derived A means differ from frozen preparation")
    np.savez_compressed(folder/"A_training_post_means.npz", genes=data["all_genes"], targets=data["all_targets"],
                        folds=data["all_folds"], post=post)
    pd.DataFrame(audits).to_csv(folder/"measurement_input_audit.csv", index=False)
    write_json(folder/"measurement_input_manifest.json", dict(q=q, thinning_seed=seed,
        q1_A_post_reconstruction_max_absolute_error=error,
        normalization="log1p(counts/sum_all_8749_observed_genes*10000) per cell before gene selection",
        Dm_measurement_reuse="One whole-gene thinning draw shared by native, Pearson and linear Dm reference",
        A_measurement_reuse="One draw per target shared by every fold; each fit then selects only other-fold target columns",
        rng_rule=cfg["rng_rule"], counts_generated_once=True,
        zero_library_policy="Any zero library fails; no resampling or cell replacement"))
    return dm, dm_normalized, post


def run_condition(data, q, seed, cfg, index, numerical):
    print(f"Preparing model measurement q={q} seed={seed}: full-gene Dm and each of 505 A inputs", flush=True)
    dm, expression, post = condition_inputs(data, q, seed, cfg)
    model_expression = expression[:, data["model_indices"]].T
    pearson, zero_targets, zero_genes = absolute_pearson(model_expression, data["evaluation_target_model_indices"])
    if q == 1:
        with np.errstate(invalid="ignore", divide="ignore"):
            independent = np.nan_to_num(np.abs(np.corrcoef(model_expression))[data["evaluation_target_model_indices"]], nan=0.0)
        max_error = float(np.max(np.abs(independent-pearson)))
        if not np.allclose(independent, pearson, atol=1e-12, rtol=1e-12):
            raise ValueError("Pearson q1 independent formula disagrees")
        write_json(OUT/"q1_pearson_formula_check.json", dict(max_absolute_error=max_error, tolerance=1e-12))
    save_model("control_only_Pearson", pearson[:, data["candidate_model_indices"]], data, q, seed, 0, index, numerical,
        target_zero_variance=zero_targets, candidate_zero_variance=zero_genes[data["candidate_model_indices"]])
    permutations = cfg["training_label_permutation_seeds"] if q == 1 else []
    print(f"Fitting bilinear q={q} seed={seed}: three held-out folds", flush=True)
    predictions, fits = fit_selected_targets(post, data["all_folds"], data["target_gene_indices"],
        data["evaluation_indices"], expression.mean(axis=0), cfg["linear"]["components"],
        cfg["linear"]["ridge_penalty_each_side"], permutations, cfg["protocol_id"])
    fit_folder = OUT/f"q{q:g}_thin{seed}"/"linear_fits"
    fit_folder.mkdir(parents=True, exist_ok=True)
    permutation_rows = []
    for fit in fits:
        train = fit["training_indices"]; test = fit["test_indices"]
        np.savez_compressed(fit_folder/f"fold_{fit['fold']}.npz", **fit,
            genes=data["all_genes"], training_targets=data["all_targets"][train], test_targets=data["all_targets"][test])
        for p, permutation_seed in enumerate(permutations):
            permutation_rows.extend(dict(fold=fit["fold"], permutation_seed=permutation_seed,
                condition_target=data["all_targets"][train[j]],
                assigned_embedding_target=data["all_targets"][train[fit["permutation_orders"][p,j]]]) for j in range(len(train)))
        print(f"Completed bilinear q={q} seed={seed} fold={fit['fold']}, train={len(train)}, test={len(test)}", flush=True)
    if permutation_rows:
        pd.DataFrame(permutation_rows).to_csv(fit_folder/"training_label_permutations.csv", index=False)
    for name, prediction in predictions.items():
        signed = prediction[:, data["candidate_indices"]]
        save_model(name, np.abs(signed), data, q, seed, 0, index, numerical, signed=signed)
    algorithm_seeds = [cfg["native"]["baseline_seed"]]
    if q == 1:
        algorithm_seeds += cfg["native"]["additional_q1_algorithm_seeds"]
    for algorithm_seed in algorithm_seeds:
        native = native_fit(dm[:, data["model_indices"]], data["model_genes"], data["evaluation_targets"],
                            q, seed, algorithm_seed, cfg, data["preparation_hashes"])
        save_model("native_scTenifoldKnk", native[:, data["candidate_model_indices"]], data, q, seed,
                   algorithm_seed, index, numerical)
    print(f"Completed all model outputs q={q} seed={seed}", flush=True)


def main():
    cfg = json.loads(CONFIG.read_text()); protocol_hash = sha256(CONFIG)
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT/"manifest.json").exists():
        raise SystemExit("Completed RPE1 models exist; refusing overwrite")
    assert_native_configuration(cfg)
    code_files = [Path(__file__), RUNNER, HERE/"evaluate_heldout_linear.py", HERE/"evaluate_pearson_measurement.py",
                  HERE/"measurement_diagnostic_pilot.py", HERE/"carf_core.py", HERE/"io_utils.py"]
    code_hashes = {p.name: sha256(p) for p in code_files}
    start = datetime.now(timezone.utc).isoformat()
    try:
        data = load_preparation(cfg, protocol_hash)
        index = []; numerical = []
        conditions = [(1.0, 0)] + [(q, s) for q in cfg["thinning_probabilities"] for s in cfg["thinning_seeds"]]
        for q, seed in conditions:
            run_condition(data, q, seed, cfg, index, numerical)
        pd.DataFrame(index).to_csv(OUT/"model_index.csv", index=False)
        pd.DataFrame(numerical).to_csv(OUT/"numerical_information_audit.csv", index=False)
        if len(index) != 35:
            raise ValueError("Not all 35 prespecified model/measurement fits were exported")
        if sha256(CONFIG) != protocol_hash:
            raise ValueError("Frozen protocol changed during fitting")
        checked_hash(PREP/"manifest.json", data["preparation_manifest_sha256"])
        for name, expected in data["preparation_hashes"].items():
            checked_hash(PREP/name, expected)
        for path in code_files:
            checked_hash(path, code_hashes[path.name])
        outputs = {str(p.relative_to(OUT)): sha256(p) for p in sorted(OUT.rglob("*"))
                   if p.is_file() and p.name != "manifest.json"}
        write_json(OUT/"manifest.json", dict(status="completed_models_without_outcome_evaluation",
            started_utc=start, finished_utc=datetime.now(timezone.utc).isoformat(),
            protocol_sha256=protocol_hash, protocol=cfg, code_sha256=code_hashes,
            preparation_manifest_sha256=data["preparation_manifest_sha256"], preparation_inputs=data["preparation_hashes"],
            n_model_measurement_outputs=len(index), output_sha256=outputs,
            input_sharing="Full-gene Dm thinning shared across models; target-specific A thinning shared across folds",
            evaluation_boundary="No B/C, control-null or empirical-A outcomes, and no Dc adjustment covariates, loaded or evaluated",
            linear_adaptation="Published two-sided bilinear ridge, no ctrl training column, other-fold 505-target A training; current-condition Dm reference subtraction",
            feature_budgets="Native/Pearson use frozen Dm-selected model genes; linear uses all 8749 observed genes. Feature/supervision budgets differ.",
            environment=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__)))
        print("Completed all 35 prespecified RPE1 model outputs; no B/C evaluation performed", flush=True)
    except Exception as exc:
        write_json(OUT/"failure.json", dict(time_utc=datetime.now(timezone.utc).isoformat(), error=repr(exc),
            protocol_sha256=protocol_hash, code_sha256=code_hashes))
        raise


if __name__ == "__main__":
    main()
