"""Auditable ranking diagnostics. No causal interpretation or coverage guarantee.

Reference membership means supported evidence, not that its complement is a
biologically validated negative. Gene-level randomization is a conditional
reference model; dependence between genes can invalidate its inferential use.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def topk_weights(scores, k):
    """Expected selection indicators under uniform tie-breaking at the boundary."""
    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 1 or not len(scores) or not np.isfinite(scores).all():
        raise ValueError("scores must be a nonempty finite vector")
    if isinstance(k, bool) or int(k) != k or not 0 < k <= len(scores):
        raise ValueError("k must be an integer in [1, scope size]")
    boundary = np.partition(scores, len(scores) - int(k))[len(scores) - int(k)]
    weights = (scores > boundary).astype(float)
    tied = scores == boundary
    weights[tied] = (int(k) - weights.sum()) / tied.sum()
    return weights


def rank_metrics(scores, supported, k, total_supported=None):
    y = np.asarray(supported, dtype=bool)
    w = topk_weights(scores, k)
    if y.shape != w.shape:
        raise ValueError("scores and reference labels have different shapes")
    m, n = int(y.sum()), len(y)
    M = m if total_supported is None else total_supported
    if int(M) != M or M < m:
        raise ValueError("total reference count must be an integer >= in-scope count")
    hits = float(w @ y)
    return {
        "n_scope": n, "n_supported_in_scope": m, "n_supported_total": int(M),
        "k": int(k), "expected_hits": hits, "precision": hits / k,
        "scope_enrichment": hits * n / (k * m) if m else np.nan,
        "full_reference_recall": hits / M if M else np.nan,
        "reference_coverage": m / M if M else np.nan,
        "rankable": bool(np.ptp(scores) > 0),
    }


def conditional_moments(weights, labels, strata=None):
    """Exact moments of weighted hits when labels permute within fixed strata.

    These moments describe the permutation reference, not a biological null.
    Zero variance includes constant weights within every stratum: such a
    ranking cannot supply within-stratum evidence.
    """
    w, y = np.asarray(weights, float), np.asarray(labels, bool)
    s = np.zeros(len(w), dtype=int) if strata is None else np.asarray(strata)
    if w.shape != y.shape or s.shape != w.shape or not np.isfinite(w).all():
        raise ValueError("weights, labels, and strata must align and be finite")
    expected = variance = 0.0
    for group in np.unique(s):
        keep = s == group
        wi, yi = w[keep], y[keep]
        n, m = len(wi), int(yi.sum())
        expected += wi.sum() * m / n
        if n > 1:
            variance += m * (n - m) / (n * (n - 1)) * np.square(wi - wi.mean()).sum()
    observed = float(w @ y)
    excess = observed - expected
    return {
        "reference_mean_hits": float(expected),
        "reference_variance_hits": float(variance),
        "excess_hits": float(excess),
        "standardized_excess": float(excess / np.sqrt(variance)) if variance > 1e-14 else 0.0,
        "within_stratum_information": bool(variance > 1e-14),
    }


def permutation_pvalue(weights, labels, strata, rng, n_permutations=999):
    """One-sided Monte Carlo conditional randomization p, plus-one correction."""
    w, y, s = np.asarray(weights, float), np.asarray(labels, bool), np.asarray(strata)
    moments = conditional_moments(w, y, s)
    if not moments["within_stratum_information"]:
        return 1.0
    if n_permutations < 1:
        raise ValueError("at least one permutation required")
    observed = float(w @ y)
    values = np.zeros(n_permutations)
    for group in np.unique(s):
        idx = np.flatnonzero(s == group)
        for b in range(n_permutations):
            values[b] += w[idx] @ rng.permutation(y[idx])
    return float((1 + (values >= observed - 1e-12).sum()) / (n_permutations + 1))


def quantile_strata(covariates, n_bins=5):
    """Marginal quantile grid; identical values stay together, no label access."""
    x = np.asarray(covariates, float)
    if x.ndim == 1:
        x = x[:, None]
    if not np.isfinite(x).all():
        raise ValueError("missing covariates require explicit treatment")
    codes = np.zeros(len(x), dtype=np.int64)
    for j in range(x.shape[1]):
        edges = np.unique(np.quantile(x[:, j], np.arange(1, n_bins) / n_bins))
        codes = codes * (n_bins + 1) + np.searchsorted(edges, x[:, j], side="right")
    return np.unique(codes, return_inverse=True)[1]


def residualize(response, covariates):
    y, x = np.asarray(response, float), np.asarray(covariates, float)
    if x.ndim == 1:
        x = x[:, None]
    if len(x) != len(y) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("regression inputs must align and be finite")
    sd = x.std(axis=0)
    x = (x[:, sd > 1e-12] - x[:, sd > 1e-12].mean(axis=0)) / sd[sd > 1e-12]
    design = np.column_stack([np.ones(len(y)), x])
    resid = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    # A mathematically constant residual must not become a numerical ranking.
    if np.ptp(resid) < 1e-12 * max(1.0, np.max(np.abs(y))):
        resid[:] = 0
    return resid


def bh_adjust(pvalues):
    p = np.asarray(pvalues, float)
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be finite and in [0, 1]")
    order = np.argsort(p)
    q_sorted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.minimum(q_sorted, 1)
    return q


def welch_reference(a, b, alpha=0.05, min_abs_delta=0.1):
    """Cell-split exploratory reference, not biological replicate inference."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    delta = a.mean(axis=0) - b.mean(axis=0)
    va, vb = a.var(axis=0, ddof=1) / len(a), b.var(axis=0, ddof=1) / len(b)
    v = va + vb
    denom = va * va / (len(a) - 1) + vb * vb / (len(b) - 1)
    df = np.divide(v * v, denom, out=np.ones_like(v), where=denom > 0)
    t = np.divide(delta, np.sqrt(v), out=np.zeros_like(v), where=v > 0)
    p = 2 * stats.t.sf(np.abs(t), df)
    # Deterministic nonzero shifts with zero variance are explicitly supported.
    p[(v == 0) & (delta != 0)] = 0
    q = bh_adjust(p)
    return delta, q, (q <= alpha) & (np.abs(delta) >= min_abs_delta)


def bootstrap_mean(values, rng, n_bootstrap=2000):
    """Resample rows with multiplicity. Caller defines independent sampling unit."""
    x = np.asarray(values, float)
    if not len(x) or not np.isfinite(x).all():
        raise ValueError("finite nonempty values required")
    means = x[rng.integers(0, len(x), size=(n_bootstrap, len(x)))].mean(axis=1)
    return np.quantile(means, [0.025, 0.975])
