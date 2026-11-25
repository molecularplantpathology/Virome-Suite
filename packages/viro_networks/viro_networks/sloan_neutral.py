"""Domain-Clipped Sloan Neutral Community Model Solver (`viro_networks.sloan_neutral`)."""

from __future__ import annotations

import dataclasses
import math
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import betainc
from scipy.stats import binomtest


@dataclasses.dataclass
class SloanTaxonFit:
    """Neutral community occurrence-abundance fit for an individual viral taxon."""

    taxon: str
    mean_rel_abundance: float
    observed_occurrence: float
    predicted_occurrence: float
    wilson_95_low: float
    wilson_95_high: float
    p_value: float
    q_value: float
    partition_class: str


def _wilson_interval(p_pred: float, n_samples: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    p = min(max(float(p_pred), 1e-9), 1.0 - 1e-9)
    n = float(max(1, n_samples))
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) / n) + (z2 / (4.0 * n * n)))) / denom
    return max(0.0, center - half), min(1.0, center + half)


def sloan_predicted_occurrence(
    p_vec: np.ndarray,
    ne_m: float,
    detection_limit: float,
) -> np.ndarray:
    """Domain-clipped Sloan (2006) incomplete Beta prediction (`P0-VN-1`)."""
    p_c = np.clip(np.asarray(p_vec, dtype=float), 1e-9, 1.0 - 1e-9)
    d_c = min(max(float(detection_limit), 1e-9), 1.0 - 1e-9)
    nm = min(max(float(ne_m), 1e-3), 1e4)
    a = np.maximum(2.0 * nm * p_c, 1e-9)
    b = np.maximum(2.0 * nm * (1.0 - p_c), 1e-9)
    return 1.0 - betainc(a, b, d_c)


def fit_sloan_neutral_model(
    counts: np.ndarray,
    taxa_names: Sequence[str],
    alpha_fdr: float = 0.05,
) -> Tuple[List[SloanTaxonFit], Dict[str, float]]:
    """Fit Sloan Neutral Community Model on untransformed relative abundances (`P0-VN-1`)."""
    mat = np.asarray(counts, dtype=float)
    n_samples, _p_taxa = mat.shape
    if n_samples < 5:
        raise ValueError("Sloan neutral model requires at least 5 samples.")

    row_totals = np.sum(mat, axis=1)
    valid_rows = row_totals > 0
    if int(np.sum(valid_rows)) < 5:
        raise ValueError("Fewer than 5 non-empty samples for Sloan neutral model.")

    mat_valid = mat[valid_rows, :]
    totals_valid = row_totals[valid_rows]
    n_eff = mat_valid.shape[0]

    # Retain taxa observed at least once across the cohort
    present_mask = np.sum(mat_valid > 0, axis=0) > 0
    kept_idx = [i for i, flag in enumerate(present_mask) if flag]
    if not kept_idx:
        return [], {"Ne_m": 0.0, "r_squared": 0.0, "detection_limit": 0.0}

    sub_counts = mat_valid[:, kept_idx]
    sub_names = [taxa_names[i] for i in kept_idx]

    rel_abund = sub_counts / totals_valid[:, None]
    mean_p = np.clip(np.mean(rel_abund, axis=0), 1e-9, 1.0 - 1e-9)
    obs_f = np.mean(sub_counts > 0, axis=0)
    d_lim = min(max(1.0 / float(np.median(totals_valid)), 1e-9), 0.5)

    def sse_obj(nm_val: float) -> float:
        pred = sloan_predicted_occurrence(mean_p, nm_val, d_lim)
        return float(np.sum((obs_f - pred) ** 2))

    opt = minimize_scalar(sse_obj, bounds=(1e-3, 1e4), method="bounded")
    ne_m_hat = float(opt.x)
    pred_f = sloan_predicted_occurrence(mean_p, ne_m_hat, d_lim)

    # Zero-variance R^2 guard (`P0-VN-1`)
    ss_res = float(np.sum((obs_f - pred_f) ** 2))
    f_bar = float(np.mean(obs_f))
    ss_tot = float(np.sum((obs_f - f_bar) ** 2))
    p_len = float(len(obs_f))
    if ss_tot >= 1e-12:
        r2 = max(-1.0, min(1.0, 1.0 - (ss_res / ss_tot)))
    else:
        denom_guard = p_len * max(f_bar * (1.0 - f_bar), 0.05)
        r2 = max(0.0, min(1.0, 1.0 - (ss_res / denom_guard)))

    # Binomial departure p-values + Benjamini-Hochberg adjustment
    raw_pvals: List[float] = []
    obs_k_list = np.sum(sub_counts > 0, axis=0).astype(int)
    for k_obs, pf in zip(obs_k_list, pred_f):
        pf_c = min(max(float(pf), 1e-9), 1.0 - 1e-9)
        pv = float(binomtest(int(k_obs), n=n_eff, p=pf_c).pvalue)
        raw_pvals.append(pv)

    m_taxa = len(raw_pvals)
    order = sorted(range(m_taxa), key=lambda i: raw_pvals[i])
    q_vals = [1.0] * m_taxa
    running = 1.0
    for rank_back in range(m_taxa - 1, -1, -1):
        orig_i = order[rank_back]
        cand_q = (raw_pvals[orig_i] * float(m_taxa)) / float(rank_back + 1)
        running = min(running, cand_q)
        q_vals[orig_i] = running

    results: List[SloanTaxonFit] = []
    for idx, name in enumerate(sub_names):
        w_lo, w_hi = _wilson_interval(float(pred_f[idx]), n_eff)
        of = float(obs_f[idx])
        pf = float(pred_f[idx])
        qv = q_vals[idx]
        if qv <= alpha_fdr and of > w_hi:
            p_class = "candidate_above_neutral_partition"
        elif qv <= alpha_fdr and of < w_lo:
            p_class = "candidate_below_neutral_partition"
        else:
            p_class = "neutral_drift_consistent"

        results.append(
            SloanTaxonFit(
                taxon=name,
                mean_rel_abundance=float(mean_p[idx]),
                observed_occurrence=of,
                predicted_occurrence=pf,
                wilson_95_low=w_lo,
                wilson_95_high=w_hi,
                p_value=raw_pvals[idx],
                q_value=qv,
                partition_class=p_class,
            )
        )

    summary = {
        "Ne_m": ne_m_hat,
        "r_squared": r2,
        "detection_limit": d_lim,
        "n_samples_used": n_eff,
        "n_taxa_evaluated": m_taxa,
    }
    return results, summary
