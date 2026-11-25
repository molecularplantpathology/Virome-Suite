"""Compositional CLR Normalization & Three-Estimator Bootstrap Consensus Network (`viro_networks.clr_consensus`)."""

from __future__ import annotations

import dataclasses
import math
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.stats import spearmanr


class CohortInsufficientError(ValueError):
    """Raised when sample count N < 5 requires safe abstention (Exit Code 4)."""


@dataclasses.dataclass
class ConsensusEdge:
    """Validated non-random candidate co-occurrence edge."""

    taxon_a: str
    taxon_b: str
    sign: int  # +1 or -1
    sparcc_r: float
    partial_omega: float
    spearman_rho: float
    estimators_supporting: int
    bootstrap_stability: float
    pseudocount_concordance: float


def clr_transform(counts: np.ndarray, pseudocount: float = 0.5) -> np.ndarray:
    """Apply pseudocount zero-replacement and Centered Log-Ratio (CLR) per sample row."""
    c_adj = np.asarray(counts, dtype=float) + float(pseudocount)
    log_c = np.log(np.maximum(c_adj, 1e-12))
    return log_c - np.mean(log_c, axis=1, keepdims=True)


def compute_sparcc_matrix(counts: np.ndarray, pseudocount: float = 0.5) -> np.ndarray:
    """Positive-floored SparCC compositional correlation estimator (`P0-VN-2`)."""
    c_adj = np.asarray(counts, dtype=float) + float(pseudocount)
    log_c = np.log(np.maximum(c_adj, 1e-12))
    clr_mat = log_c - np.mean(log_c, axis=1, keepdims=True)

    n_samples, p_taxa = log_c.shape
    if p_taxa < 3 or n_samples < 2:
        corr = np.corrcoef(clr_mat, rowvar=False)
        return np.nan_to_num(corr, nan=0.0)

    # Aitchison log-ratio variation matrix t_ij = Var(ln(x_i / x_j))
    diff = log_c[:, :, None] - log_c[:, None, :]
    t_mat = np.var(diff, axis=0, ddof=1)
    t_row_sums = np.sum(t_mat, axis=1)
    total_t = float(np.sum(t_row_sums))

    clr_vars = np.var(clr_mat, axis=0, ddof=1)
    raw_omega2 = (t_row_sums - (total_t / (2.0 * float(p_taxa - 1)))) / float(p_taxa - 2)
    # Positive floor prevents negative basis variance -> NaN (`P0-VN-2`)
    omega2 = np.maximum(np.maximum(0.10 * clr_vars, 1e-6), raw_omega2)
    omega = np.sqrt(omega2)

    denom = 2.0 * np.outer(omega, omega)
    r_sparcc = (omega2[:, None] + omega2[None, :] - t_mat) / np.maximum(denom, 1e-12)
    r_sparcc = np.clip(r_sparcc, -1.0, 1.0)
    np.fill_diagonal(r_sparcc, 1.0)
    return r_sparcc


def compute_shrinkage_partial_corr(clr_mat: np.ndarray, alpha: float = 0.15) -> np.ndarray:
    """Correlation-scaled shrinkage precision matrix partial correlation (`P1-VN-3`)."""
    p_taxa = clr_mat.shape[1]
    if p_taxa < 2:
        return np.ones((p_taxa, p_taxa), dtype=float)

    stds = np.std(clr_mat, axis=0, ddof=1)
    safe_stds = np.where(stds > 1e-12, stds, 1.0)
    z_mat = (clr_mat - np.mean(clr_mat, axis=0, keepdims=True)) / safe_stds
    r_corr = (z_mat.T @ z_mat) / float(max(1, clr_mat.shape[0] - 1))
    np.fill_diagonal(r_corr, 1.0)

    shrunk = (1.0 - alpha) * r_corr + alpha * np.eye(p_taxa)
    theta = np.linalg.pinv(shrunk)
    diag_sqrt = np.sqrt(np.maximum(np.diag(theta), 1e-12))
    omega = -theta / np.outer(diag_sqrt, diag_sqrt)
    omega = np.clip(omega, -1.0, 1.0)
    np.fill_diagonal(omega, 1.0)
    return omega


def compute_spearman_clr(clr_mat: np.ndarray) -> np.ndarray:
    """Rank-based Spearman correlation matrix on CLR coordinates."""
    p_taxa = clr_mat.shape[1]
    if p_taxa < 2:
        return np.ones((p_taxa, p_taxa), dtype=float)
    res = spearmanr(clr_mat, axis=0).statistic
    if np.ndim(res) == 0:
        mat = np.array([[1.0, float(res)], [float(res), 1.0]])
    else:
        mat = np.asarray(res, dtype=float)
    mat = np.nan_to_num(mat, nan=0.0)
    np.fill_diagonal(mat, 1.0)
    return mat


def _ensemble_vote_matrix(
    counts: np.ndarray,
    pseudocount: float = 0.5,
    thresh_sparcc: float = 0.35,
    thresh_partial: float = 0.25,
    thresh_rank: float = 0.35,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (signed_support_matrix, vote_count_matrix, r_sparcc, omega_part, rho_rank)."""
    clr_mat = clr_transform(counts, pseudocount=pseudocount)
    r_sp = compute_sparcc_matrix(counts, pseudocount=pseudocount)
    om_p = compute_shrinkage_partial_corr(clr_mat)
    rh_r = compute_spearman_clr(clr_mat)

    s1 = np.where(np.abs(r_sp) >= thresh_sparcc, np.sign(r_sp), 0.0)
    s2 = np.where(np.abs(om_p) >= thresh_partial, np.sign(om_p), 0.0)
    s3 = np.where(np.abs(rh_r) >= thresh_rank, np.sign(rh_r), 0.0)

    pos_votes = (s1 > 0).astype(int) + (s2 > 0).astype(int) + (s3 > 0).astype(int)
    neg_votes = (s1 < 0).astype(int) + (s2 < 0).astype(int) + (s3 < 0).astype(int)

    signed_consensus = np.zeros_like(r_sp, dtype=int)
    vote_counts = np.zeros_like(r_sp, dtype=int)

    pos_mask = pos_votes >= 2
    neg_mask = neg_votes >= 2
    signed_consensus[pos_mask] = 1
    signed_consensus[neg_mask] = -1
    vote_counts[pos_mask] = pos_votes[pos_mask]
    vote_counts[neg_mask] = neg_votes[neg_mask]
    np.fill_diagonal(signed_consensus, 0)
    np.fill_diagonal(vote_counts, 0)
    return signed_consensus, vote_counts, r_sp, om_p, rh_r


def infer_consensus_network(
    counts: np.ndarray,
    taxa_names: Sequence[str],
    n_boot: int = 100,
    stability_threshold: float = 0.80,
    pseudocount_grid: Tuple[float, ...] = (0.1, 0.5, 1.0),
    seed: int = 42,
) -> Tuple[List[ConsensusEdge], Dict[str, float]]:
    """Infer bootstrap-stabilized 3-estimator consensus co-occurrence network."""
    mat = np.asarray(counts, dtype=float)
    n_samples, p_taxa = mat.shape
    if n_samples < 5:
        raise CohortInsufficientError(
            f"Cohort has N={n_samples} samples (< 5 minimum required); abstaining safely."
        )

    # Filter taxa present in >= max(2, ceil(0.10 * N)) samples
    min_prev = max(2, int(math.ceil(0.10 * n_samples)))
    keep_mask = np.sum(mat > 0, axis=0) >= min_prev
    kept_indices = [i for i, k in enumerate(keep_mask) if k]
    if len(kept_indices) < 2:
        return [], {"retained_taxa": len(kept_indices), "pseudocount_grid_jaccard": 1.0}

    sub_mat = mat[:, kept_indices]
    sub_names = [taxa_names[i] for i in kept_indices]
    p_kept = len(sub_names)

    # Primary fit at default pseudocount = 0.5
    base_sign, base_votes, r_sp, om_p, rh_r = _ensemble_vote_matrix(sub_mat, pseudocount=0.5)

    # Pseudocount sensitivity evaluation (`T5.2`)
    grid_signs = {}
    for pc in pseudocount_grid:
        gs, _, _, _, _ = _ensemble_vote_matrix(sub_mat, pseudocount=pc)
        grid_signs[pc] = gs

    # Sample-level bootstrap resampling (`B = 100`)
    rng = np.random.default_rng(seed)
    pos_boot_counts = np.zeros((p_kept, p_kept), dtype=int)
    neg_boot_counts = np.zeros((p_kept, p_kept), dtype=int)
    for _ in range(n_boot):
        b_idx = rng.integers(0, n_samples, size=n_samples)
        b_sign, _, _, _, _ = _ensemble_vote_matrix(sub_mat[b_idx, :], pseudocount=0.5)
        pos_boot_counts += (b_sign == 1).astype(int)
        neg_boot_counts += (b_sign == -1).astype(int)

    edges: List[ConsensusEdge] = []
    for i in range(p_kept):
        for j in range(i + 1, p_kept):
            sgn = int(base_sign[i, j])
            if sgn == 0:
                continue
            stab = (
                float(pos_boot_counts[i, j]) / float(n_boot)
                if sgn == 1
                else float(neg_boot_counts[i, j]) / float(n_boot)
            )
            if stab < stability_threshold:
                continue
            agree_pc = sum(1 for pc in pseudocount_grid if int(grid_signs[pc][i, j]) == sgn)
            pc_conc = float(agree_pc) / float(len(pseudocount_grid))
            edges.append(
                ConsensusEdge(
                    taxon_a=sub_names[i],
                    taxon_b=sub_names[j],
                    sign=sgn,
                    sparcc_r=float(r_sp[i, j]),
                    partial_omega=float(om_p[i, j]),
                    spearman_rho=float(rh_r[i, j]),
                    estimators_supporting=int(base_votes[i, j]),
                    bootstrap_stability=stab,
                    pseudocount_concordance=pc_conc,
                )
            )

    # Compute overall Jaccard concordance between pseudocount extremes
    e_lo = set(zip(*np.where(np.triu(grid_signs[pseudocount_grid[0]], k=1) != 0)))
    e_hi = set(zip(*np.where(np.triu(grid_signs[pseudocount_grid[-1]], k=1) != 0)))
    union_len = len(e_lo | e_hi)
    jacc = float(len(e_lo & e_hi)) / float(union_len) if union_len > 0 else 1.0

    diagnostics = {
        "n_samples": n_samples,
        "retained_taxa": p_kept,
        "consensus_edge_count": len(edges),
        "pseudocount_extremes_jaccard": round(jacc, 4),
    }
    return edges, diagnostics
