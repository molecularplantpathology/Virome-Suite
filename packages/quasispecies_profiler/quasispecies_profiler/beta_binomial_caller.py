"""Hierarchical Empirical Bayes Beta-Binomial Variant Caller (`quasispecies_profiler.beta_binomial_caller`)."""

from __future__ import annotations

import dataclasses
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import betaln, gammaln
from scipy.stats import chi2


@dataclasses.dataclass
class VariantCallRecord:
    """Evaluated intra-host candidate site output."""

    pos: int
    ref: str
    alt: str
    mean_depth: float
    vaf_hat: float
    rho_hat: float
    ci_95_low: float
    ci_95_high: float
    lrt_stat: float
    p_value: float
    q_value: float
    status: str
    is_called: bool
    replicate_vafs: List[float]


def log_beta_binomial_pmf(
    k: np.ndarray,
    n: np.ndarray,
    mu: float,
    rho: float,
) -> np.ndarray:
    """Numerically stable log-PMF of Beta-Binomial(n, mu, rho) in pure log-space (`EC-3.1`)."""
    mu_c = min(max(float(mu), 1e-7), 1.0 - 1e-7)
    rho_c = min(max(float(rho), 1e-7), 0.5)
    phi = (1.0 / rho_c) - 1.0
    alpha = mu_c * phi
    beta_p = (1.0 - mu_c) * phi
    log_comb = gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)
    return log_comb + betaln(k + alpha, n - k + beta_p) - betaln(alpha, beta_p)


def penalized_log_likelihood(
    mu: float,
    rho: float,
    k: np.ndarray,
    n: np.ndarray,
    mu_rho_prior: float,
    sigma_rho_prior: float,
) -> float:
    """Penalized Empirical Bayes Beta-Binomial log-likelihood."""
    ll = float(np.sum(log_beta_binomial_pmf(k, n, mu, rho)))
    penalty = ((math.log(max(rho, 1e-7)) - mu_rho_prior) ** 2) / (2.0 * (sigma_rho_prior ** 2))
    return ll - penalty


def estimate_background_error_and_prior(
    sites: Sequence[Dict[str, Any]],
    max_bg_vaf: float = 0.002,
) -> Tuple[Dict[str, float], float, float, float]:
    """Select monomorphic background sites via MAX replicate VAF < 0.002 (`P0-QP-1`).

    Returns (context_eps_map, global_eps_bg, mu_rho_prior, sigma_rho_prior).
    """
    bg_sites: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for s in sites:
        k_arr = np.asarray(s["k"], dtype=float)
        n_arr = np.asarray(s["n"], dtype=float)
        if np.any(n_arr <= 0):
            continue
        vafs = k_arr / n_arr
        if float(np.max(vafs)) < max_bg_vaf:
            sub_ctx = f"{s.get('ref', 'N')}>{s.get('alt', 'N')}"
            bg_sites.append((sub_ctx, k_arr, n_arr))

    if not bg_sites:
        return {}, 5e-4, math.log(1e-4), 3.5

    total_k = sum(float(np.sum(k)) for _, k, _ in bg_sites)
    total_n = sum(float(np.sum(n)) for _, _, n in bg_sites)
    global_eps = max(5e-4, total_k / max(1.0, total_n))

    ctx_counts: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = {}
    log_rhos: List[float] = []
    for ctx, k_arr, n_arr in bg_sites:
        ctx_counts.setdefault(ctx, []).append((k_arr, n_arr))
        # Method-of-moments overdispersion estimate bounded to [1e-6, 0.5]
        p_bar = float(np.sum(k_arr)) / max(1.0, float(np.sum(n_arr)))
        n_bar = float(np.mean(n_arr))
        if p_bar > 0 and n_bar > 1.0 and len(k_arr) >= 2:
            s2 = float(np.var(k_arr, ddof=1))
            bin_var = n_bar * p_bar * (1.0 - p_bar)
            rho_mm = (s2 - bin_var) / max(1e-9, bin_var * (n_bar - 1.0))
            rho_est = min(max(rho_mm, 1e-6), 0.5)
        else:
            rho_est = 1e-5
        log_rhos.append(math.log(rho_est))

    ctx_eps_map: Dict[str, float] = {}
    for ctx, items in ctx_counts.items():
        if len(items) >= 5:
            ck = sum(float(np.sum(x[0])) for x in items)
            cn = sum(float(np.sum(x[1])) for x in items)
            ctx_eps_map[ctx] = max(5e-4, ck / max(1.0, cn))
        else:
            ctx_eps_map[ctx] = global_eps

    mu_rho = float(np.mean(log_rhos)) if log_rhos else math.log(1e-4)
    raw_std = float(np.std(log_rhos, ddof=1)) if len(log_rhos) >= 2 else 0.0
    # Weakly informative diffuse log-dispersion prior scale (`P0-QP-2`)
    sigma_rho = max(3.5, raw_std)
    return ctx_eps_map, global_eps, mu_rho, sigma_rho


def _fit_penalized_mle(
    k: np.ndarray,
    n: np.ndarray,
    mu_rho_prior: float,
    sigma_rho_prior: float,
) -> Tuple[float, float, float]:
    """Bounded L-BFGS-B + multi-start initialization for (mu_hat, rho_hat, max_ll)."""
    p_init = min(max(float(np.sum(k)) / max(1.0, float(np.sum(n))), 1e-5), 1.0 - 1e-5)

    def neg_obj(params: np.ndarray) -> float:
        mu_val = float(params[0])
        rho_val = math.exp(float(params[1]))
        return -penalized_log_likelihood(mu_val, rho_val, k, n, mu_rho_prior, sigma_rho_prior)

    best_mu = p_init
    best_rho = math.exp(mu_rho_prior)
    best_ll = penalized_log_likelihood(best_mu, best_rho, k, n, mu_rho_prior, sigma_rho_prior)

    for rho_start in (1e-5, 1e-3, 0.02):
        res = minimize(
            neg_obj,
            x0=np.array([p_init, math.log(rho_start)]),
            method="L-BFGS-B",
            bounds=[(1e-5, 1.0 - 1e-5), (math.log(1e-6), math.log(0.5))],
        )
        cand_mu = float(res.x[0])
        cand_rho = math.exp(float(res.x[1]))
        cand_ll = -float(res.fun)
        if cand_ll > best_ll:
            best_ll = cand_ll
            best_mu = cand_mu
            best_rho = cand_rho

    return best_mu, best_rho, best_ll


def _profile_ll_at_mu(
    mu_fixed: float,
    k: np.ndarray,
    n: np.ndarray,
    mu_rho_prior: float,
    sigma_rho_prior: float,
) -> Tuple[float, float]:
    """Maximize penalized log-likelihood over rho for a fixed mu."""
    res = minimize_scalar(
        lambda log_r: -penalized_log_likelihood(
            mu_fixed, math.exp(float(log_r)), k, n, mu_rho_prior, sigma_rho_prior
        ),
        bounds=(math.log(1e-6), math.log(0.5)),
        method="bounded",
    )
    rho_opt = math.exp(float(res.x))
    return -float(res.fun), rho_opt


def compute_profile_ci_95(
    mu_hat: float,
    max_ll: float,
    k: np.ndarray,
    n: np.ndarray,
    mu_rho_prior: float,
    sigma_rho_prior: float,
) -> Tuple[float, float]:
    """Compute sign-bracketed 95% Profile Likelihood Confidence Interval (`P1-QP-4`)."""
    cutoff = 3.841458820694124  # chi2_{1, 0.95}
    mu_min = 1e-5
    mu_max = 1.0 - 1e-5

    def diff_fn(m: float) -> float:
        pll, _ = _profile_ll_at_mu(m, k, n, mu_rho_prior, sigma_rho_prior)
        return 2.0 * (max_ll - pll) - cutoff

    # Lower bound in [mu_min, mu_hat]
    if mu_hat <= mu_min or diff_fn(mu_min) <= 0.0:
        ci_low = 0.0
    else:
        ci_low = float(brentq(diff_fn, mu_min, mu_hat, xtol=1e-6))

    # Upper bound in [mu_hat, mu_max]
    if mu_hat >= mu_max or diff_fn(mu_max) <= 0.0:
        ci_high = 1.0
    else:
        ci_high = float(brentq(diff_fn, mu_hat, mu_max, xtol=1e-6))

    return ci_low, ci_high


def _vectorized_null_bootstrap_pvalue(
    obs_lrt: float,
    n: np.ndarray,
    eps_bg: float,
    rho_null: float,
    n_boot: int = 2000,
    seed: int = 42,
) -> float:
    """Vectorized null parametric bootstrap p-value with asymptotic tail resolution (`P0-QP-3`)."""
    if obs_lrt <= 0.0:
        return 1.0

    rng = np.random.default_rng(seed)
    m_reps = len(n)
    rho_c = min(max(rho_null, 1e-6), 0.5)
    phi = (1.0 / rho_c) - 1.0
    a_null = max(eps_bg * phi, 1e-6)
    b_null = max((1.0 - eps_bg) * phi, 1e-6)

    # Simulate (n_boot, m_reps) Beta-Binomial null counts in one vectorized NumPy draw
    p_draws = rng.beta(a_null, b_null, size=(n_boot, m_reps))
    k_boot = rng.binomial(n.astype(int)[None, :], p_draws).astype(float)

    n_row = n[None, :]
    mu_boot = np.clip(np.sum(k_boot, axis=1) / np.sum(n_row, axis=1), eps_bg, 1.0 - 1e-5)

    # Evaluate vectorized Beta-Binomial LRT under null dispersion
    a_alt = mu_boot[:, None] * phi
    b_alt = (1.0 - mu_boot[:, None]) * phi
    ll_alt = np.sum(betaln(k_boot + a_alt, n_row - k_boot + b_alt) - betaln(a_alt, b_alt), axis=1)
    ll_0 = np.sum(betaln(k_boot + a_null, n_row - k_boot + b_null) - betaln(a_null, b_null), axis=1)
    lrt_boot = np.maximum(0.0, 2.0 * (ll_alt - ll_0))

    b_exceed = int(np.sum(lrt_boot >= obs_lrt))
    if b_exceed == 0 and obs_lrt > 15.0:
        # Extreme tail resolution calibrated by 0.5 * chi2_1 survival function capped by 1/(n_boot+1)
        asymp_p = 0.5 * float(chi2.sf(obs_lrt, df=1))
        return min(1.0 / float(n_boot + 1), max(asymp_p, 1e-16))

    return float(b_exceed + 1) / float(n_boot + 1)


def call_variants_across_replicates(
    sites: Sequence[Dict[str, Any]],
    replicate_type: str = "technical",
    alpha_fdr: float = 0.01,
    n_boot: int = 2000,
    seed: int = 42,
) -> List[VariantCallRecord]:
    """Execute full Empirical Bayes Beta-Binomial variant calling across M replicates."""
    ctx_eps_map, global_eps, mu_rho_prior, sigma_rho_prior = estimate_background_error_and_prior(sites)

    prelim_records: List[VariantCallRecord] = []
    evaluable_indices: List[int] = []

    for idx, site in enumerate(sites):
        pos = int(site["pos"])
        ref = str(site.get("ref", "N")).upper()
        alt = str(site.get("alt", "N")).upper()
        k = np.asarray(site["k"], dtype=float)
        n = np.asarray(site["n"], dtype=float)
        m_reps = len(k)
        mean_depth = float(np.mean(n)) if m_reps > 0 else 0.0
        rep_vafs = [float(kv / max(1.0, nv)) for kv, nv in zip(k, n)]
        pooled_vaf = float(np.sum(k)) / max(1.0, float(np.sum(n)))

        # Short-circuit zero alternate counts (`EC-3.4`)
        if float(np.sum(k)) <= 0.0:
            prelim_records.append(
                VariantCallRecord(
                    pos=pos, ref=ref, alt=alt, mean_depth=mean_depth,
                    vaf_hat=0.0, rho_hat=1e-6, ci_95_low=0.0, ci_95_high=0.0,
                    lrt_stat=0.0, p_value=1.0, q_value=1.0,
                    status="MONOMORPHIC_ZERO", is_called=False, replicate_vafs=rep_vafs,
                )
            )
            continue

        # Step 1: Depth-dependent Limit of Detection Pre-Gate (`EC-3.3`)
        if mean_depth < 100.0:
            prelim_records.append(
                VariantCallRecord(
                    pos=pos, ref=ref, alt=alt, mean_depth=mean_depth,
                    vaf_hat=pooled_vaf, rho_hat=1e-6, ci_95_low=0.0, ci_95_high=0.0,
                    lrt_stat=0.0, p_value=1.0, q_value=1.0,
                    status="UNASSESSED_LOW_DEPTH", is_called=False, replicate_vafs=rep_vafs,
                )
            )
            continue

        if mean_depth < 1000.0 and pooled_vaf < 0.030:
            prelim_records.append(
                VariantCallRecord(
                    pos=pos, ref=ref, alt=alt, mean_depth=mean_depth,
                    vaf_hat=pooled_vaf, rho_hat=1e-6, ci_95_low=0.0, ci_95_high=0.0,
                    lrt_stat=0.0, p_value=1.0, q_value=1.0,
                    status="BELOW_LOD_AT_DEPTH", is_called=False, replicate_vafs=rep_vafs,
                )
            )
            continue

        ctx_key = f"{ref}>{alt}"
        eps_bg = ctx_eps_map.get(ctx_key, global_eps)

        # Fit penalized Beta-Binomial MLE
        mu_hat, rho_hat, max_ll = _fit_penalized_mle(k, n, mu_rho_prior, sigma_rho_prior)

        # Pure monomorphic / sub-threshold background check (`P0-QP-1`)
        active_thresh = max(3.0 * eps_bg, 0.0015)
        m_active = sum(1 for v in rep_vafs if v > active_thresh)
        cv_vaf = (float(np.std(rep_vafs, ddof=1)) / max(1e-9, float(np.mean(rep_vafs)))) if m_reps >= 2 else 0.0

        if float(np.max(rep_vafs)) < 0.002 or (m_active == 0 and (mu_hat <= eps_bg or pooled_vaf < 0.005)):
            prelim_records.append(
                VariantCallRecord(
                    pos=pos, ref=ref, alt=alt, mean_depth=mean_depth,
                    vaf_hat=mu_hat, rho_hat=rho_hat, ci_95_low=0.0, ci_95_high=min(1.0, mu_hat * 2.0),
                    lrt_stat=0.0, p_value=1.0, q_value=1.0,
                    status="BACKGROUND_NOISE", is_called=False, replicate_vafs=rep_vafs,
                )
            )
            continue

        # Technical Replicate Jackpot Veto (`P0-QP-2`)
        if replicate_type == "technical" and m_reps >= 3 and (
            m_active <= 1 or (rho_hat > 0.005 and cv_vaf > 1.0)
        ):
            prelim_records.append(
                VariantCallRecord(
                    pos=pos, ref=ref, alt=alt, mean_depth=mean_depth,
                    vaf_hat=mu_hat, rho_hat=rho_hat, ci_95_low=0.0, ci_95_high=min(1.0, max(rep_vafs)),
                    lrt_stat=0.0, p_value=1.0, q_value=1.0,
                    status="REJECTED_PCR_JACKPOT", is_called=False, replicate_vafs=rep_vafs,
                )
            )
            continue

        if mu_hat <= eps_bg or pooled_vaf < 0.005:
            prelim_records.append(
                VariantCallRecord(
                    pos=pos, ref=ref, alt=alt, mean_depth=mean_depth,
                    vaf_hat=mu_hat, rho_hat=rho_hat, ci_95_low=0.0, ci_95_high=min(1.0, mu_hat * 2.0),
                    lrt_stat=0.0, p_value=1.0, q_value=1.0,
                    status="BACKGROUND_NOISE", is_called=False, replicate_vafs=rep_vafs,
                )
            )
            continue

        null_ll, rho_0 = _profile_ll_at_mu(eps_bg, k, n, mu_rho_prior, sigma_rho_prior)
        lrt_stat = max(0.0, 2.0 * (max_ll - null_ll))
        p_val = _vectorized_null_bootstrap_pvalue(
            obs_lrt=lrt_stat, n=n, eps_bg=eps_bg, rho_null=max(rho_0, math.exp(mu_rho_prior)),
            n_boot=n_boot, seed=seed + idx,
        )
        ci_low, ci_high = compute_profile_ci_95(mu_hat, max_ll, k, n, mu_rho_prior, sigma_rho_prior)

        evaluable_indices.append(idx)
        prelim_records.append(
            VariantCallRecord(
                pos=pos, ref=ref, alt=alt, mean_depth=mean_depth,
                vaf_hat=mu_hat, rho_hat=rho_hat, ci_95_low=ci_low, ci_95_high=ci_high,
                lrt_stat=lrt_stat, p_value=p_val, q_value=1.0,
                status="CANDIDATE_EVALUATED", is_called=False, replicate_vafs=rep_vafs,
            )
        )

    # Benjamini-Hochberg FDR multiplicity adjustment over evaluable candidate family (`P0-QP-3`)
    m_eval = len(evaluable_indices)
    if m_eval > 0:
        sorted_eval = sorted(evaluable_indices, key=lambda i: prelim_records[i].p_value)
        q_raw = [0.0] * m_eval
        for rank_idx, rec_idx in enumerate(sorted_eval, start=1):
            q_raw[rank_idx - 1] = (prelim_records[rec_idx].p_value * float(m_eval)) / float(rank_idx)
        # Enforce monotonicity from largest rank down
        running_min = 1.0
        for r_back in range(m_eval - 1, -1, -1):
            running_min = min(running_min, q_raw[r_back])
            rec_i = sorted_eval[r_back]
            prelim_records[rec_i].q_value = running_min
            if running_min <= alpha_fdr:
                prelim_records[rec_i].is_called = True
                prelim_records[rec_i].status = "PASS_ISNV"
            else:
                prelim_records[rec_i].status = "FILTERED_FDR"

    return prelim_records
