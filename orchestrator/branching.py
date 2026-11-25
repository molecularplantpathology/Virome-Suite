"""Morphology and Topology Dispatch Router for Virome_Suite."""

from __future__ import annotations

import math
from typing import Optional, Tuple

from orchestrator.schemas import (
    ContigEvidenceRecord,
    TopologyState,
    compute_coverage_uniformity,
    compute_host_flank_ratio,
    expected_spanning_pairs,
    expected_split_reads,
)


def compute_e_junction_score(
    n_split: int,
    n_spanning_pairs: int,
    d_origin: float,
    d_median: float,
    landmark_score: float,
    read_length: int = 150,
    mean_insert_size: float = 350.0,
    k_min_anchor: int = 20,
    d_min: float = 10.0,
) -> Tuple[float, dict]:
    """Compute the Physical Circular Junction Evidence Score (E_junction in [0, 100]).

    Uses read-geometry-scaled densities for split reads and spanning pairs,
    a low-depth guard (D_median >= 10x) on coverage uniformity, and a logistic model.
    """
    exp_split = expected_split_reads(d_median, read_length=read_length, k_min_anchor=k_min_anchor)
    exp_span = expected_spanning_pairs(
        d_median, read_length=read_length, mean_insert_size=mean_insert_size
    )

    s_split = min(1.0, float(n_split) / max(1.0, exp_split)) if d_median >= d_min else 0.0
    s_span = min(1.0, float(n_spanning_pairs) / max(1.0, exp_span)) if d_median >= d_min else 0.0
    s_cov = compute_coverage_uniformity(d_origin=d_origin, d_median=d_median, d_min=d_min)
    s_landmark = max(0.0, min(1.0, float(landmark_score)))

    # Calibrated logistic weights: intercept b0 = -2.2, positive weights on junction evidence
    b0 = -2.20
    w_split = 2.10
    w_span = 1.60
    w_cov = 1.20
    w_landmark = 0.90

    logit = b0 + (w_split * s_split) + (w_span * s_span) + (w_cov * s_cov) + (w_landmark * s_landmark)
    prob = 1.0 / (1.0 + math.exp(-logit))
    e_score = 100.0 * prob

    features = {
        "s_split": round(s_split, 4),
        "s_span": round(s_span, 4),
        "s_cov": round(s_cov, 4),
        "s_landmark": round(s_landmark, 4),
        "expected_split_reads": round(exp_split, 2),
        "expected_spanning_pairs": round(exp_span, 2),
    }
    return e_score, features


def classify_topology_route(
    e_junction_score: float,
    d_median: float,
    p_flank: Optional[float],
    d_min: float = 10.0,
    flank_conflict_threshold: float = 0.15,
    has_tir: bool = False,
) -> TopologyState:
    """Route a candidate contig to its appropriate pipeline branch with explicit Precedence Rules.

    Precedence Rules:
    0. Terminal Inverted Repeat (TIR) Veto -> LINEAR_CANDIDATE
    1. Low depth (<10x) or unassessed flank data -> UNRESOLVED_TOPOLOGY (Abstain)
    2. Host flank evidence >= 0.15:
       - If E_junction >= 70 -> MIXED_OR_CONFLICTED (both circular and host-flank signals)
       - Else -> INTEGRATED_CANDIDATE (EVE / chromosomal integration)
    3. High Circular Junction Support: E_junction >= 70, D_median >= 10x, P_flank < 0.15
       -> CIRCULAR_CANDIDATE
    4. Intermediate score (30 <= E_junction < 70) -> UNRESOLVED_TOPOLOGY (Abstain)
    5. Low circular evidence (E_junction < 30) and low flank (<0.15) -> LINEAR_CANDIDATE
    """
    if has_tir:
        return TopologyState.LINEAR_CANDIDATE

    if d_median < d_min or p_flank is None:
        return TopologyState.UNRESOLVED_TOPOLOGY

    if p_flank >= flank_conflict_threshold:
        if e_junction_score >= 70.0:
            return TopologyState.MIXED_OR_CONFLICTED
        return TopologyState.INTEGRATED_CANDIDATE

    if e_junction_score >= 70.0:
        return TopologyState.CIRCULAR_CANDIDATE

    if e_junction_score >= 30.0:
        return TopologyState.UNRESOLVED_TOPOLOGY

    return TopologyState.LINEAR_CANDIDATE


def build_contig_evidence_record(
    contig_id: str,
    length_bp: int,
    d_median: float,
    d_origin: float,
    n_split_reads: int,
    n_spanning_pairs: int,
    n_host_viral_chimeras: int,
    landmark_score: float,
    landmark_name: str,
    raw_sha256: str,
    curated_sha256: str,
    read_length: int = 150,
    mean_insert_size: float = 350.0,
) -> ContigEvidenceRecord:
    """Evaluate a candidate contig and package its routing decision into a ContigEvidenceRecord."""
    e_score, _ = compute_e_junction_score(
        n_split=n_split_reads,
        n_spanning_pairs=n_spanning_pairs,
        d_origin=d_origin,
        d_median=d_median,
        landmark_score=landmark_score,
        read_length=read_length,
        mean_insert_size=mean_insert_size,
    )
    p_flank = compute_host_flank_ratio(
        n_split=n_split_reads,
        n_host_viral_chimera=n_host_viral_chimeras,
    )
    state = classify_topology_route(
        e_junction_score=e_score,
        d_median=d_median,
        p_flank=p_flank,
    )
    return ContigEvidenceRecord(
        contig_id=contig_id,
        length_bp=length_bp,
        d_median=d_median,
        d_origin=d_origin,
        n_split_reads=n_split_reads,
        n_spanning_pairs=n_spanning_pairs,
        n_host_viral_chimeras=n_host_viral_chimeras,
        landmark_score=landmark_score,
        landmark_name=landmark_name,
        e_junction_score=e_score,
        p_flank=p_flank,
        topology_state=state,
        raw_sha256=raw_sha256,
        curated_sha256=curated_sha256,
    )
