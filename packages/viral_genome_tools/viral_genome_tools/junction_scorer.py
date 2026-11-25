"""Read-Geometry Physical Circular Junction Scorer & Precedence Router (`viral_genome_tools.junction_scorer`)."""

from __future__ import annotations

import dataclasses
import enum
import math
from typing import Dict, Iterable, List, Optional, Set, Tuple

from viral_genome_tools.monomer_curator import revcomp


class TopologyDecision(str, enum.Enum):
    """Self-contained topology classification states for `viral-genome-tools`."""

    CIRCULAR_CANDIDATE = "circular_candidate"
    LINEAR_CANDIDATE = "linear_candidate"
    INTEGRATED_CANDIDATE = "integrated_candidate"
    MIXED_OR_CONFLICTED = "mixed_or_conflicted"
    UNRESOLVED_TOPOLOGY = "unresolved_topology"


def expected_split_reads(
    d_median: float,
    read_length: int = 150,
    k_min_anchor: int = 20,
) -> float:
    """E[N_split] = D_median * (1 - 2 * k_min / L_r)."""
    if d_median <= 0 or read_length <= 2 * k_min_anchor:
        return 0.0
    return float(d_median) * (1.0 - (2.0 * float(k_min_anchor)) / float(read_length))


def expected_spanning_pairs(
    d_median: float,
    read_length: int = 150,
    mean_insert_size: float = 350.0,
) -> float:
    """E[N_span] = (D_median / (2 * L_r)) * max(0, mu_ins - 2 * L_r)."""
    if d_median <= 0 or read_length <= 0:
        return 0.0
    window = max(0.0, float(mean_insert_size) - 2.0 * float(read_length))
    return (float(d_median) / (2.0 * float(read_length))) * window


def compute_coverage_uniformity(
    d_origin: float,
    d_median: float,
    d_min: float = 10.0,
    epsilon: float = 1.0,
) -> float:
    """Bounded coverage uniformity across circular origin with <10x low-depth guard."""
    if d_median < d_min:
        return 0.0
    ratio = (float(d_origin) + epsilon) / (float(d_median) + epsilon)
    return math.exp(-abs(math.log(ratio)))


def compute_host_flank_ratio(
    n_split: int,
    n_host_viral_chimera: int,
) -> Optional[float]:
    """Orthogonal host-integration evidence ratio P_flank (None when 0 informative reads)."""
    total = int(n_split) + int(n_host_viral_chimera)
    if total <= 0:
        return None
    return float(n_host_viral_chimera) / float(total)


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
) -> Tuple[float, Dict[str, float]]:
    """Compute full-precision Physical Circular Junction Evidence Score E_junction in [0, 100] (`P1-VGT-3`)."""
    exp_split = expected_split_reads(d_median, read_length=read_length, k_min_anchor=k_min_anchor)
    exp_span = expected_spanning_pairs(
        d_median, read_length=read_length, mean_insert_size=mean_insert_size
    )

    s_split = min(1.0, float(n_split) / max(1.0, exp_split)) if d_median >= d_min else 0.0
    s_span = min(1.0, float(n_spanning_pairs) / max(1.0, exp_span)) if d_median >= d_min else 0.0
    s_cov = compute_coverage_uniformity(d_origin=d_origin, d_median=d_median, d_min=d_min)
    s_landmark = max(0.0, min(1.0, float(landmark_score)))

    b0 = -2.20
    w_split = 2.10
    w_span = 1.60
    w_cov = 1.20
    w_landmark = 0.90

    logit = b0 + (w_split * s_split) + (w_span * s_span) + (w_cov * s_cov) + (w_landmark * s_landmark)
    prob = 1.0 / (1.0 + math.exp(-logit))
    e_score = 100.0 * prob

    features = {
        "s_split": s_split,
        "s_span": s_span,
        "s_cov": s_cov,
        "s_landmark": s_landmark,
        "expected_split_reads": exp_split,
        "expected_spanning_pairs": exp_span,
    }
    return e_score, features


def classify_topology(
    e_junction_score: float,
    d_median: float,
    p_flank: Optional[float],
    has_tir: bool = False,
    d_min: float = 10.0,
    flank_conflict_threshold: float = 0.15,
) -> TopologyDecision:
    """Classify candidate contig topology with Precedence Rule 0 (TIR Veto -> LINEAR_CANDIDATE)."""
    if has_tir:
        return TopologyDecision.LINEAR_CANDIDATE
    if d_median < d_min or p_flank is None:
        return TopologyDecision.UNRESOLVED_TOPOLOGY
    if p_flank >= flank_conflict_threshold:
        if e_junction_score >= 70.0:
            return TopologyDecision.MIXED_OR_CONFLICTED
        return TopologyDecision.INTEGRATED_CANDIDATE
    if e_junction_score >= 70.0:
        return TopologyDecision.CIRCULAR_CANDIDATE
    if e_junction_score >= 30.0:
        return TopologyDecision.UNRESOLVED_TOPOLOGY
    return TopologyDecision.LINEAR_CANDIDATE


@dataclasses.dataclass
class JunctionReadEvidence:
    """Measured read support across a proposed circular monomer junction."""

    d_median: float
    d_origin: float
    n_split: int
    n_spanning_pairs: int
    n_host_viral_chimera: int
    read_length: int


def profile_reads_across_seam(
    monomer_seq: str,
    read_pairs: Iterable[Tuple[str, str, Optional[str]]],
    host_kmers: Optional[Set[str]] = None,
    k_anchor: int = 20,
) -> JunctionReadEvidence:
    """Measure split reads, spanning pairs, depth, and host-viral chimeras directly from reads."""
    s = monomer_seq.upper()
    g_len = len(s)
    pairs = list(read_pairs)
    if not pairs or g_len < 2 * k_anchor:
        return JunctionReadEvidence(0.0, 0.0, 0, 0, 0, 150)

    doubled = s + s
    host_set = host_kmers or set()

    total_aligned_bases = 0
    origin_cov_count = 0
    n_split = 0
    n_span = 0
    n_chimera = 0
    read_lengths: List[int] = []

    def _locate_on_doubled(read_seq: str) -> Optional[Tuple[int, int]]:
        r = read_seq.upper()
        for candidate in (r, revcomp(r)):
            anchor = candidate[:k_anchor]
            pos = doubled.find(anchor)
            if pos != -1 and pos < g_len:
                end = pos + len(candidate)
                if end <= len(doubled):
                    # Verify >=90% match
                    m = sum(1 for a, b in zip(candidate, doubled[pos:end]) if a == b)
                    if m >= 0.90 * len(candidate):
                        return pos, end
        return None

    def _has_host_anchor(read_seq: str) -> bool:
        if not host_set or len(read_seq) < 21:
            return False
        u = read_seq.upper()
        for i in range(0, len(u) - 20, 5):
            sub = u[i : i + 21]
            rc = revcomp(sub)
            if (sub if sub <= rc else rc) in host_set:
                return True
        return False

    for _rid, r1, r2 in pairs:
        read_lengths.append(len(r1))
        loc1 = _locate_on_doubled(r1)
        loc2 = _locate_on_doubled(r2) if r2 else None

        for r_seq, loc in ((r1, loc1), (r2, loc2)):
            if r_seq is None:
                continue
            if loc is not None:
                st, en = loc
                total_aligned_bases += (en - st)
                # Seam is at coordinate g_len on `doubled`!
                if st <= (g_len - k_anchor) and en >= (g_len + k_anchor):
                    n_split += 1
                    origin_cov_count += 1
                elif st < g_len < en:
                    origin_cov_count += 1
                if _has_host_anchor(r_seq):
                    n_chimera += 1

        if loc1 is not None and loc2 is not None:
            st1, en1 = loc1
            st2, en2 = loc2
            # Spanning pair: one mate strictly before g_len (in 3' end of monomer), one mate in Revolution 2 (>= g_len)
            if (en1 <= g_len and st2 >= g_len and (en2 - st1) <= 1000) or (
                en2 <= g_len and st1 >= g_len and (en1 - st2) <= 1000
            ):
                n_span += 1
        elif (loc1 is not None and r2 is not None and _has_host_anchor(r2)) or (
            loc2 is not None and _has_host_anchor(r1)
        ):
            n_chimera += 1

    mean_rl = int(round(sum(read_lengths) / float(len(read_lengths)))) if read_lengths else 150
    d_med = total_aligned_bases / float(g_len)
    d_ori = float(origin_cov_count) if origin_cov_count > 0 else (d_med if n_span > 0 else 0.0)
    return JunctionReadEvidence(
        d_median=d_med,
        d_origin=d_ori,
        n_split=n_split,
        n_spanning_pairs=n_span,
        n_host_viral_chimera=n_chimera,
        read_length=mean_rl,
    )
