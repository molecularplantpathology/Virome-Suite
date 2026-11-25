"""Strand-Aware Canonical K-mer Iterative Paired Read Recruitment (`viro_drift.recruitment`)."""

from __future__ import annotations

import dataclasses
from typing import Dict, Iterable, List, Optional, Set, Tuple

from viro_drift.host_partition import canonical_kmer
from viro_drift.translation_bait import shannon_entropy_nt


@dataclasses.dataclass
class RecruitmentResult:
    """Summary and recruited read pairs from iterative recruitment."""

    recruited_pairs: List[Tuple[str, str, Optional[str]]]
    rounds_executed: int
    stop_reason: str
    per_round_counts: List[int]
    filtered_host_decoys_count: int


def _extract_terminal_canonical_kmers(
    seq: str,
    k: int,
    host_mask: Set[str],
    min_entropy: float = 1.30,
    end_window: int = 35,
) -> Set[str]:
    """Extract canonical k-mers from 5' and 3' read ends (`P0-VD-3`)."""
    s = seq.upper()
    if len(s) < k:
        return set()
    out: Set[str] = set()
    w = min(len(s), max(k, end_window))
    regions = [s[:w], s[len(s) - w :]]
    for reg in regions:
        for i in range(0, len(reg) - k + 1):
            sub = reg[i : i + k]
            if "N" in sub:
                continue
            if shannon_entropy_nt(sub) < min_entropy:
                continue
            ck = canonical_kmer(sub)
            if ck not in host_mask:
                out.add(ck)
    return out


def _all_canonical_kmers(seq: str, k: int) -> Set[str]:
    s = seq.upper()
    if len(s) < k:
        return set()
    return {
        canonical_kmer(s[i : i + k])
        for i in range(0, len(s) - k + 1)
        if "N" not in s[i : i + k]
    }


def _estimate_frontier_branching(seqs: Iterable[str], k: int = 21) -> float:
    """Compute mean out-degree on k-mer prefix -> suffix transitions to guard repeat runaway."""
    adj: Dict[str, Set[str]] = {}
    km1 = k - 1
    if km1 < 5:
        return 1.0
    for seq in seqs:
        s = seq.upper()
        for i in range(0, len(s) - k + 1):
            sub = s[i : i + k]
            if "N" not in sub:
                adj.setdefault(sub[:km1], set()).add(sub[1:])
    if not adj:
        return 0.0
    return sum(len(v) for v in adj.values()) / float(len(adj))


def iterative_recruit_reads(
    all_candidate_pairs: Iterable[Tuple[str, str, Optional[str]]],
    initial_seed_ids: Set[str],
    host_kmer_mask: Optional[Set[str]] = None,
    k_nt: int = 21,
    max_rounds: int = 5,
    min_growth_rate: float = 0.01,
    max_branching_degree: float = 4.0,
) -> RecruitmentResult:
    """Iteratively recruit paired mates and strand-aware overlapping reads (`P0-VD-3`)."""
    pairs_map: Dict[str, Tuple[str, str, Optional[str]]] = {
        rid: (rid, r1, r2) for rid, r1, r2 in all_candidate_pairs
    }
    host_mask = host_kmer_mask or set()

    if not initial_seed_ids:
        return RecruitmentResult(
            recruited_pairs=[],
            rounds_executed=0,
            stop_reason="NO_SEED_READS",
            per_round_counts=[],
            filtered_host_decoys_count=0,
        )

    # Filter out any initial seeds that collide with host mask (>20% host k-mers)
    recruited: Set[str] = set()
    filtered_decoys = 0
    for rid in initial_seed_ids:
        if rid not in pairs_map:
            continue
        _, r1, r2 = pairs_map[rid]
        c_kmers = _all_canonical_kmers(r1, k_nt) | (_all_canonical_kmers(r2, k_nt) if r2 else set())
        if host_mask and c_kmers and (len(c_kmers & host_mask) / float(len(c_kmers))) > 0.20:
            filtered_decoys += 1
            continue
        recruited.add(rid)

    if not recruited:
        return RecruitmentResult(
            recruited_pairs=[],
            rounds_executed=0,
            stop_reason="ALL_SEEDS_FILTERED_AS_HOST_DECOYS",
            per_round_counts=[],
            filtered_host_decoys_count=filtered_decoys,
        )

    # Pre-index unrecruited candidate reads by canonical k-mers for fast lookup
    kmer_to_rids: Dict[str, Set[str]] = {}
    read_host_flagged: Set[str] = set()
    counted_decoy_rids: Set[str] = set()
    for rid, (_, r1, r2) in pairs_map.items():
        cks = _all_canonical_kmers(r1, k_nt)
        if r2:
            cks |= _all_canonical_kmers(r2, k_nt)
        if host_mask and cks and (len(cks & host_mask) / float(len(cks))) > 0.20:
            read_host_flagged.add(rid)
        for ck in cks:
            kmer_to_rids.setdefault(ck, set()).add(rid)

    frontier: Set[str] = set(recruited)
    per_round_counts: List[int] = [len(recruited)]
    stop_reason = "MAX_ROUNDS_REACHED"
    rounds_done = 0

    for _round_idx in range(1, max_rounds + 1):
        rounds_done = _round_idx
        frontier_seqs: List[str] = []
        bait_kmers: Set[str] = set()
        for rid in frontier:
            _, r1, r2 = pairs_map[rid]
            frontier_seqs.append(r1)
            bait_kmers |= _extract_terminal_canonical_kmers(r1, k_nt, host_mask)
            if r2:
                frontier_seqs.append(r2)
                bait_kmers |= _extract_terminal_canonical_kmers(r2, k_nt, host_mask)

        if _estimate_frontier_branching(frontier_seqs, k=min(k_nt, 17)) > max_branching_degree:
            stop_reason = "BRANCHING_CEILING_EXCEEDED"
            break

        newly_added: Set[str] = set()
        for bk in bait_kmers:
            for cand_rid in kmer_to_rids.get(bk, ()):
                if cand_rid in recruited:
                    continue
                if cand_rid in read_host_flagged:
                    if cand_rid not in counted_decoy_rids:
                        counted_decoy_rids.add(cand_rid)
                        filtered_decoys += 1
                    continue
                newly_added.add(cand_rid)

        prev_count = len(recruited)
        recruited |= newly_added
        per_round_counts.append(len(recruited))
        growth = len(newly_added) / float(max(1, prev_count))

        if growth < min_growth_rate:
            stop_reason = "GROWTH_PLATEAU_REACHED"
            break
        frontier = newly_added

    ordered_pairs = [pairs_map[rid] for rid in sorted(recruited)]
    return RecruitmentResult(
        recruited_pairs=ordered_pairs,
        rounds_executed=rounds_done,
        stop_reason=stop_reason,
        per_round_counts=per_round_counts,
        filtered_host_decoys_count=filtered_decoys,
    )
