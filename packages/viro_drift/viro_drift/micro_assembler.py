"""Two-Phase Cycle-Safe Local De Bruijn Micro-Assembler (`viro_drift.micro_assembler`)."""

from __future__ import annotations

import collections
import dataclasses
from typing import Dict, Iterable, List, Optional, Set, Tuple

from viro_drift.host_partition import revcomp


@dataclasses.dataclass
class AssembledContig:
    """Locally assembled viral candidate contig."""

    contig_id: str
    sequence: str
    mean_kmer_depth: float
    is_circular_cycle: bool


def _choose_kmer_cascade(reads: List[str]) -> Tuple[int, ...]:
    if not reads:
        return (21, 33, 55)
    min_len = min(len(r) for r in reads)
    if min_len < 60:
        return tuple(k for k in (15, 21, 31) if k < min_len) or (11,)
    return tuple(k for k in (21, 33, 55) if k < min_len) or (21,)


def _build_dbg_and_extract_unitigs(
    reads: List[str],
    k: int,
    min_edge_cov: int = 1,
) -> List[Tuple[str, float, bool]]:
    """Build directed de Bruijn graph and run Two-Phase Cycle-Safe Traversal (`P0-VD-3`).

    Phase 1: Walk all non-branching paths starting from nodes with in_deg != 1 or out_deg != 1.
    Phase 2: Walk remaining unvisited closed circular cycles starting from lexicographically
             smallest k-mer node, traversing each cycle edge once and breaking cleanly.
    """
    edge_counts: Dict[Tuple[str, str], int] = collections.Counter()
    for r in reads:
        s = r.upper()
        if len(s) < k + 1:
            continue
        for i in range(0, len(s) - k):
            sub = s[i : i + k + 1]
            if "N" not in sub:
                u, v = sub[:k], sub[1:]
                edge_counts[(u, v)] += 1

    # Filter ultra-low singleton sequencing error tips if graph has deep coverage
    max_cov = max(edge_counts.values(), default=0)
    eff_min_cov = max(min_edge_cov, 2 if max_cov >= 8 else 1)
    active_edges = {e: c for e, c in edge_counts.items() if c >= eff_min_cov}
    if not active_edges and edge_counts:
        active_edges = dict(edge_counts)

    out_adj: Dict[str, List[str]] = collections.defaultdict(list)
    in_deg: Dict[str, int] = collections.Counter()
    nodes: Set[str] = set()

    for (u, v), _cov in sorted(active_edges.items()):
        out_adj[u].append(v)
        in_deg[v] += 1
        nodes.add(u)
        nodes.add(v)

    # Sort outgoing edges by descending coverage for greedy primary path selection at bubbles
    for u in out_adj:
        out_adj[u].sort(key=lambda v: (-active_edges[(u, v)], v))

    # Prune short dead-end tips (< 2*k bp) branching off higher-coverage nodes
    for u in list(out_adj.keys()):
        if len(out_adj[u]) > 1:
            best_v = out_adj[u][0]
            best_cov = active_edges[(u, best_v)]
            kept = [best_v]
            for alt_v in out_adj[u][1:]:
                if active_edges[(u, alt_v)] * 4 >= best_cov or len(out_adj.get(alt_v, [])) > 0:
                    kept.append(alt_v)
                else:
                    in_deg[alt_v] -= 1
                    active_edges.pop((u, alt_v), None)
            out_adj[u] = kept

    visited_edges: Set[Tuple[str, str]] = set()
    contigs: List[Tuple[str, float, bool]] = []

    # Phase 1: Start walks from all non-1-in-1-out nodes
    start_nodes = sorted(
        n for n in nodes if (in_deg[n] != 1 or len(out_adj.get(n, [])) != 1) and len(out_adj.get(n, [])) > 0
    )
    for start in start_nodes:
        for nxt in out_adj[start]:
            if (start, nxt) in visited_edges:
                continue
            seq_chars = [start, nxt[-1]]
            covs = [active_edges[(start, nxt)]]
            visited_edges.add((start, nxt))
            curr = nxt
            while in_deg[curr] == 1 and len(out_adj.get(curr, [])) == 1:
                w = out_adj[curr][0]
                if (curr, w) in visited_edges:
                    break
                visited_edges.add((curr, w))
                seq_chars.append(w[-1])
                covs.append(active_edges[(curr, w)])
                curr = w
            contigs.append(("".join(seq_chars), sum(covs) / float(len(covs)), False))

    # Phase 2: Pure closed circular cycles (`in_deg == 1` and `out_deg == 1` everywhere)
    remaining_edges = sorted(e for e in active_edges if e not in visited_edges)
    for u0, _v0 in remaining_edges:
        if (u0, out_adj[u0][0]) in visited_edges:
            continue
        # Find lexicographically smallest node in this cycle for deterministic canonical start
        cycle_nodes = [u0]
        curr = out_adj[u0][0]
        safety = 0
        while curr != u0 and curr not in cycle_nodes and len(out_adj.get(curr, [])) > 0 and safety <= len(active_edges):
            cycle_nodes.append(curr)
            curr = out_adj[curr][0]
            safety += 1
        start = min(cycle_nodes)
        curr = start
        seq_chars = [start]
        covs = []
        is_closed = False
        while len(out_adj.get(curr, [])) > 0:
            nxt = out_adj[curr][0]
            if (curr, nxt) in visited_edges:
                is_closed = (nxt == start or curr == start)
                break
            visited_edges.add((curr, nxt))
            covs.append(active_edges[(curr, nxt)])
            if nxt == start:
                # Returning to start completes exact 1.0x circular k-mer loop (plus k-1 overlap if needed)
                is_closed = True
                break
            seq_chars.append(nxt[-1])
            curr = nxt
        if covs:
            contigs.append(("".join(seq_chars), sum(covs) / float(len(covs)), is_closed))

    return contigs


def _deduplicate_rc_contigs(
    raw_contigs: List[Tuple[str, float, bool]],
    min_len: int,
    k: int = 21,
) -> List[Tuple[str, float, bool]]:
    """Collapse reverse-complement twin contigs (including rotated circular cycles)."""
    filtered = [c for c in raw_contigs if len(c[0]) >= min_len]
    filtered.sort(key=lambda x: (-len(x[0]), -x[1], x[0]))
    kept: List[Tuple[str, float, bool]] = []
    seen_linear: List[str] = []
    seen_circular_doubled: List[str] = []
    for seq, cov, is_circ in filtered:
        rc = revcomp(seq)
        if any(seq in s or rc in s for s in seen_linear):
            continue
        if is_circ and len(seq) > k:
            core = seq[: -(k - 1)] if seq[: k - 1] == seq[-(k - 1) :] else seq
            rc_core = revcomp(core)
            if any(core in dc or rc_core in dc for dc in seen_circular_doubled):
                continue
            seen_circular_doubled.append(core + core)
        kept.append((seq, cov, is_circ))
        seen_linear.append(seq)
    return kept


def assemble_viral_contigs(
    recruited_pairs: Iterable[Tuple[str, str, Optional[str]]],
    k_values: Optional[Tuple[int, ...]] = None,
    min_contig_len: int = 120,
) -> List[AssembledContig]:
    """Assemble isolated viral reads across an adaptive odd-k cascade and select optimal contiguity."""
    strand_expanded_reads: List[str] = []
    raw_reads: List[str] = []
    for _rid, r1, r2 in recruited_pairs:
        raw_reads.append(r1)
        strand_expanded_reads.append(r1.upper())
        strand_expanded_reads.append(revcomp(r1.upper()))
        if r2:
            raw_reads.append(r2)
            strand_expanded_reads.append(r2.upper())
            strand_expanded_reads.append(revcomp(r2.upper()))

    if not raw_reads:
        return []

    cascade = k_values or _choose_kmer_cascade(raw_reads)
    best_unitigs: List[Tuple[str, float, bool]] = []
    best_max_len = -1

    for k in cascade:
        if k % 2 == 0:
            k += 1  # Enforce odd k to prevent palindromic self-complement k-mers (`P0-VD-3`)
        unitigs = _build_dbg_and_extract_unitigs(strand_expanded_reads, k=k)
        dedup = _deduplicate_rc_contigs(unitigs, min_len=min(min_contig_len, max(len(r) for r in raw_reads)), k=k)
        top_len = max((len(u[0]) for u in dedup), default=0)
        if top_len > best_max_len:
            best_max_len = top_len
            best_unitigs = dedup

    results: List[AssembledContig] = []
    for idx, (seq, cov, is_circ) in enumerate(best_unitigs, start=1):
        results.append(
            AssembledContig(
                contig_id=f"viro_drift_contig_{idx:03d}",
                sequence=seq,
                mean_kmer_depth=round(cov / 2.0, 2),
                is_circular_cycle=is_circ,
            )
        )
    return results
