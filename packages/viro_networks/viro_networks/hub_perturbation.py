"""Zero-`networkx` Pure-Python Centrality & Hub Perturbation Engine (`viro_networks.hub_perturbation`)."""

from __future__ import annotations

import collections
import dataclasses
import math
from typing import Dict, List, Sequence, Set, Tuple

import numpy as np


@dataclasses.dataclass
class HubVulnerabilityRecord:
    """Topological centrality and node-knockout structural impact (non-causal)."""

    taxon: str
    degree: int
    degree_centrality: float
    betweenness_centrality: float
    natural_connectivity_drop_pct: float
    lcc_fraction_after_removal: float


def _natural_connectivity(adj: np.ndarray) -> float:
    """Compute Estrada/Wu natural connectivity ln((1/N) * sum(exp(lambda_i)))."""
    n = adj.shape[0]
    if n == 0:
        return 0.0
    eigvals = np.linalg.eigvalsh(adj)
    max_ev = float(np.max(eigvals))
    # Log-sum-exp for overflow safety
    lse = max_ev + math.log(float(np.mean(np.exp(eigvals - max_ev))))
    return lse


def _largest_connected_component_size(adj: np.ndarray) -> int:
    n = adj.shape[0]
    if n == 0:
        return 0
    visited: Set[int] = set()
    best = 0
    for start in range(n):
        if start in visited:
            continue
        q = collections.deque([start])
        visited.add(start)
        comp = 0
        while q:
            u = q.popleft()
            comp += 1
            for v in range(n):
                if adj[u, v] > 0 and v not in visited:
                    visited.add(v)
                    q.append(v)
        if comp > best:
            best = comp
    return best


def _brandes_betweenness(adj: np.ndarray) -> List[float]:
    """Unweighted Brandes (2001) shortest-path betweenness centrality in pure Python (`P1-VN-3`)."""
    n = adj.shape[0]
    cb = [0.0] * n
    if n <= 2:
        return cb

    neighbors = [[v for v in range(n) if adj[u, v] > 0] for u in range(n)]

    for s in range(n):
        stack: List[int] = []
        preds: List[List[int]] = [[] for _ in range(n)]
        sigma = [0.0] * n
        sigma[s] = 1.0
        dist = [-1] * n
        dist[s] = 0
        q = collections.deque([s])

        while q:
            v = q.popleft()
            stack.append(v)
            for w in neighbors[v]:
                if dist[w] < 0:
                    q.append(w)
                    dist[w] = dist[v] + 1
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    preds[w].append(v)

        delta = [0.0] * n
        while stack:
            w = stack.pop()
            for v in preds[w]:
                if sigma[w] > 0:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                cb[w] += delta[w]

    # Normalize for undirected graph: divide by (n - 1)(n - 2)
    scale = 1.0 / float((n - 1) * (n - 2))
    return [val * scale for val in cb]


def evaluate_hub_perturbation(
    taxa_names: Sequence[str],
    edge_pairs: Sequence[Tuple[str, str]],
) -> List[HubVulnerabilityRecord]:
    """Rank candidate topological hubs via centrality and in silico node deletion (`P1-VN-3`)."""
    names = list(taxa_names)
    n = len(names)
    if n == 0:
        return []

    idx_map: Dict[str, int] = {t: i for i, t in enumerate(names)}
    adj = np.zeros((n, n), dtype=float)
    for a, b in edge_pairs:
        if a in idx_map and b in idx_map and a != b:
            i, j = idx_map[a], idx_map[b]
            adj[i, j] = 1.0
            adj[j, i] = 1.0

    base_nc = _natural_connectivity(adj)
    btw = _brandes_betweenness(adj)
    records: List[HubVulnerabilityRecord] = []

    for i, taxon in enumerate(names):
        deg = int(np.sum(adj[i, :] > 0))
        deg_c = float(deg) / float(max(1, n - 1))

        if n > 1:
            keep = [k for k in range(n) if k != i]
            sub_adj = adj[np.ix_(keep, keep)]
            rem_nc = _natural_connectivity(sub_adj)
            drop_pct = (
                max(0.0, (base_nc - rem_nc) / base_nc * 100.0)
                if base_nc > 1e-9
                else 0.0
            )
            lcc_frac = float(_largest_connected_component_size(sub_adj)) / float(n - 1)
        else:
            drop_pct = 0.0
            lcc_frac = 0.0

        records.append(
            HubVulnerabilityRecord(
                taxon=taxon,
                degree=deg,
                degree_centrality=deg_c,
                betweenness_centrality=btw[i],
                natural_connectivity_drop_pct=drop_pct,
                lcc_fraction_after_removal=lcc_frac,
            )
        )

    records.sort(key=lambda r: (-r.degree, -r.natural_connectivity_drop_pct, r.taxon))
    return records
