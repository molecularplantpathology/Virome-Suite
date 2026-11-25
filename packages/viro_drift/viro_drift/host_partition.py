"""Evidence-Preserving Host Read Partitioning (`viro_drift.host_partition`)."""

from __future__ import annotations

import dataclasses
from typing import Dict, Iterable, List, Optional, Set, Tuple

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def canonical_kmer(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def build_host_kmer_index(host_sequences: Iterable[str], k: int = 21) -> Set[str]:
    """Index canonical nucleotide k-mers from host reference sequences."""
    index: Set[str] = set()
    for seq in host_sequences:
        s = seq.upper()
        n = len(s)
        for i in range(0, n - k + 1):
            sub = s[i : i + k]
            if "N" not in sub:
                index.add(canonical_kmer(sub))
    return index


def _host_kmer_hit_profile(seq: str, host_kmers: Set[str], k: int = 21) -> Tuple[float, bool, bool]:
    """Return (fraction_host_kmers, has_host_segment_20bp, has_nonhost_segment_20bp)."""
    s = seq.upper()
    if len(s) < k or not host_kmers:
        return 0.0, False, len(s) >= 20

    total = len(s) - k + 1
    hits = [
        1 if ("N" not in s[i : i + k] and canonical_kmer(s[i : i + k]) in host_kmers) else 0
        for i in range(total)
    ]
    frac = sum(hits) / float(total)

    # Check if read is chimeric: contains at least one host k-mer (>=21bp) AND at least 20bp non-host stretch
    has_host = any(h == 1 for h in hits)
    # Count contiguous zero run of length >= 20 - k + 1 (Wait: if hits[i]==0 for m consecutive positions,
    # the uncovered span is m + k - 1 bp; even 2 consecutive non-host 21-mers without overlap means >=22bp non-host)
    zero_run = 0
    max_zero_run = 0
    for h in hits:
        if h == 0:
            zero_run += 1
            if zero_run > max_zero_run:
                max_zero_run = zero_run
        else:
            zero_run = 0
    has_nonhost_20bp = max_zero_run >= k  # >= k consecutive non-host k-mer starts guarantees >= 21bp non-host free of host k-mer
    return frac, has_host, has_nonhost_20bp


@dataclasses.dataclass
class PartitionedReads:
    """Three-stream evidence-preserving read partition output."""

    host_concordant: List[Tuple[str, str, Optional[str]]]
    host_spanning_ambiguous: List[Tuple[str, str, Optional[str]]]
    viral_candidate: List[Tuple[str, str, Optional[str]]]

    def summary_counts(self) -> Dict[str, int]:
        return {
            "host_concordant_pairs": len(self.host_concordant),
            "host_spanning_ambiguous_pairs": len(self.host_spanning_ambiguous),
            "viral_candidate_pairs": len(self.viral_candidate),
        }


def partition_host_reads(
    read_pairs: Iterable[Tuple[str, str, Optional[str]]],
    host_kmers: Set[str],
    k: int = 21,
    concordant_threshold: float = 0.70,
) -> PartitionedReads:
    """Partition paired or single-end reads into 3 mutually exclusive evidence streams.

    1. host_concordant: Both mates align to host (>= concordant_threshold k-mer identity) without >=21bp non-host tails.
    2. host_spanning_ambiguous: Discordant mate pairs (one host, one non-host) OR split chimeric reads.
    3. viral_candidate: Unmapped non-host reads forwarded to 6-frame viral baiting.
    """
    concordant: List[Tuple[str, str, Optional[str]]] = []
    ambiguous: List[Tuple[str, str, Optional[str]]] = []
    viral: List[Tuple[str, str, Optional[str]]] = []

    for read_id, r1, r2 in read_pairs:
        if not host_kmers:
            viral.append((read_id, r1, r2))
            continue

        f1, h1, nh1 = _host_kmer_hit_profile(r1, host_kmers, k=k)
        if r2 is not None:
            f2, h2, nh2 = _host_kmer_hit_profile(r2, host_kmers, k=k)
            r1_pure_host = f1 >= concordant_threshold and not (h1 and nh1)
            r2_pure_host = f2 >= concordant_threshold and not (h2 and nh2)
            r1_pure_nonhost = not h1
            r2_pure_nonhost = not h2

            if r1_pure_host and r2_pure_host:
                concordant.append((read_id, r1, r2))
            elif r1_pure_nonhost and r2_pure_nonhost:
                viral.append((read_id, r1, r2))
            else:
                ambiguous.append((read_id, r1, r2))
        else:
            if f1 >= concordant_threshold and not (h1 and nh1):
                concordant.append((read_id, r1, None))
            elif not h1:
                viral.append((read_id, r1, None))
            else:
                ambiguous.append((read_id, r1, None))

    return PartitionedReads(
        host_concordant=concordant,
        host_spanning_ambiguous=ambiguous,
        viral_candidate=viral,
    )
