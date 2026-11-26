"""TIR Veto & Indel-Tolerant Revolution-2 Concatemer Monomer Curator (`viral_genome_tools.monomer_curator`)."""

from __future__ import annotations

import dataclasses
from typing import List, Set

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclasses.dataclass
class MonomerCurationResult:
    """Output of Stage 1 concatemer trimming and TIR evaluation."""

    monomer_seq: str
    original_length: int
    monomer_length: int
    was_trimmed: bool
    has_terminal_inverted_repeat: bool
    overlap_length: int
    overlap_identity: float


def detect_terminal_inverted_repeat(
    seq: str,
    min_arm: int = 50,
    max_arm: int = 500,
    min_identity: float = 0.90,
) -> bool:
    """Precedence Rule 0: Detect Terminal Inverted Repeats (TIRs) to veto false circularization."""
    s = seq.upper()
    n = len(s)
    upper_w = min(max_arm, n // 3)
    if upper_w < min_arm:
        return False

    # Compare 5' prefix S[0:w] against reverse-complement of 3' suffix rc(S[n-w:n])
    rc_full = revcomp(s)
    # Note: rc(S[n-w:n]) == rc_full[0:w]!
    # If S[0:w] matches rc_full[0:w] from index 0 across >= min_arm bases with >= 90% identity:
    best_tir_id = 0.0
    best_dtr_id = 0.0
    for w in range(min_arm, upper_w + 1, 10):
        matches_tir = sum(1 for a, b in zip(s[:w], rc_full[:w]) if a == b)
        tir_id = matches_tir / float(w)
        if tir_id > best_tir_id:
            best_tir_id = tir_id

        matches_dtr = sum(1 for a, b in zip(s[:w], s[n - w :]) if a == b)
        dtr_id = matches_dtr / float(w)
        if dtr_id > best_dtr_id:
            best_dtr_id = dtr_id

    # Also check exact sliding window near 0 in case of 1-2 bp terminal adapter trim
    if best_tir_id < min_identity and upper_w >= min_arm:
        w = min_arm
        matches_tir = sum(1 for a, b in zip(s[:w], rc_full[:w]) if a == b)
        best_tir_id = max(best_tir_id, matches_tir / float(w))

    return best_tir_id >= min_identity and best_tir_id > best_dtr_id


def _banded_prefix_suffix_identity(
    prefix: str,
    suffix: str,
    band: int = 16,
) -> float:
    """Banded semi-global Levenshtein/alignment identity anchored at (0, 0) through end of suffix (`P0-VGT-1`).

    Computes minimum edit distance between `suffix` (length m) and `prefix[:m_aligned]`
    allowing free end on `prefix` if `len(prefix) > len(suffix)` (or vice-versa atcontig end),
    within diagonal band `[-band, +band]` so small indels inside the overlap tail never shift P*.
    """
    m = len(suffix)
    n = len(prefix)
    if m == 0 or n == 0:
        return 0.0

    # Fast exact check first (O(m))
    cmp_len = min(m, n)
    if abs(m - n) == 0 or n >= m:
        exact_matches = sum(1 for a, b in zip(prefix[:cmp_len], suffix[:cmp_len]) if a == b)
        if exact_matches / float(cmp_len) >= 0.96 and cmp_len == m:
            return exact_matches / float(cmp_len)

    # Banded Myers/Needleman-Wunsch edit distance anchored at (0, 0)
    inf = m + n + 10
    # dp[j] = min edits to align suffix[:i] with prefix[:j] for j in [max(0, i-band)..min(n, i+band)]
    prev_lo = 0
    prev_hi = min(n, band)
    prev_dp = {j: j for j in range(prev_lo, prev_hi + 1)}

    for i in range(1, m + 1):
        lo = max(0, i - band)
        hi = min(n, i + band)
        curr_dp = {}
        s_ch = suffix[i - 1]
        for j in range(lo, hi + 1):
            # Deletion from suffix (vertical: i-1, j)
            c_del = prev_dp.get(j, inf) + 1
            # Insertion into suffix (horizontal: i, j-1)
            c_ins = curr_dp.get(j - 1, inf) + 1
            # Match/substitution (diagonal: i-1, j-1)
            if j >= 1:
                cost_sub = 0 if prefix[j - 1] == s_ch else 1
                c_diag = prev_dp.get(j - 1, inf) + cost_sub
            else:
                c_diag = inf
            curr_dp[j] = min(c_del, c_ins, c_diag)
        prev_dp = curr_dp

    # Best alignment ending at end of suffix (row m) within band around m (or at end of prefix n if n < m)
    if n >= m:
        candidates = [prev_dp[j] for j in prev_dp if abs(j - m) <= band]
        denom = float(m)
    else:
        candidates = [prev_dp[n]] if n in prev_dp else list(prev_dp.values())
        denom = float(max(m, n))

    if not candidates:
        return 0.0
    min_edits = min(candidates)
    return max(0.0, 1.0 - (float(min_edits) / denom))


def trim_concatemer_to_monomer(
    seq: str,
    min_monomer_len: int = 500,
    min_overlap: int = 40,
    min_identity: float = 0.95,
    band: int = 16,
) -> MonomerCurationResult:
    """Resolve 1.05x-3.0x circular concatemers to exact 1.0x monomer via Revolution-2 anchor (`P0-VGT-1`)."""
    s = seq.upper().strip()
    n = len(s)

    if detect_terminal_inverted_repeat(s):
        return MonomerCurationResult(
            monomer_seq=s,
            original_length=n,
            monomer_length=n,
            was_trimmed=False,
            has_terminal_inverted_repeat=True,
            overlap_length=0,
            overlap_identity=0.0,
        )

    eff_min_monomer = min(min_monomer_len, max(100, n // 4))
    if n < eff_min_monomer + min_overlap:
        return MonomerCurationResult(
            monomer_seq=s,
            original_length=n,
            monomer_length=n,
            was_trimmed=False,
            has_terminal_inverted_repeat=False,
            overlap_length=0,
            overlap_identity=0.0,
        )

    # Seed k-mers from 5' prefix S[0 : probe_len] to propose Revolution-2 start coordinate p2 (`P0-VGT-1`)
    probe_len = min(250, n // 2)
    candidate_periods: Set[int] = set()
    for seed_k in (15, 11):
        if probe_len < seed_k:
            continue
        seed_offsets = {}
        for s_off in range(0, probe_len - seed_k + 1, 3):
            kmer = s[s_off : s_off + seed_k]
            if "N" not in kmer and kmer not in seed_offsets:
                seed_offsets[kmer] = s_off

        search_start = eff_min_monomer
        for pos in range(search_start, n - seed_k + 1):
            kmer = s[pos : pos + seed_k]
            if kmer in seed_offsets:
                p2 = pos - seed_offsets[kmer]
                if eff_min_monomer <= p2 <= (n - min_overlap):
                    candidate_periods.add(p2)
        if candidate_periods:
            break

    # Group nearby candidate periods within 2*band (indel-shifted hits of the same revolution)
    # so we always pick the fundamental 1.0x cluster first, and within that cluster select the
    # exact Revolution-2 start p2 with strongest opening anchor + full suffix identity (`P0-VGT-1`).
    sorted_periods = sorted(candidate_periods)
    clusters: List[List[int]] = []
    for p2 in sorted_periods:
        if not clusters or (p2 - clusters[-1][-1]) > (2 * band):
            clusters.append([p2])
        else:
            clusters[-1].append(p2)

    for cluster in clusters:
        best_choice = None
        best_key = (-1.0, -1.0)
        for p2 in cluster:
            suffix = s[p2:]
            prefix = s[: min(n, len(suffix) + band)]
            ident = _banded_prefix_suffix_identity(prefix=prefix, suffix=suffix, band=band)
            if ident < min_identity:
                continue
            if n >= 2 * p2:
                rev_ident = _banded_prefix_suffix_identity(s[:p2], s[p2 : 2 * p2], band=band)
                if rev_ident < min_identity:
                    continue
            anchor_w = min(30, len(suffix), p2)
            anchor_exact = sum(1 for a, b in zip(s[:anchor_w], suffix[:anchor_w]) if a == b) / float(max(1, anchor_w))
            score_key = (round(ident, 6), round(anchor_exact, 6))
            if score_key > best_key:
                best_key = score_key
                best_choice = (p2, ident)

        if best_choice is not None:
            p2_win, ident_win = best_choice
            monomer = s[:p2_win]
            return MonomerCurationResult(
                monomer_seq=monomer,
                original_length=n,
                monomer_length=p2_win,
                was_trimmed=True,
                has_terminal_inverted_repeat=False,
                overlap_length=n - p2_win,
                overlap_identity=round(ident_win, 4),
            )

    return MonomerCurationResult(
        monomer_seq=s,
        original_length=n,
        monomer_length=n,
        was_trimmed=False,
        has_terminal_inverted_repeat=False,
        overlap_length=0,
        overlap_identity=0.0,
    )
