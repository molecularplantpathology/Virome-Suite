"""Intra-Host Viral Population Genetics Metrics (`quasispecies_profiler.diversity_metrics`)."""

from __future__ import annotations

import math
from typing import Dict, Optional, Sequence, Tuple

_STANDARD_CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def compute_shannon_entropy(vafs: Sequence[float], genome_length: Optional[int] = None) -> float:
    """Compute mean normalized per-site Shannon entropy H in bits across genome length L."""
    l_denom = int(genome_length) if (genome_length is not None and genome_length > 0) else max(1, len(vafs))
    total_h = 0.0
    ln2 = math.log(2.0)
    for mu in vafs:
        p = min(max(float(mu), 0.0), 1.0)
        if 0.0 < p < 1.0:
            total_h += -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p)) / ln2
    return total_h / float(l_denom)


def compute_nucleotide_diversity_pi(
    vafs: Sequence[float],
    depths: Sequence[float],
    genome_length: Optional[int] = None,
) -> float:
    """Compute unbiased intra-host nucleotide diversity pi with n_i < 2 zero-division guard (`P1-QP-4`)."""
    l_denom = int(genome_length) if (genome_length is not None and genome_length > 0) else max(1, len(vafs))
    total_pi = 0.0
    for mu, n_val in zip(vafs, depths):
        n_f = float(n_val)
        if n_f < 2.0:
            continue
        p = min(max(float(mu), 0.0), 1.0)
        total_pi += (n_f / (n_f - 1.0)) * 2.0 * p * (1.0 - p)
    return total_pi / float(l_denom)


def _codon_syn_nonsyn_sites(codon: str) -> Tuple[float, float]:
    """Count Nei-Gojobori (1986) synonymous (s) and non-synonymous (n) sites for a sense codon (`P1-QP-4`)."""
    c = codon.upper()
    aa = _STANDARD_CODON_TABLE.get(c, "*")
    if aa == "*" or len(c) != 3:
        return 0.0, 0.0

    s_sites = 0.0
    n_sites = 0.0
    for pos in range(3):
        syn_neighbors = 0
        valid_neighbors = 0
        for base in "ACGT":
            if base == c[pos]:
                continue
            mut = c[:pos] + base + c[pos + 1 :]
            mut_aa = _STANDARD_CODON_TABLE.get(mut, "*")
            if mut_aa == "*":
                continue  # Exclude nonsense stop-codon pathways (`P1-QP-4`)
            valid_neighbors += 1
            if mut_aa == aa:
                syn_neighbors += 1
        if valid_neighbors > 0:
            frac_s = syn_neighbors / float(valid_neighbors)
            s_sites += frac_s
            n_sites += (1.0 - frac_s)
    return s_sites, n_sites


def _jukes_cantor_distance(p: float) -> float:
    """Guarded Jukes-Cantor saturation correction (`P1-QP-4`)."""
    if p <= 0.0:
        return 0.0
    if p >= 0.749999:
        return -0.75 * math.log(1.0 - (4.0 / 3.0) * 0.749999)
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)


def compute_nei_gojobori_dnds(
    cds_seq: str,
    coding_variants: Sequence[Tuple[int, str, float]],
) -> Dict[str, Optional[float]]:
    """Compute frequency-weighted Nei-Gojobori dN, dS, and dN/dS across a coding sequence (`P1-QP-4`).

    `coding_variants`: sequence of `(zero_based_cds_nt_offset, alt_base, vaf)`.
    """
    seq = cds_seq.upper()
    n_codons = len(seq) // 3
    total_s_sites = 0.0
    total_n_sites = 0.0

    for c_idx in range(n_codons):
        codon = seq[3 * c_idx : 3 * c_idx + 3]
        s_c, n_c = _codon_syn_nonsyn_sites(codon)
        total_s_sites += s_c
        total_n_sites += n_c

    s_diffs = 0.0
    n_diffs = 0.0
    for nt_offset, alt_base, vaf in coding_variants:
        c_idx = int(nt_offset) // 3
        pos_in_codon = int(nt_offset) % 3
        if c_idx < 0 or c_idx >= n_codons:
            continue
        ref_codon = seq[3 * c_idx : 3 * c_idx + 3]
        ref_aa = _STANDARD_CODON_TABLE.get(ref_codon, "*")
        if ref_aa == "*":
            continue
        mut_codon = ref_codon[:pos_in_codon] + alt_base.upper() + ref_codon[pos_in_codon + 1 :]
        mut_aa = _STANDARD_CODON_TABLE.get(mut_codon, "*")
        if mut_aa == "*":
            continue
        w = min(max(float(vaf), 0.0), 1.0)
        if mut_aa == ref_aa:
            s_diffs += w
        else:
            n_diffs += w

    p_s = s_diffs / total_s_sites if total_s_sites > 0 else 0.0
    p_n = n_diffs / total_n_sites if total_n_sites > 0 else 0.0
    d_s = _jukes_cantor_distance(p_s)
    d_n = _jukes_cantor_distance(p_n)
    dnds: Optional[float] = (d_n / d_s) if d_s > 0.0 else None

    return {
        "S_sites": total_s_sites,
        "N_sites": total_n_sites,
        "p_S": p_s,
        "p_N": p_n,
        "d_S": d_s,
        "d_N": d_n,
        "dN_dS": dnds,
    }
