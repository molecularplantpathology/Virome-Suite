"""2x Circular Seam Wraparound ORF Annotator & Type-Aware CDS Validator (`viralbio_reporting.circular_annotator`)."""

from __future__ import annotations

import dataclasses
from typing import Dict, List, Tuple

from Bio.Seq import Seq

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
_STOP_CODONS = {"TAA", "TAG", "TGA"}


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclasses.dataclass
class AnnotatedCDS:
    """Biologically validated complete or partial viral coding feature."""

    feature_id: str
    strand: str  # "+" or "-"
    intervals_1based: List[Tuple[int, int]]  # 1-based inclusive NCBI .tbl intervals
    partial_5prime: bool
    partial_3prime: bool
    codon_start: int  # 1, 2, or 3
    spans_origin: bool
    aa_sequence: str
    nt_length: int
    product: str = "hypothetical viral protein"


def _project_circular_intervals(
    s0: int,
    e0: int,
    genome_len: int,
    strand: str,
) -> List[Tuple[int, int]]:
    """Project 0-based [s0, e0) on S_2x (or rc(S_2x)) into 1-based inclusive NCBI .tbl intervals (`P0-VBR-2`)."""
    l_nt = genome_len
    if strand == "+":
        if e0 <= l_nt:
            return [(s0 + 1, e0)]
        return [(s0 + 1, l_nt), (1, e0 - l_nt)]

    # Negative strand (`-`) on rc(S + S) = rc(S) + rc(S) (`P0-VBR-2`)
    if e0 <= l_nt:
        return [(l_nt - s0, l_nt - e0 + 1)]
    return [(l_nt - s0, 1), (l_nt, 2 * l_nt - e0 + 1)]


def reconstruct_cds_nt_from_intervals(
    genome_seq: str,
    feature: AnnotatedCDS,
) -> str:
    """Independently splice nucleotide sequence from parent genome using 1-based .tbl intervals & codon_start."""
    s = genome_seq.upper()
    parts: List[str] = []
    for st_1, en_1 in feature.intervals_1based:
        if feature.strand == "+":
            parts.append(s[st_1 - 1 : en_1])
        else:
            # Descending 1-based inclusive interval [st_1, en_1] where st_1 >= en_1
            lo = min(st_1, en_1)
            hi = max(st_1, en_1)
            parts.append(revcomp(s[lo - 1 : hi]))
    spliced = "".join(parts)
    offset = max(0, int(feature.codon_start) - 1)
    return spliced[offset:]


def validate_cds_feature(genome_seq: str, feature: AnnotatedCDS) -> Tuple[bool, str]:
    """Type-aware biological invariant verification (coordinates, frame, start/stop, exact translation)."""
    g_len = len(genome_seq)
    if not feature.intervals_1based:
        return False, "Empty interval list"

    for st_1, en_1 in feature.intervals_1based:
        if st_1 < 1 or st_1 > g_len or en_1 < 1 or en_1 > g_len:
            return False, f"Interval ({st_1}, {en_1}) out of genome bounds [1, {g_len}]"
        if feature.strand == "+" and st_1 > en_1:
            return False, f"Forward interval ({st_1}, {en_1}) must be ascending"
        if feature.strand == "-" and st_1 < en_1:
            return False, f"Reverse interval ({st_1}, {en_1}) must be descending"

    coding_nt = reconstruct_cds_nt_from_intervals(genome_seq, feature)
    trim_len = (len(coding_nt) // 3) * 3
    if trim_len < 3:
        return False, "Coding sequence shorter than 1 codon"

    if not feature.partial_5prime and coding_nt[:3] not in ("ATG", "CTG", "GTG"):
        return False, f"Complete CDS lacks valid start codon (found {coding_nt[:3]})"

    if not feature.partial_3prime:
        if len(coding_nt) % 3 != 0:
            return False, "Complete CDS nucleotide length is not a multiple of 3"
        if coding_nt[-3:] not in _STOP_CODONS:
            return False, f"Complete CDS lacks terminal stop codon (found {coding_nt[-3:]})"
        aa_full = str(Seq(coding_nt).translate(table=1))
        if not aa_full.endswith("*"):
            return False, "Translated complete CDS does not end with '*'"
        aa_body = aa_full[:-1]
    else:
        aa_full = str(Seq(coding_nt[:trim_len]).translate(table=1))
        aa_body = aa_full.rstrip("*")

    if "*" in aa_body:
        return False, "Premature internal stop codon detected"

    if aa_body != feature.aa_sequence:
        return False, "Reconstructed translation does not match annotated aa_sequence"

    return True, "VALID"


def annotate_genome_orfs(
    genome_seq: str,
    topology: str = "circular",
    min_aa: int = 60,
    include_partial_linear: bool = True,
) -> List[AnnotatedCDS]:
    """Discover ORFs with 360-degree rotation invariance on circular genomes (`P0-VBR-1`, `P0-VBR-2`, `P1-VBR-3`)."""
    s = genome_seq.upper().strip()
    l_nt = len(s)
    if l_nt < 3:
        return []

    is_circular = (topology.lower() == "circular")
    rc = revcomp(s)

    if is_circular:
        # Group by modulo-L biological stop coordinate `(strand, (e0 - 1) % L)` (`P0-VBR-1`)
        best_by_stop: Dict[Tuple[str, int], Tuple[int, int, str]] = {}

        for strand, template in (("+", s + s), ("-", rc + rc)):
            for frame in (0, 1, 2):
                # Scan codons step by step in frame
                pos = frame
                active_start: int | None = None
                while pos + 3 <= 2 * l_nt:
                    codon = template[pos : pos + 3]
                    if active_start is None:
                        if codon == "ATG":
                            if pos < l_nt:
                                active_start = pos
                            else:
                                # Past first revolution with no active ORF started in [0, L)
                                break
                    elif codon in _STOP_CODONS:
                        e0 = pos + 3  # Exclusive 0-based end INCLUDING 3 bp stop codon (`P0-VBR-1`)
                        nt_len = e0 - active_start
                        # Rule 1 & Rule 2: 0 <= s0 < L and nt_len <= L (`EC-4.3`)
                        if nt_len <= l_nt:
                            aa_len = (nt_len // 3) - 1
                            if aa_len >= min_aa:
                                aa_seq = str(Seq(template[active_start : e0 - 3]).translate(table=1))
                                if "*" not in aa_seq:
                                    stop_key = (strand, (e0 - 1) % l_nt)
                                    prev = best_by_stop.get(stop_key)
                                    if prev is None or nt_len > (prev[1] - prev[0]):
                                        best_by_stop[stop_key] = (active_start, e0, aa_seq)
                        active_start = None
                    pos += 3

        features: List[AnnotatedCDS] = []
        sorted_items = sorted(
            best_by_stop.items(),
            key=lambda kv: (-len(kv[1][2]), kv[0][0], kv[1][0]),
        )
        for idx, ((strand, _stop_mod), (s0, e0, aa_seq)) in enumerate(sorted_items, start=1):
            intervals = _project_circular_intervals(s0, e0, l_nt, strand)
            features.append(
                AnnotatedCDS(
                    feature_id=f"ORF_{idx}",
                    strand=strand,
                    intervals_1based=intervals,
                    partial_5prime=False,
                    partial_3prime=False,
                    codon_start=1,
                    spans_origin=(e0 > l_nt),
                    aa_sequence=aa_seq,
                    nt_length=e0 - s0,
                    product="polyprotein" if len(aa_seq) >= 500 else f"hypothetical protein ORF{idx}",
                )
            )
        return features

    # Linear genome annotation (supporting complete ORFs and optional 5'/3' boundary partial ORFs `P1-VBR-3`)
    linear_candidates: List[AnnotatedCDS] = []
    orf_counter = 0

    for strand, template in (("+", s), ("-", rc)):
        for frame in (0, 1, 2):
            pos = frame
            first_stop_seen = False
            active_start: int | None = None
            while pos + 3 <= l_nt:
                codon = template[pos : pos + 3]
                if codon in _STOP_CODONS:
                    if active_start is not None:
                        e0 = pos + 3
                        aa_len = ((e0 - active_start) // 3) - 1
                        if aa_len >= min_aa:
                            aa_seq = str(Seq(template[active_start : e0 - 3]).translate(table=1))
                            orf_counter += 1
                            ivs = [(active_start + 1, e0)] if strand == "+" else [(l_nt - active_start, l_nt - e0 + 1)]
                            linear_candidates.append(
                                AnnotatedCDS(
                                    feature_id=f"ORF_{orf_counter}",
                                    strand=strand,
                                    intervals_1based=ivs,
                                    partial_5prime=False,
                                    partial_3prime=False,
                                    codon_start=1,
                                    spans_origin=False,
                                    aa_sequence=aa_seq,
                                    nt_length=e0 - active_start,
                                )
                            )
                        active_start = None
                    elif not first_stop_seen and include_partial_linear:
                        # 5'-partial ORF extending from contig boundary to first stop codon (`P1-VBR-3`)
                        e0 = pos + 3
                        aa_len = (pos - frame) // 3
                        if aa_len >= min_aa:
                            aa_seq = str(Seq(template[frame : e0 - 3]).translate(table=1))
                            orf_counter += 1
                            # Anchor interval start at 1 (or L on '-' strand) and set codon_start = frame + 1 (`P1-VBR-3`)
                            ivs = [(1, e0)] if strand == "+" else [(l_nt, l_nt - e0 + 1)]
                            linear_candidates.append(
                                AnnotatedCDS(
                                    feature_id=f"ORF_{orf_counter}",
                                    strand=strand,
                                    intervals_1based=ivs,
                                    partial_5prime=True,
                                    partial_3prime=False,
                                    codon_start=frame + 1,
                                    spans_origin=False,
                                    aa_sequence=aa_seq,
                                    nt_length=e0,
                                )
                            )
                    first_stop_seen = True
                elif active_start is None and codon == "ATG":
                    active_start = pos
                pos += 3

            # Check trailing 3'-partial ORF running off 3' end of linear contig (`P1-VBR-3`)
            if include_partial_linear:
                if active_start is not None:
                    rem_codons = (l_nt - active_start) // 3
                    if rem_codons >= min_aa:
                        aa_seq = str(Seq(template[active_start : active_start + 3 * rem_codons]).translate(table=1))
                        orf_counter += 1
                        ivs = [(active_start + 1, l_nt)] if strand == "+" else [(l_nt - active_start, 1)]
                        linear_candidates.append(
                            AnnotatedCDS(
                                feature_id=f"ORF_{orf_counter}",
                                strand=strand,
                                intervals_1based=ivs,
                                partial_5prime=False,
                                partial_3prime=True,
                                codon_start=1,
                                spans_origin=False,
                                aa_sequence=aa_seq,
                                nt_length=l_nt - active_start,
                            )
                        )
                elif not first_stop_seen:
                    # Both 5'-partial and 3'-partial (internal fragment of a large polyprotein)
                    rem_codons = (l_nt - frame) // 3
                    if rem_codons >= min_aa:
                        aa_seq = str(Seq(template[frame : frame + 3 * rem_codons]).translate(table=1))
                        orf_counter += 1
                        ivs = [(1, l_nt)] if strand == "+" else [(l_nt, 1)]
                        linear_candidates.append(
                            AnnotatedCDS(
                                feature_id=f"ORF_{orf_counter}",
                                strand=strand,
                                intervals_1based=ivs,
                                partial_5prime=True,
                                partial_3prime=True,
                                codon_start=frame + 1,
                                spans_origin=False,
                                aa_sequence=aa_seq,
                                nt_length=l_nt,
                            )
                        )

    linear_candidates.sort(key=lambda f: -len(f.aa_sequence))
    for i, feat in enumerate(linear_candidates, start=1):
        feat.feature_id = f"ORF_{i}"
    return linear_candidates
