"""Two-Pass Priority Multi-Family Circular Landmark Canonicalization (`viral_genome_tools.landmark_rotator`)."""

from __future__ import annotations

import dataclasses
from typing import List, Optional, Tuple

from Bio.Seq import Seq

from viral_genome_tools.monomer_curator import revcomp


@dataclasses.dataclass(frozen=True)
class LandmarkMotif:
    """Curated biological landmark for circular coordinate-1 canonicalization."""

    family: str
    genus: str
    landmark_name: str
    motif_seq: str
    description: str


MULTI_FAMILY_LANDMARK_REGISTRY: List[LandmarkMotif] = [
    LandmarkMotif(
        family="Caulimoviridae",
        genus="Caulimovirus",
        landmark_name="tRNA-Met_PBS_Caulimo",
        motif_seq="TGGTATCAGAGCC",
        description="Caulimovirus 35S pregenomic RNA tRNA-Met primer binding site",
    ),
    LandmarkMotif(
        family="Caulimoviridae",
        genus="Badnavirus",
        landmark_name="tRNA-Met_PBS",
        motif_seq="TGGTATCAGAGC",
        description="Cytosolic plant tRNA-Met primer binding site (ICTV coordinate 1 anchor)",
    ),
    LandmarkMotif(
        family="Geminiviridae",
        genus="Begomovirus",
        landmark_name="StemLoop_Nonanucleotide",
        motif_seq="TAATATTAC",
        description="Conserved rolling-circle origin stem-loop nonanucleotide (TAATATT/AC)",
    ),
    LandmarkMotif(
        family="Geminiviridae",
        genus="Mastrevirus",
        landmark_name="StemLoop_Nonanucleotide",
        motif_seq="TAATATTAC",
        description="Conserved rolling-circle origin stem-loop nonanucleotide (TAATATT/AC)",
    ),
    LandmarkMotif(
        family="Circoviridae",
        genus="Circovirus",
        landmark_name="Origin_Nonanucleotide",
        motif_seq="TAGTATTAC",
        description="Conserved circoviral rolling-circle replication origin nonanucleotide",
    ),
    LandmarkMotif(
        family="Nanoviridae",
        genus="Nanovirus",
        landmark_name="Origin_Nonanucleotide",
        motif_seq="TAGTATTAC",
        description="Conserved multipartite circular ssDNA origin nonanucleotide",
    ),
]


@dataclasses.dataclass
class LandmarkRotationResult:
    """Result of Stage 2 strand normalization and coordinate-1 rotation."""

    rotated_seq: str
    strand_flipped: bool
    rotation_offset_bp: int
    landmark_name: str
    matched_family: str
    landmark_score: float
    orientation_mode: str


def _find_motif_on_circular_strand(
    seq: str,
    motif: str,
    max_mismatches: int = 0,
) -> Optional[int]:
    """Scan circular wrap `seq + seq[:len(motif)-1]` restricting start index to 0 <= pos < len(seq) (`P0-VGT-2`)."""
    p_len = len(seq)
    m_len = len(motif)
    if p_len < m_len or m_len == 0:
        return None

    wrap = seq + seq[: m_len - 1]
    if max_mismatches == 0:
        idx = wrap.find(motif)
        return idx if 0 <= idx < p_len else None

    for pos in range(p_len):
        mism = 0
        for k in range(m_len):
            if wrap[pos + k] != motif[k]:
                mism += 1
                if mism > max_mismatches:
                    break
        if mism <= max_mismatches:
            return pos
    return None


def _find_longest_circular_orf_start(seq: str) -> Tuple[str, bool, int]:
    """Fallback for uncharacterized circular viruses: orient to (+) strand of longest ORF and rotate to ATG."""
    p_len = len(seq)
    rc = revcomp(seq)
    best_len = -1
    best_strand_seq = seq
    best_flipped = False
    best_start = 0

    for flipped, template in ((False, seq), (True, rc)):
        doubled = template + template
        for frame in (0, 1, 2):
            sub = doubled[frame:]
            trim = (len(sub) // 3) * 3
            if trim < 30:
                continue
            prot = str(Seq(sub[:trim]).translate(table=1))
            n_aa = len(prot)
            i = 0
            while i < n_aa:
                if prot[i] == "M":
                    start_nt = frame + i * 3
                    if start_nt < p_len:
                        stop_idx = prot.find("*", i + 1)
                        if stop_idx != -1:
                            end_nt = frame + (stop_idx + 1) * 3
                            orf_nt_len = end_nt - start_nt
                            if orf_nt_len <= p_len and orf_nt_len > best_len:
                                best_len = orf_nt_len
                                best_strand_seq = template
                                best_flipped = flipped
                                best_start = start_nt
                            i = stop_idx + 1
                            continue
                i += 1

    rotated = best_strand_seq[best_start:] + best_strand_seq[:best_start]
    return rotated, best_flipped, best_start


def rotate_to_landmark(
    monomer_seq: str,
    custom_motif: Optional[str] = None,
    registry: Optional[List[LandmarkMotif]] = None,
) -> LandmarkRotationResult:
    """Standardize circular monomer to (+) strand and rotate Coordinate 1 via Two-Pass Priority (`P0-VGT-2`)."""
    s = monomer_seq.upper().strip()
    rc = revcomp(s)

    motifs: List[LandmarkMotif] = []
    if custom_motif:
        motifs.append(
            LandmarkMotif(
                family="Custom",
                genus="Custom",
                landmark_name="Custom_Motif",
                motif_seq=custom_motif.upper().strip(),
                description="User-supplied landmark motif",
            )
        )
    base_reg = list(registry if registry is not None else MULTI_FAMILY_LANDMARK_REGISTRY)
    # Sort longest motif first on ties (`P0-VGT-2` Badnavirus 12bp vs Caulimovirus 13bp)
    base_reg.sort(key=lambda m: -len(m.motif_seq))
    motifs.extend(base_reg)

    # Pass 1: Exact 0-mismatch priority scan across all motifs (+ strand, then - strand)
    for lm in motifs:
        pos_f = _find_motif_on_circular_strand(s, lm.motif_seq, max_mismatches=0)
        if pos_f is not None:
            rot = s[pos_f:] + s[:pos_f]
            return LandmarkRotationResult(
                rotated_seq=rot,
                strand_flipped=False,
                rotation_offset_bp=pos_f,
                landmark_name=lm.landmark_name,
                matched_family=lm.family,
                landmark_score=1.0,
                orientation_mode="biological_landmark_exact",
            )
        pos_r = _find_motif_on_circular_strand(rc, lm.motif_seq, max_mismatches=0)
        if pos_r is not None:
            rot = rc[pos_r:] + rc[:pos_r]
            return LandmarkRotationResult(
                rotated_seq=rot,
                strand_flipped=True,
                rotation_offset_bp=pos_r,
                landmark_name=lm.landmark_name,
                matched_family=lm.family,
                landmark_score=1.0,
                orientation_mode="biological_landmark_exact",
            )

    # Pass 2: 1-mismatch fallback strictly for motifs >= 12 bp (e.g. degeneracies in tRNA-Met PBS)
    for lm in motifs:
        if len(lm.motif_seq) < 12:
            continue
        pos_f = _find_motif_on_circular_strand(s, lm.motif_seq, max_mismatches=1)
        if pos_f is not None:
            rot = s[pos_f:] + s[:pos_f]
            return LandmarkRotationResult(
                rotated_seq=rot,
                strand_flipped=False,
                rotation_offset_bp=pos_f,
                landmark_name=lm.landmark_name,
                matched_family=lm.family,
                landmark_score=0.5,
                orientation_mode="biological_landmark_1mm",
            )
        pos_r = _find_motif_on_circular_strand(rc, lm.motif_seq, max_mismatches=1)
        if pos_r is not None:
            rot = rc[pos_r:] + rc[:pos_r]
            return LandmarkRotationResult(
                rotated_seq=rot,
                strand_flipped=True,
                rotation_offset_bp=pos_r,
                landmark_name=lm.landmark_name,
                matched_family=lm.family,
                landmark_score=0.5,
                orientation_mode="biological_landmark_1mm",
            )

    # Deterministic computational convention fallback for novel uncharacterized viruses
    rot_seq, flipped, orf_pos = _find_longest_circular_orf_start(s)
    return LandmarkRotationResult(
        rotated_seq=rot_seq,
        strand_flipped=flipped,
        rotation_offset_bp=orf_pos,
        landmark_name="none",
        matched_family="Unclassified",
        landmark_score=0.0,
        orientation_mode="canonical_orf_orientation",
    )
