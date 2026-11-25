"""6-Frame Translated Catalytic Seed Baiting (`viro_drift.translation_bait`)."""

from __future__ import annotations

import math
import pathlib
import shutil
import subprocess
from typing import Dict, Iterable, List, Optional, Set, Tuple

from Bio.Seq import Seq

from viro_drift.host_partition import revcomp

# 11-Class Chemical Equivalence Reduced Amino Acid Alphabet (`P0-VD-2`)
_AA_REDUCED_GROUPS = {
    "I": "L", "L": "L", "M": "L", "V": "L",
    "F": "F", "Y": "F", "W": "F",
    "K": "K", "R": "K",
    "D": "D", "E": "D",
    "S": "S", "T": "S",
    "N": "N", "Q": "N",
    "A": "A",
    "C": "C",
    "G": "G",
    "H": "H",
    "P": "P",
}


def shannon_entropy_nt(window: str) -> float:
    """Compute nucleotide Shannon entropy in bits over {A, C, G, T}."""
    counts = [0, 0, 0, 0]
    total = 0
    for ch in window.upper():
        if ch == "A":
            counts[0] += 1
            total += 1
        elif ch == "C":
            counts[1] += 1
            total += 1
        elif ch == "G":
            counts[2] += 1
            total += 1
        elif ch == "T":
            counts[3] += 1
            total += 1
    if total == 0:
        return 0.0
    ent = 0.0
    for c in counts:
        if c > 0:
            p = c / float(total)
            ent -= p * math.log2(p)
    return ent


def mask_low_complexity_nt(seq: str, window_size: int = 24, min_entropy_bits: float = 1.20) -> str:
    """Mask low-complexity sliding windows (H < 1.20 bits) to 'N' (`P0-VD-1`)."""
    s = seq.upper()
    n = len(s)
    if n == 0:
        return ""
    if n < window_size:
        return "N" * n if shannon_entropy_nt(s) < min_entropy_bits else s

    mask = [False] * n
    for i in range(0, n - window_size + 1):
        if shannon_entropy_nt(s[i : i + window_size]) < min_entropy_bits:
            for j in range(i, i + window_size):
                mask[j] = True
    return "".join("N" if mask[i] else s[i] for i in range(n))


def translate_6frames(seq: str) -> List[Tuple[str, str]]:
    """Translate nucleotide sequence across all 6 reading frames (+0,+1,+2,-0,-1,-2)."""
    s = seq.upper().replace("U", "T")
    rc = revcomp(s)
    frames: List[Tuple[str, str]] = []
    for strand_label, template in (("+", s), ("-", rc)):
        for offset in (0, 1, 2):
            sub = template[offset:]
            trim_len = (len(sub) // 3) * 3
            if trim_len < 15:
                frames.append((f"{strand_label}{offset}", ""))
                continue
            aa = str(Seq(sub[:trim_len]).translate(table=1))
            frames.append((f"{strand_label}{offset}", aa))
    return frames


def _reduce_peptide(peptide: str) -> str:
    return "".join(_AA_REDUCED_GROUPS.get(ch, "X") for ch in peptide.upper())


class ReducedAlphabetSeedIndex:
    """Collinear reduced-alphabet and exact hexapeptide viral seed index (`P0-VD-1`, `P0-VD-2`)."""

    def __init__(self, seeds: Dict[str, str], k_reduced: int = 5, k_exact: int = 6) -> None:
        self.k_reduced = k_reduced
        self.k_exact = k_exact
        # reduced_kmer -> list of (seed_id, pos)
        self.reduced_map: Dict[str, List[Tuple[str, int]]] = {}
        # exact_kmer -> set of seed_ids
        self.exact_map: Dict[str, Set[str]] = {}

        for seed_id, aa_seq in seeds.items():
            clean_aa = aa_seq.upper().strip()
            for segment in self._valid_segments(clean_aa):
                seg_seq, seg_offset = segment
                for i in range(0, len(seg_seq) - k_exact + 1):
                    ek = seg_seq[i : i + k_exact]
                    self.exact_map.setdefault(ek, set()).add(seed_id)
                red_seg = _reduce_peptide(seg_seq)
                for i in range(0, len(red_seg) - k_reduced + 1):
                    rk = red_seg[i : i + k_reduced]
                    if "X" not in rk:
                        self.reduced_map.setdefault(rk, []).append((seed_id, seg_offset + i))

    @staticmethod
    def _valid_segments(aa_seq: str) -> List[Tuple[str, int]]:
        """Split translated polypeptide on '*' and 'X' so masked/stop residues never match (`P0-VD-1`)."""
        segments: List[Tuple[str, int]] = []
        start = 0
        n = len(aa_seq)
        i = 0
        while i < n:
            if aa_seq[i] in ("*", "X"):
                if i > start:
                    segments.append((aa_seq[start:i], start))
                start = i + 1
            i += 1
        if n > start:
            segments.append((aa_seq[start:n], start))
        return segments

    def match_read(self, nt_seq: str) -> Optional[Tuple[str, str]]:
        """Return (matched_seed_id, frame_label) if read matches any seed with collinear consistency."""
        masked = mask_low_complexity_nt(nt_seq)
        if masked.count("N") == len(masked):
            return None

        kr = self.k_reduced
        ke = self.k_exact

        for frame_label, aa_full in translate_6frames(masked):
            if not aa_full:
                continue
            for seg_seq, seg_offset in self._valid_segments(aa_full):
                if len(seg_seq) < kr:
                    continue
                # 1. Exact contiguous match >= 6 aa (`P0-VD-2` Criterion 1)
                for i in range(0, len(seg_seq) - ke + 1):
                    ek = seg_seq[i : i + ke]
                    if ek in self.exact_map:
                        return next(iter(self.exact_map[ek])), frame_label

                # 2. Collinear reduced-alphabet double 5-mer match (`P0-VD-2` Criterion 2)
                red_seg = _reduce_peptide(seg_seq)
                hits_by_seed: Dict[str, List[Tuple[int, int]]] = {}
                for i in range(0, len(red_seg) - kr + 1):
                    rk = red_seg[i : i + kr]
                    if "X" in rk or rk not in self.reduced_map:
                        continue
                    q_pos = seg_offset + i
                    for seed_id, p_pos in self.reduced_map[rk]:
                        prev_list = hits_by_seed.setdefault(seed_id, [])
                        for q1, p1 in prev_list:
                            if (q_pos - q1) >= kr and abs((q_pos - q1) - (p_pos - p1)) <= 3:
                                return seed_id, frame_label
                        prev_list.append((q_pos, p_pos))
        return None


def check_external_bait_tool(mode: str) -> str:
    """Verify external binary availability when diamond/hmmer mode is explicitly requested."""
    bin_name = "diamond" if mode == "diamond" else "hmmsearch"
    exe = shutil.which(bin_name)
    if exe is None:
        raise FileNotFoundError(
            f"Requested external mode '{mode}' requires '{bin_name}' in PATH."
        )
    return exe


def bait_reads_with_seeds(
    candidate_pairs: Iterable[Tuple[str, str, Optional[str]]],
    seeds: Dict[str, str],
    mode: str = "kmer_seed",
    seeds_path: Optional[pathlib.Path] = None,
    work_dir: Optional[pathlib.Path] = None,
) -> Tuple[Set[str], Dict[str, str]]:
    """Bait divergent viral candidate reads in 6-frame protein space.

    Returns (recruited_read_ids, seed_hit_map).
    """
    pairs_list = list(candidate_pairs)
    if mode in ("diamond", "hmmer"):
        exe = check_external_bait_tool(mode)
        if work_dir is not None and seeds_path is not None and mode == "diamond" and pairs_list:
            db_prefix = work_dir / "seed_db"
            query_fa = work_dir / "candidates_query.fasta"
            out_tsv = work_dir / "diamond_hits.tsv"
            with query_fa.open("w", encoding="utf-8") as qf:
                for rid, r1, r2 in pairs_list:
                    qf.write(f">{rid}/1\n{r1}\n")
                    if r2:
                        qf.write(f">{rid}/2\n{r2}\n")
            subprocess.run([exe, "makedb", "--in", str(seeds_path), "-d", str(db_prefix), "--quiet"], check=True)
            subprocess.run(
                [
                    exe, "blastx", "-d", str(db_prefix), "-q", str(query_fa),
                    "-o", str(out_tsv), "--ultra-sensitive", "-e", "1e-4",
                    "--outfmt", "6", "qseqid", "sseqid", "evalue", "--quiet",
                ],
                check=True,
            )
            recruited: Set[str] = set()
            hit_map: Dict[str, str] = {}
            if out_tsv.exists():
                for line in out_tsv.read_text(encoding="utf-8").splitlines():
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        base_id = parts[0].rsplit("/", 1)[0]
                        recruited.add(base_id)
                        hit_map[base_id] = parts[1]
            return recruited, hit_map

    index = ReducedAlphabetSeedIndex(seeds)
    recruited_ids: Set[str] = set()
    seed_hits: Dict[str, str] = {}
    for rid, r1, r2 in pairs_list:
        m1 = index.match_read(r1)
        if m1 is not None:
            recruited_ids.add(rid)
            seed_hits[rid] = m1[0]
            continue
        if r2 is not None:
            m2 = index.match_read(r2)
            if m2 is not None:
                recruited_ids.add(rid)
                seed_hits[rid] = m2[0]
    return recruited_ids, seed_hits
