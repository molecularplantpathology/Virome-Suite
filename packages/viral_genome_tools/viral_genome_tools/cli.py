"""CLI Entrypoint and Pipeline Runner for Package 2: `viral-genome-tools`."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import sys
from typing import List, Optional, Set, Tuple

from Bio import SeqIO

from viral_genome_tools import __version__
from viral_genome_tools.junction_scorer import (
    TopologyDecision,
    classify_topology,
    compute_e_junction_score,
    compute_host_flank_ratio,
    profile_reads_across_seam,
)
from viral_genome_tools.landmark_rotator import rotate_to_landmark
from viral_genome_tools.monomer_curator import revcomp, trim_concatemer_to_monomer


def _load_reads(
    r1_path: Optional[pathlib.Path],
    r2_path: Optional[pathlib.Path],
) -> List[Tuple[str, str, Optional[str]]]:
    if r1_path is None:
        return []
    fmt1 = "fastq" if r1_path.suffix.lower() in (".fq", ".fastq") else "fasta"
    recs1 = list(SeqIO.parse(str(r1_path), fmt1))
    if r2_path is None:
        return [(r.id, str(r.seq), None) for r in recs1]
    fmt2 = "fastq" if r2_path.suffix.lower() in (".fq", ".fastq") else "fasta"
    recs2 = list(SeqIO.parse(str(r2_path), fmt2))
    return [(a.id, str(a.seq), str(b.seq)) for a, b in zip(recs1, recs2)]


def _build_host_kmers(host_path: Optional[pathlib.Path]) -> Set[str]:
    if host_path is None or not host_path.exists():
        return set()
    kmers: Set[str] = set()
    for rec in SeqIO.parse(str(host_path), "fasta"):
        s = str(rec.seq).upper()
        for i in range(0, len(s) - 20):
            sub = s[i : i + 21]
            if "N" not in sub:
                rc = revcomp(sub)
                kmers.add(sub if sub <= rc else rc)
    return kmers


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viral-genome-tools",
        description="viral_genome_tools: Tandem concatemer monomerization, multi-family landmark rotation, and E_junction scoring.",
    )
    parser.add_argument("--version", action="version", version=f"viral-genome-tools {__version__}")
    sub = parser.add_subparsers(dest="subcommand", required=True)

    cur_p = sub.add_parser("curate", help="Curate candidate circular viral contigs to 1.0x canonical monomers.")
    cur_p.add_argument("--assembly", type=pathlib.Path, required=True, help="Input assembled contigs FASTA.")
    cur_p.add_argument("--r1", type=pathlib.Path, required=False, help="Forward reads for junction support.")
    cur_p.add_argument("--r2", type=pathlib.Path, required=False, help="Reverse reads for junction support.")
    cur_p.add_argument("--host", type=pathlib.Path, required=False, help="Optional host FASTA for P_flank check.")
    cur_p.add_argument("--custom-motif", type=str, default=None, help="Optional custom coordinate-1 motif.")
    cur_p.add_argument("--insert-size", type=float, default=350.0, help="Mean paired-end insert size.")
    cur_p.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Output directory.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code) if exc.code is not None else 2

    if args.subcommand == "curate":
        if not args.assembly.exists():
            print(f"Error: Assembly file not found: {args.assembly}", file=sys.stderr)
            return 2

        records = list(SeqIO.parse(str(args.assembly), "fasta"))
        if not records:
            print("Error: Assembly FASTA contains zero sequences.", file=sys.stderr)
            return 2

        read_pairs = _load_reads(args.r1, args.r2)
        host_kmers = _build_host_kmers(args.host)
        args.output.mkdir(parents=True, exist_ok=True)

        out_fasta = args.output / "curated_genomes.fasta"
        report_tsv = args.output / "curation_report.tsv"
        summary_json = args.output / "curation_summary.json"

        rows = []
        with out_fasta.open("w", encoding="utf-8") as fa_out:
            for rec in records:
                raw_seq = str(rec.seq).upper()
                raw_hash = hashlib.sha256(raw_seq.encode("utf-8")).hexdigest()[:16]

                # Stage 1: TIR Veto & Concatemer Trimming (`P0-VGT-1`)
                trim_res = trim_concatemer_to_monomer(raw_seq)

                if trim_res.has_terminal_inverted_repeat:
                    final_seq = trim_res.monomer_seq
                    lm_name = "none_TIR_veto"
                    lm_family = "Linear_TIR"
                    lm_score = 0.0
                    ori_mode = "linear_unmodified"
                else:
                    # Stage 2: Multi-Family Landmark Rotation (`P0-VGT-2`)
                    rot_res = rotate_to_landmark(trim_res.monomer_seq, custom_motif=args.custom_motif)
                    final_seq = rot_res.rotated_seq
                    lm_name = rot_res.landmark_name
                    lm_family = rot_res.matched_family
                    lm_score = rot_res.landmark_score
                    ori_mode = rot_res.orientation_mode

                cur_hash = hashlib.sha256(final_seq.encode("utf-8")).hexdigest()[:16]

                # Stage 3: Read-Geometry E_junction & Orthogonal P_flank (`P1-VGT-3`)
                if read_pairs and not trim_res.has_terminal_inverted_repeat:
                    ev = profile_reads_across_seam(final_seq, read_pairs, host_kmers=host_kmers)
                    e_score, _ = compute_e_junction_score(
                        n_split=ev.n_split,
                        n_spanning_pairs=ev.n_spanning_pairs,
                        d_origin=ev.d_origin,
                        d_median=ev.d_median,
                        landmark_score=lm_score,
                        read_length=ev.read_length,
                        mean_insert_size=args.insert_size,
                    )
                    p_flank = compute_host_flank_ratio(ev.n_split, ev.n_host_viral_chimera)
                    decision = classify_topology(
                        e_junction_score=e_score,
                        d_median=ev.d_median,
                        p_flank=p_flank,
                        has_tir=False,
                    )
                    d_med = round(ev.d_median, 2)
                elif trim_res.has_terminal_inverted_repeat:
                    e_score = 0.0
                    p_flank = None
                    decision = TopologyDecision.LINEAR_CANDIDATE
                    d_med = 0.0
                else:
                    # Assembly-only mode (no FASTQ reads provided): score structural DTR + landmark, abstain on read-level topology
                    e_score, _ = compute_e_junction_score(
                        n_split=0,
                        n_spanning_pairs=0,
                        d_origin=0.0,
                        d_median=0.0,
                        landmark_score=lm_score,
                    )
                    p_flank = None
                    decision = TopologyDecision.UNRESOLVED_TOPOLOGY
                    d_med = 0.0

                fa_out.write(
                    f">{rec.id} length={len(final_seq)} topology={decision.value} "
                    f"landmark={lm_name} e_junction={e_score:.4f}\n{final_seq}\n"
                )
                rows.append(
                    {
                        "contig_id": rec.id,
                        "original_length": trim_res.original_length,
                        "curated_length": len(final_seq),
                        "was_trimmed": trim_res.was_trimmed,
                        "has_tir": trim_res.has_terminal_inverted_repeat,
                        "matched_family": lm_family,
                        "landmark_name": lm_name,
                        "orientation_mode": ori_mode,
                        "d_median": d_med,
                        "e_junction_score": f"{e_score:.6f}",
                        "p_flank": "NA" if p_flank is None else f"{p_flank:.4f}",
                        "topology_decision": decision.value,
                        "raw_sha256": raw_hash,
                        "curated_sha256": cur_hash,
                    }
                )

        with report_tsv.open("w", encoding="utf-8", newline="") as tf:
            writer = csv.DictWriter(tf, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        summary_json.write_text(json.dumps({"records": rows}, indent=2), encoding="utf-8")
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
