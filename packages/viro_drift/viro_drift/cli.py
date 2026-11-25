"""CLI Entrypoint and Pipeline Runner for Package 1: `viro-drift`."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from viro_drift import __version__
from viro_drift.host_partition import build_host_kmer_index, partition_host_reads
from viro_drift.micro_assembler import assemble_viral_contigs
from viro_drift.recruitment import iterative_recruit_reads
from viro_drift.translation_bait import bait_reads_with_seeds, check_external_bait_tool


def _detect_seq_format(path: pathlib.Path) -> str:
    suffix = path.suffix.lower()
    if suffix in (".fq", ".fastq"):
        return "fastq"
    with path.open("r", encoding="utf-8") as f:
        first = f.readline().strip()
    return "fastq" if first.startswith("@") else "fasta"


def _load_read_pairs(
    r1_path: pathlib.Path,
    r2_path: Optional[pathlib.Path] = None,
) -> List[Tuple[str, str, Optional[str]]]:
    fmt1 = _detect_seq_format(r1_path)
    recs1 = list(SeqIO.parse(str(r1_path), fmt1))
    if r2_path is None:
        return [(r.id.split("/")[0], str(r.seq), None) for r in recs1]
    fmt2 = _detect_seq_format(r2_path)
    recs2 = list(SeqIO.parse(str(r2_path), fmt2))
    if len(recs1) != len(recs2):
        raise ValueError(f"Mismatched paired read counts: R1={len(recs1)}, R2={len(recs2)}")
    return [
        (a.id.split("/")[0], str(a.seq), str(b.seq))
        for a, b in zip(recs1, recs2)
    ]


def _load_fasta_dict(path: pathlib.Path) -> Dict[str, str]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No valid FASTA records found in {path}")
    return {r.id: str(r.seq) for r in records}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viro-drift",
        description="viro_drift: 6-frame translated protein-seeded viral read baiting and local micro-assembly.",
    )
    parser.add_argument("--version", action="version", version=f"viro-drift {__version__}")
    sub = parser.add_subparsers(dest="subcommand", required=True)

    run_p = sub.add_parser("run", help="Bait divergent viral reads and execute local micro-assembly.")
    run_p.add_argument("--r1", type=pathlib.Path, required=True, help="Forward reads (FASTA/FASTQ).")
    run_p.add_argument("--r2", type=pathlib.Path, required=False, help="Reverse reads (FASTA/FASTQ).")
    run_p.add_argument("--seeds", type=pathlib.Path, required=True, help="Core viral protein seeds (FAA).")
    run_p.add_argument("--host", type=pathlib.Path, required=False, help="Optional host reference FASTA for masking.")
    run_p.add_argument("--mode", choices=["kmer_seed", "diamond", "hmmer"], default="kmer_seed")
    run_p.add_argument("--max-rounds", type=int, default=5, help="Maximum iterative recruitment rounds.")
    run_p.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Output directory.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code) if exc.code is not None else 2

    if args.subcommand == "run":
        if not args.r1.exists() or not args.seeds.exists():
            print("Error: Input reads or seed file does not exist.", file=sys.stderr)
            return 2
        if args.r2 is not None and not args.r2.exists():
            print("Error: Specified --r2 file does not exist.", file=sys.stderr)
            return 2

        if args.mode in ("diamond", "hmmer"):
            try:
                check_external_bait_tool(args.mode)
            except FileNotFoundError as exc:
                print(f"External binary missing: {exc}", file=sys.stderr)
                return 3

        try:
            read_pairs = _load_read_pairs(args.r1, args.r2)
            seeds = _load_fasta_dict(args.seeds)
            host_seqs = list(_load_fasta_dict(args.host).values()) if (args.host and args.host.exists()) else []
        except Exception as exc:
            print(f"Input validation error: {exc}", file=sys.stderr)
            return 2

        args.output.mkdir(parents=True, exist_ok=True)
        host_kmers = build_host_kmer_index(host_seqs, k=21) if host_seqs else set()
        partition = partition_host_reads(read_pairs, host_kmers, k=21)

        seed_ids, seed_hits = bait_reads_with_seeds(
            partition.viral_candidate,
            seeds=seeds,
            mode=args.mode,
            seeds_path=args.seeds,
            work_dir=args.output,
        )

        out_fasta = args.output / "assembled_contigs.fasta"
        summary_json = args.output / "viro_drift_summary.json"

        if not seed_ids:
            out_fasta.write_text("", encoding="utf-8")
            summary = {
                "status": "NO_VIRAL_CANDIDATES",
                "partition": partition.summary_counts(),
                "initial_seed_reads": 0,
                "recruited_pairs": 0,
                "assembled_contigs": 0,
            }
            summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
            return 0

        rec = iterative_recruit_reads(
            partition.viral_candidate,
            initial_seed_ids=seed_ids,
            host_kmer_mask=host_kmers,
            max_rounds=args.max_rounds,
        )
        contigs = assemble_viral_contigs(rec.recruited_pairs)

        with out_fasta.open("w", encoding="utf-8") as f:
            for c in contigs:
                f.write(
                    f">{c.contig_id} length={len(c.sequence)} depth={c.mean_kmer_depth} "
                    f"circular_cycle={c.is_circular_cycle}\n{c.sequence}\n"
                )

        summary = {
            "status": "CANDIDATES_ASSEMBLED" if contigs else "SEEDS_RECRUITED_NO_CONTIG",
            "partition": partition.summary_counts(),
            "initial_seed_reads": len(seed_ids),
            "matched_seed_domains": sorted(set(seed_hits.values())),
            "recruited_pairs": len(rec.recruited_pairs),
            "rounds_executed": rec.rounds_executed,
            "stop_reason": rec.stop_reason,
            "filtered_host_decoys_count": rec.filtered_host_decoys_count,
            "assembled_contigs": len(contigs),
        }
        summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
