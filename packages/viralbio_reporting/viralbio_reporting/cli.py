"""CLI Entrypoint and Pipeline Runner for Package 4: `viralbio-reporting`."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import List, Optional

from Bio import SeqIO

from viralbio_reporting import __version__
from viralbio_reporting.circular_annotator import annotate_genome_orfs
from viralbio_reporting.ncbi_packager import build_ncbi_submission_package


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viralbio-reporting",
        description="viralbio_reporting: Modulo-L circular seam wraparound ORF annotation and NCBI table2asn package builder.",
    )
    parser.add_argument("--version", action="version", version=f"viralbio-reporting {__version__}")
    sub = parser.add_subparsers(dest="subcommand", required=True)

    ann_p = sub.add_parser("annotate", help="Annotate circular or linear viral genomes and generate GenBank .tbl/.fsa.")
    ann_p.add_argument("--fasta", type=pathlib.Path, required=True, help="Curated viral genome FASTA.")
    ann_p.add_argument("--topology", choices=["circular", "linear"], default="circular", help="Genome topology.")
    ann_p.add_argument("--min-aa", type=int, default=60, help="Minimum ORF length in amino acids.")
    ann_p.add_argument("--organism", type=str, default="Unclassified virus", help="Organism name for GenBank header.")
    ann_p.add_argument("--host", type=str, default="Plant host", help="Host organism qualifier.")
    ann_p.add_argument("--country", type=str, default="Unspecified", help="Geographic location qualifier.")
    ann_p.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Output directory.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code) if exc.code is not None else 2

    if args.subcommand == "annotate":
        if not args.fasta.exists():
            print(f"Error: Input FASTA not found: {args.fasta}", file=sys.stderr)
            return 2

        records = list(SeqIO.parse(str(args.fasta), "fasta"))
        if not records:
            print("Error: Input FASTA contains zero sequences.", file=sys.stderr)
            return 2

        args.output.mkdir(parents=True, exist_ok=True)
        summaries = []
        for rec in records:
            seq_str = str(rec.seq).upper()
            orfs = annotate_genome_orfs(seq_str, topology=args.topology, min_aa=args.min_aa)
            pkg_res = build_ncbi_submission_package(
                seq_id=rec.id,
                sequence=seq_str,
                features=orfs,
                output_dir=args.output / rec.id if len(records) > 1 else args.output,
                topology=args.topology,
                organism=args.organism,
                host=args.host,
                country=args.country,
            )
            summaries.append(pkg_res)

        (args.output / "annotation_summary.json").write_text(
            json.dumps({"records": summaries}, indent=2),
            encoding="utf-8",
        )
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
