"""Unified Command-Line Interface for Virome_Suite (`virome-suite`)."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
from typing import List, Optional

from orchestrator import __version__
from orchestrator.preflight import format_preflight_report, run_preflight_check
from orchestrator.runner import execute_end_to_end_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="virome-suite",
        description=(
            "Virome_Suite: Modular Python Ecosystem for Divergent Viral Discovery, "
            "Circular Genome Curation, Replicate Quasispecies Profiling, Seam-Aware "
            "NCBI Annotation, and Ecological Co-Occurrence Networks."
        ),
    )
    parser.add_argument("--version", action="version", version=f"virome-suite {__version__}")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # 1. check-env
    p_check = subparsers.add_parser(
        "check-env",
        help="Verify Python dependencies, external bioinformatics tools, and foundation artifacts.",
    )
    p_check.add_argument("--json", action="store_true", help="Emit preflight report as JSON.")

    # 2. audit-summary
    p_audit = subparsers.add_parser(
        "audit-summary",
        help="Summarize the 263-script inventory in scripts_manifest.tsv.",
    )
    p_audit.add_argument(
        "--manifest",
        type=pathlib.Path,
        default=pathlib.Path("scripts_manifest.tsv"),
        help="Path to scripts_manifest.tsv",
    )

    # 3. run
    p_run = subparsers.add_parser(
        "run",
        help="Execute the branching end-to-end Virome_Suite pipeline on input reads or assemblies.",
    )
    p_run.add_argument("--r1", type=pathlib.Path, required=True, help="Forward reads or contigs FASTA/FASTQ.")
    p_run.add_argument("--r2", type=pathlib.Path, required=False, help="Reverse paired-end FASTA/FASTQ.")
    p_run.add_argument("--seeds", type=pathlib.Path, required=False, help="Core viral protein seeds FAA.")
    p_run.add_argument("--host", type=pathlib.Path, required=False, help="Host reference FASTA.")
    p_run.add_argument("--replicates", type=pathlib.Path, required=False, help="Optional replicate variant manifest TSV.")
    p_run.add_argument("--counts", type=pathlib.Path, required=False, help="Optional cohort count matrix TSV.")
    p_run.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        required=True,
        help="Output directory for curated genomes, variants, annotations, and provenance log.",
    )

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code) if exc.code is not None else 2

    if args.subcommand == "check-env":
        report = run_preflight_check()
        if args.json:
            print(json.dumps(report, indent=2))
        else:
            print(format_preflight_report(report))
        return 0

    if args.subcommand == "audit-summary":
        if not args.manifest.exists():
            print(f"Error: Manifest not found at {args.manifest}", file=sys.stderr)
            return 2
        with args.manifest.open("r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        pkg_counts: dict[str, int] = {}
        status_counts: dict[str, int] = {}
        for r in rows:
            pkg = r["target_package"]
            st = "internal_duplicate" if r["status"].startswith("duplicate_of_") else r["status"]
            pkg_counts[pkg] = pkg_counts.get(pkg, 0) + 1
            status_counts[st] = status_counts.get(st, 0) + 1
        print(f"Total Audited Scripts: {len(rows)}")
        print("Breakdown by Target Package:")
        for k, v in sorted(pkg_counts.items()):
            print(f"  - {k:30s}: {v}")
        print("Breakdown by Provenance Status:")
        for k, v in sorted(status_counts.items()):
            print(f"  - {k:30s}: {v}")
        return 0

    if args.subcommand == "run":
        if not args.r1.exists():
            print(f"Error: Input file not found: {args.r1}", file=sys.stderr)
            return 2
        execute_end_to_end_pipeline(
            r1_path=args.r1,
            output_dir=args.output,
            r2_path=args.r2,
            seeds_path=args.seeds,
            host_path=args.host,
            replicate_manifest=args.replicates,
            cohort_counts=args.counts,
        )
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
