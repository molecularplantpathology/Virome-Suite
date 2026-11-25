"""CLI Entrypoint and Pipeline Runner for Package 3: `quasispecies-profiler`."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
from typing import Any, Dict, List, Optional

from quasispecies_profiler import __version__
from quasispecies_profiler.beta_binomial_caller import call_variants_across_replicates
from quasispecies_profiler.diversity_metrics import (
    compute_nucleotide_diversity_pi,
    compute_shannon_entropy,
)


def _parse_replicate_tsv(manifest_path: pathlib.Path) -> List[Dict[str, Any]]:
    with manifest_path.open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"Replicate table is empty: {manifest_path}")

    cols = list(rows[0].keys())
    k_cols = sorted([c for c in cols if c.endswith("_k")])
    n_cols = sorted([c for c in cols if c.endswith("_n")])
    if not k_cols or len(k_cols) != len(n_cols):
        raise ValueError("Manifest TSV must contain matched rep*_k and rep*_n columns.")

    sites: List[Dict[str, Any]] = []
    for r in rows:
        sites.append(
            {
                "pos": int(r["pos"]),
                "ref": r.get("ref", "N"),
                "alt": r.get("alt", "N"),
                "k": [float(r[c]) for c in k_cols],
                "n": [float(r[c]) for c in n_cols],
            }
        )
    return sites


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="quasispecies-profiler",
        description="quasispecies_profiler: Hierarchical Empirical Bayes Beta-Binomial intra-host variant caller.",
    )
    parser.add_argument("--version", action="version", version=f"quasispecies-profiler {__version__}")
    sub = parser.add_subparsers(dest="subcommand", required=True)

    call_p = sub.add_parser("call", help="Call replicate-aware intra-host variants (iSNVs) with PCR noise filtering.")
    call_p.add_argument("--manifest", type=pathlib.Path, required=True, help="Replicate count table TSV.")
    call_p.add_argument("--ref", type=pathlib.Path, required=False, help="Optional viral reference genome FASTA.")
    call_p.add_argument(
        "--replicate-type",
        choices=["technical", "biological"],
        default="technical",
        help="Replication model.",
    )
    call_p.add_argument("--alpha", type=float, default=0.01, help="Target Benjamini-Hochberg FDR threshold.")
    call_p.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Output directory.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code) if exc.code is not None else 2

    if args.subcommand == "call":
        if not args.manifest.exists():
            print(f"Error: Manifest file not found: {args.manifest}", file=sys.stderr)
            return 2

        try:
            sites = _parse_replicate_tsv(args.manifest)
        except Exception as exc:
            print(f"Manifest parse error: {exc}", file=sys.stderr)
            return 2

        genome_len: Optional[int] = None
        if args.ref is not None and args.ref.exists():
            raw_lines = [
                ln.strip() for ln in args.ref.read_text(encoding="utf-8").splitlines()
                if ln.strip() and not ln.startswith(">")
            ]
            genome_len = sum(len(ln) for ln in raw_lines)

        calls = call_variants_across_replicates(
            sites=sites,
            replicate_type=args.replicate_type,
            alpha_fdr=args.alpha,
        )

        args.output.mkdir(parents=True, exist_ok=True)
        out_tsv = args.output / "variants_called.tsv"
        summary_json = args.output / "quasispecies_summary.json"

        fieldnames = [
            "pos", "ref", "alt", "mean_depth", "vaf_hat", "rho_hat",
            "ci_95_low", "ci_95_high", "lrt_stat", "p_value", "q_value",
            "status", "is_called",
        ]
        with out_tsv.open("w", encoding="utf-8", newline="") as tf:
            writer = csv.DictWriter(tf, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for c in calls:
                writer.writerow(
                    {
                        "pos": c.pos,
                        "ref": c.ref,
                        "alt": c.alt,
                        "mean_depth": f"{c.mean_depth:.1f}",
                        "vaf_hat": f"{c.vaf_hat:.6f}",
                        "rho_hat": f"{c.rho_hat:.6f}",
                        "ci_95_low": f"{c.ci_95_low:.6f}",
                        "ci_95_high": f"{c.ci_95_high:.6f}",
                        "lrt_stat": f"{c.lrt_stat:.4f}",
                        "p_value": f"{c.p_value:.6e}",
                        "q_value": f"{c.q_value:.6e}",
                        "status": c.status,
                        "is_called": c.is_called,
                    }
                )

        passed_vafs = [c.vaf_hat if c.is_called else 0.0 for c in calls]
        depths = [c.mean_depth for c in calls]
        summary = {
            "total_evaluated_loci": len(calls),
            "called_isnv_count": sum(1 for c in calls if c.is_called),
            "rejected_pcr_jackpots": sum(1 for c in calls if c.status == "REJECTED_PCR_JACKPOT"),
            "shannon_entropy_H": compute_shannon_entropy(passed_vafs, genome_length=genome_len),
            "nucleotide_diversity_pi": compute_nucleotide_diversity_pi(passed_vafs, depths, genome_length=genome_len),
        }
        summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
