"""CLI Entrypoint and Pipeline Runner for Package 5: `viro-networks`."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
from typing import List, Optional

import pandas as pd

from viro_networks import __version__
from viro_networks.clr_consensus import CohortInsufficientError, infer_consensus_network
from viro_networks.hub_perturbation import evaluate_hub_perturbation
from viro_networks.sloan_neutral import fit_sloan_neutral_model


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viro-networks",
        description="viro_networks: CLR compositional normalization, multi-estimator consensus networks, and Sloan neutral model.",
    )
    parser.add_argument("--version", action="version", version=f"viro-networks {__version__}")
    sub = parser.add_subparsers(dest="subcommand", required=True)

    inf_p = sub.add_parser("infer", help="Infer consensus co-occurrence network and fit Sloan neutral drift model.")
    inf_p.add_argument("--counts", type=pathlib.Path, required=True, help="Raw viral count matrix TSV.")
    inf_p.add_argument("--metadata", type=pathlib.Path, required=False, help="Optional sample metadata TSV.")
    inf_p.add_argument("--bootstraps", type=int, default=100, help="Bootstrap resamples for consensus stability.")
    inf_p.add_argument("--stability", type=float, default=0.80, help="Minimum bootstrap stability threshold.")
    inf_p.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Output directory.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code) if exc.code is not None else 2

    if args.subcommand == "infer":
        if not args.counts.exists():
            print(f"Error: Counts matrix not found: {args.counts}", file=sys.stderr)
            return 2

        try:
            df = pd.read_csv(args.counts, sep="\t", index_col=0)
            taxa_names = [str(c) for c in df.columns]
            mat = df.to_numpy(dtype=float)
        except Exception as exc:
            print(f"Counts table parse error: {exc}", file=sys.stderr)
            return 2

        args.output.mkdir(parents=True, exist_ok=True)
        summary_path = args.output / "network_summary.json"

        try:
            edges, net_diag = infer_consensus_network(
                counts=mat,
                taxa_names=taxa_names,
                n_boot=args.bootstraps,
                stability_threshold=args.stability,
            )
            sloan_rows, sloan_diag = fit_sloan_neutral_model(counts=mat, taxa_names=taxa_names)
            hubs = evaluate_hub_perturbation(
                taxa_names=taxa_names,
                edge_pairs=[(e.taxon_a, e.taxon_b) for e in edges],
            )
        except CohortInsufficientError as exc:
            summary_path.write_text(
                json.dumps({"status": "ABSTAIN_INSUFFICIENT_SAMPLES", "message": str(exc)}, indent=2),
                encoding="utf-8",
            )
            return 4

        edges_tsv = args.output / "consensus_edges.tsv"
        with edges_tsv.open("w", encoding="utf-8", newline="") as ef:
            fields = [
                "taxon_a", "taxon_b", "sign", "sparcc_r", "partial_omega",
                "spearman_rho", "estimators_supporting", "bootstrap_stability",
                "pseudocount_concordance",
            ]
            w = csv.DictWriter(ef, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for e in edges:
                w.writerow(
                    {
                        "taxon_a": e.taxon_a,
                        "taxon_b": e.taxon_b,
                        "sign": e.sign,
                        "sparcc_r": f"{e.sparcc_r:.4f}",
                        "partial_omega": f"{e.partial_omega:.4f}",
                        "spearman_rho": f"{e.spearman_rho:.4f}",
                        "estimators_supporting": e.estimators_supporting,
                        "bootstrap_stability": f"{e.bootstrap_stability:.4f}",
                        "pseudocount_concordance": f"{e.pseudocount_concordance:.4f}",
                    }
                )

        sloan_tsv = args.output / "sloan_neutral_fit.tsv"
        with sloan_tsv.open("w", encoding="utf-8", newline="") as sf:
            s_fields = [
                "taxon", "mean_rel_abundance", "observed_occurrence",
                "predicted_occurrence", "wilson_95_low", "wilson_95_high",
                "p_value", "q_value", "partition_class",
            ]
            sw = csv.DictWriter(sf, fieldnames=s_fields, delimiter="\t")
            sw.writeheader()
            for s in sloan_rows:
                sw.writerow(
                    {
                        "taxon": s.taxon,
                        "mean_rel_abundance": f"{s.mean_rel_abundance:.6f}",
                        "observed_occurrence": f"{s.observed_occurrence:.4f}",
                        "predicted_occurrence": f"{s.predicted_occurrence:.4f}",
                        "wilson_95_low": f"{s.wilson_95_low:.4f}",
                        "wilson_95_high": f"{s.wilson_95_high:.4f}",
                        "p_value": f"{s.p_value:.6e}",
                        "q_value": f"{s.q_value:.6e}",
                        "partition_class": s.partition_class,
                    }
                )

        hubs_tsv = args.output / "topological_hub_candidates.tsv"
        with hubs_tsv.open("w", encoding="utf-8", newline="") as hf:
            h_fields = [
                "taxon", "degree", "degree_centrality",
                "betweenness_centrality", "natural_connectivity_drop_pct",
                "lcc_fraction_after_removal",
            ]
            hw = csv.DictWriter(hf, fieldnames=h_fields, delimiter="\t")
            hw.writeheader()
            for h in hubs:
                hw.writerow(
                    {
                        "taxon": h.taxon,
                        "degree": h.degree,
                        "degree_centrality": f"{h.degree_centrality:.4f}",
                        "betweenness_centrality": f"{h.betweenness_centrality:.4f}",
                        "natural_connectivity_drop_pct": f"{h.natural_connectivity_drop_pct:.2f}",
                        "lcc_fraction_after_removal": f"{h.lcc_fraction_after_removal:.4f}",
                    }
                )

        summary_path.write_text(
            json.dumps(
                {
                    "status": "SUCCESS",
                    "network_diagnostics": net_diag,
                    "sloan_neutral_model": sloan_diag,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
