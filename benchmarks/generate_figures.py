"""Generate Explanatory Schematic Diagrams & Empirical Gate 3 Benchmark Figures (`docs/assets/*.svg`).

Outputs clean, publication-grade vector SVGs into `docs/assets/` for the GitHub README
and optional PNG previews for visual layout verification.
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
for _pkg in ("viro_networks",):
    _p = str(ROOT_DIR / "packages" / _pkg)
    if _p not in sys.path:
        sys.path.insert(0, _p)

from viro_networks.clr_consensus import infer_consensus_network
from viro_networks.sloan_neutral import fit_sloan_neutral_model, sloan_predicted_occurrence

ASSETS_DIR = ROOT_DIR / "docs" / "assets"
RESULTS_DIR = ROOT_DIR / "benchmarks" / "results"
FIXTURES_DIR = ROOT_DIR / "tests" / "fixtures"

PALETTE = {
    "navy": "#1E293B",
    "slate": "#475569",
    "light_bg": "#F8FAFC",
    "card_border": "#CBD5E1",
    "blue": "#2563EB",
    "blue_light": "#DBEAFE",
    "teal": "#0D9488",
    "teal_light": "#CCFBF1",
    "emerald": "#059669",
    "emerald_light": "#D1FAE5",
    "amber": "#D97706",
    "amber_light": "#FEF3C7",
    "purple": "#7C3AED",
    "purple_light": "#EDE9FE",
    "rose": "#E11D48",
    "rose_light": "#FFE4E6",
    "gray": "#94A3B8",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "axes.edgecolor": "#CBD5E1",
    "axes.linewidth": 1.0,
})


def _save_both(fig: plt.Figure, name: str, preview_dir: pathlib.Path | None) -> None:
    ASSETS_DIR.mkdir(parents=True, exist_ok=True)
    svg_path = ASSETS_DIR / f"{name}.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")
    if preview_dir is not None:
        preview_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview_dir / f"{name}.png", format="png", dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {svg_path}")


def _draw_box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    subtitle: str,
    bg: str,
    edge: str,
    title_color: str = "#0F172A",
    title_size: float = 10.2,
    sub_size: float = 8.2,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        facecolor=bg,
        edgecolor=edge,
        linewidth=1.8,
    )
    ax.add_patch(patch)
    if subtitle:
        ax.text(
            x + w / 2.0,
            y + h * 0.70,
            title,
            ha="center",
            va="center",
            fontsize=title_size,
            fontweight="bold",
            color=title_color,
        )
        ax.text(
            x + w / 2.0,
            y + h * 0.31,
            subtitle,
            ha="center",
            va="center",
            fontsize=sub_size,
            color="#334155",
            linespacing=1.28,
        )
    else:
        ax.text(
            x + w / 2.0,
            y + h * 0.50,
            title,
            ha="center",
            va="center",
            fontsize=title_size,
            fontweight="bold",
            color=title_color,
        )


def _draw_arrow(
    ax: plt.Axes,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    color: str = "#475569",
    label: str = "",
    label_dy: float = 0.028,
) -> None:
    arr = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>,head_length=6,head_width=4",
        mutation_scale=1.3,
        linewidth=1.8,
        color=color,
    )
    ax.add_patch(arr)
    if label:
        mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
        ax.text(
            mx,
            my + label_dy,
            label,
            ha="center",
            va="center",
            fontsize=7.6,
            fontweight="bold",
            color=color,
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": color, "lw": 0.8, "alpha": 0.96},
        )


# ============================================================================
# FIGURE 1: Visual Pipeline & Multi-Sample Cohort Blueprint
# ============================================================================


def build_fig1_pipeline_overview(preview_dir: pathlib.Path | None) -> None:
    fig, ax = plt.subplots(figsize=(14.0, 6.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.50,
        0.955,
        "Virome-Suite: End-to-End Multi-Sample & Multi-Replicate Architecture",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        color=PALETTE["navy"],
    )
    ax.text(
        0.50,
        0.905,
        "Five standalone pip-installable packages + one topology-aware cohort orchestrator (virome-suite run-cohort)",
        ha="center",
        va="center",
        fontsize=9.5,
        color=PALETTE["slate"],
    )

    # Stage 0: Input Cohort
    _draw_box(
        ax, 0.015, 0.42, 0.15, 0.28,
        "Study Cohort Input",
        "N Samples x M Replicates\nFASTQ / FASTA + Host Ref\n(study_samplesheet.tsv)",
        PALETTE["light_bg"], PALETTE["navy"], PALETTE["navy"], 9.8, 7.9,
    )

    # Stage 1: viro-drift
    _draw_box(
        ax, 0.20, 0.42, 0.17, 0.28,
        "1. viro-drift",
        "Smart Host Partitioning +\n6-Frame Enzyme Baiting\n+ Cycle-Safe Micro-Assembly",
        PALETTE["blue_light"], PALETTE["blue"], PALETTE["blue"], 10.2, 7.9,
    )
    _draw_arrow(ax, 0.170, 0.56, 0.197, 0.56, PALETTE["navy"])

    # Stage 2: viral-genome-tools (Router)
    _draw_box(
        ax, 0.405, 0.40, 0.185, 0.32,
        "2. viral-genome-tools",
        "TIR Linear Veto +\n1.0x Monomer Trimming +\nPos-1 Landmark Rotation +\nE_junction & P_flank Router",
        PALETTE["teal_light"], PALETTE["teal"], PALETTE["teal"], 10.0, 7.8,
    )
    _draw_arrow(ax, 0.373, 0.56, 0.402, 0.56, PALETTE["blue"])

    # Downstream Modules (shifted right to x=0.715 for clean arrow badge space)
    _draw_box(
        ax, 0.715, 0.64, 0.27, 0.20,
        "3. viralbio-reporting",
        "2x Seam-Wrapping ORF Finder (360° Invariant)\n+ Ready-to-Submit NCBI .tbl/.fsa + SHA-256 Chain",
        PALETTE["emerald_light"], PALETTE["emerald"], PALETTE["emerald"], 10.0, 7.7,
    )
    _draw_box(
        ax, 0.715, 0.37, 0.27, 0.20,
        "4. quasispecies-profiler",
        "Multi-Replicate Beta-Binomial Variant Caller\n+ PCR Jackpot Veto + Pi / Shannon H / dN/dS",
        PALETTE["amber_light"], PALETTE["amber"], PALETTE["amber"], 10.0, 7.7,
    )
    _draw_box(
        ax, 0.715, 0.08, 0.27, 0.21,
        "5. viro-networks  (Cohort N >= 5)",
        "Cross-Sample Count Matrix -> CLR Consensus\nNetwork + Sloan Neutral Model + Hub Ranking",
        PALETTE["purple_light"], PALETTE["purple"], PALETTE["purple"], 10.0, 7.7,
    )

    _draw_arrow(ax, 0.595, 0.63, 0.710, 0.73, PALETTE["emerald"], "Circular & Linear", label_dy=0.035)
    _draw_arrow(ax, 0.595, 0.52, 0.710, 0.48, PALETTE["amber"], "M >= 2 Replicates", label_dy=0.035)
    _draw_arrow(ax, 0.520, 0.395, 0.710, 0.19, PALETTE["purple"], "Cohort Matrix (N x P)", label_dy=-0.038)

    # Bottom-left safety callout badge wrapped cleanly inside box bounds
    _draw_box(
        ax, 0.015, 0.08, 0.575, 0.23,
        "Built-In Biological Safety Guards (Zero Silent Misclassification)",
        "• Linear TIR Veto: Locks linear terminal-repeat templates (0 bases trimmed)\n"
        "• Host-Chimera Flag (P_flank >= 0.15): Flags host-integrated EVE candidates\n"
        "• Low-Depth Abstention (<10x depth or N < 5 samples): Reports UNRESOLVED / ABSTAIN",
        "#FFF1F2", PALETTE["rose"], PALETTE["rose"], 9.5, 7.8,
    )

    _save_both(fig, "fig1_pipeline_overview", preview_dir)


# ============================================================================
# FIGURE 2: Visual "Before vs. After" Infographic of the 4 Core Challenges
# ============================================================================


def build_fig2_five_problems_solved(preview_dir: pathlib.Path | None) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 8.0))
    fig.subplots_adjust(hspace=0.35, wspace=0.18, top=0.90)
    fig.suptitle(
        "What Goes Wrong in Standard Workflows — And How Virome-Suite Fixes It",
        fontsize=14,
        fontweight="bold",
        color=PALETTE["navy"],
        y=0.97,
    )

    # Panel A: Concatemer Loop -> Exact 1.0x Monomer
    ax = axes[0, 0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("A. Stretched Circular Concatemers -> Exact 1.0x Circle", fontsize=11, fontweight="bold", color=PALETTE["teal"], loc="left")
    ax.add_patch(FancyBboxPatch((0.2, 0.3), 9.6, 5.2, boxstyle="round,pad=0.1,rounding_size=0.25", facecolor="#F8FAFC", edgecolor="#CBD5E1"))

    ax.text(2.55, 4.85, "Standard Linear Assembler", ha="center", fontsize=9.2, fontweight="bold", color=PALETTE["rose"])
    ax.add_patch(mpatches.Rectangle((0.55, 3.65), 1.75, 0.65, facecolor="#FCA5A5", edgecolor=PALETTE["rose"], lw=1.4))
    ax.add_patch(mpatches.Rectangle((2.30, 3.65), 1.75, 0.65, facecolor="#F87171", edgecolor=PALETTE["rose"], lw=1.4))
    ax.add_patch(mpatches.Rectangle((4.05, 3.65), 0.55, 0.65, facecolor="#EF4444", edgecolor=PALETTE["rose"], lw=1.4))
    ax.text(2.57, 3.97, "1.0x + 1.0x + 0.25x (2.25x)", ha="center", va="center", fontsize=7.6, fontweight="bold", color="#7F1D1D")
    ax.text(2.55, 2.75, "Duplicated genes + arbitrary start!\n(Also confuses linear TIRs as circles)", ha="center", va="center", fontsize=7.9, color="#475569")

    _draw_arrow(ax, 4.85, 3.35, 6.05, 3.35, PALETTE["teal"], "viral-genome-tools", label_dy=0.35)

    circ = Circle((7.75, 3.15), 1.30, facecolor=PALETTE["teal_light"], edgecolor=PALETTE["teal"], lw=2.6)
    ax.add_patch(circ)
    ax.plot([7.75], [4.45], marker="o", markersize=9, color=PALETTE["amber"])
    ax.text(7.75, 4.88, "Pos 1 Landmark", ha="center", fontsize=8.2, fontweight="bold", color=PALETTE["amber"])
    ax.text(7.75, 3.15, "Exact 1.0x\nMonomer\n(+/- 0 bp)", ha="center", va="center", fontsize=8.6, fontweight="bold", color=PALETTE["navy"])
    ax.text(5.0, 0.78, "Works on circular viruses such as Caulimoviridae, Geminiviridae, Circoviridae, or novel circles", ha="center", fontsize=7.4, style="italic", color="#334155")

    # Panel B: PCR Smudge vs True Replicate Variant
    ax = axes[0, 1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("B. Lab Amplification Jackpots vs. True Mutations", fontsize=11, fontweight="bold", color=PALETTE["amber"], loc="left")
    ax.add_patch(FancyBboxPatch((0.2, 0.3), 9.6, 5.2, boxstyle="round,pad=0.1,rounding_size=0.25", facecolor="#F8FAFC", edgecolor="#CBD5E1"))

    ax.text(0.6, 4.65, "Early PCR Copy Error:", fontsize=8.8, fontweight="bold", color=PALETTE["rose"])
    for idx, (lbl, vaf_str, bg_c, fg_c) in enumerate([
        ("Rep 1", "1.8% Spike", "#FEE2E2", PALETTE["rose"]),
        ("Rep 2", "0.02% (Clean)", "#E2E8F0", "#64748B"),
        ("Rep 3", "0.01% (Clean)", "#E2E8F0", "#64748B"),
    ]):
        bx = 0.6 + idx * 2.05
        ax.add_patch(FancyBboxPatch((bx, 3.55), 1.85, 0.80, boxstyle="round,pad=0.05", facecolor=bg_c, edgecolor=fg_c, lw=1.3))
        ax.text(bx + 0.92, 3.95, f"{lbl}\n{vaf_str}", ha="center", va="center", fontsize=7.6, fontweight="bold", color=fg_c)
    ax.add_patch(FancyBboxPatch((7.0, 3.55), 2.45, 0.80, boxstyle="round,pad=0.05", facecolor="#FFE4E6", edgecolor=PALETTE["rose"], lw=1.6))
    ax.text(8.22, 3.95, "REJECTED\nPCR Jackpot", ha="center", va="center", fontsize=8.0, fontweight="bold", color=PALETTE["rose"])

    ax.text(0.6, 2.35, "True Intra-Host Variant:", fontsize=8.8, fontweight="bold", color=PALETTE["emerald"])
    for idx, (lbl, vaf_str) in enumerate([("Rep 1", "1.1% Variant"), ("Rep 2", "0.9% Variant"), ("Rep 3", "1.0% Variant")]):
        bx = 0.6 + idx * 2.05
        ax.add_patch(FancyBboxPatch((bx, 1.25), 1.85, 0.80, boxstyle="round,pad=0.05", facecolor=PALETTE["emerald_light"], edgecolor=PALETTE["emerald"], lw=1.3))
        ax.text(bx + 0.92, 1.65, f"{lbl}\n{vaf_str}", ha="center", va="center", fontsize=7.6, fontweight="bold", color=PALETTE["emerald"])
    ax.add_patch(FancyBboxPatch((7.0, 1.25), 2.45, 0.80, boxstyle="round,pad=0.05", facecolor=PALETTE["emerald_light"], edgecolor=PALETTE["emerald"], lw=1.6))
    ax.text(8.22, 1.65, "CALLED (q < 0.01)\n95% Profile CI", ha="center", va="center", fontsize=8.0, fontweight="bold", color=PALETTE["emerald"])
    ax.text(5.0, 0.68, "quasispecies-profiler models cross-replicate overdispersion up to >50,000x depth", ha="center", fontsize=7.6, style="italic", color="#334155")

    # Panel C: Chopped Seam Gene -> Intact 2x Wraparound Gene
    ax = axes[1, 0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("C. Genes Cut at Circle Seam -> Intact GenBank Annotation", fontsize=11, fontweight="bold", color=PALETTE["emerald"], loc="left")
    ax.add_patch(FancyBboxPatch((0.2, 0.3), 9.6, 5.2, boxstyle="round,pad=0.1,rounding_size=0.25", facecolor="#F8FAFC", edgecolor="#CBD5E1"))

    ax.text(0.6, 4.7, "Standard Flat-File Gene Finder (Cuts Origin-Spanning ORF in Half):", fontsize=8.4, fontweight="bold", color=PALETTE["rose"])
    ax.add_patch(mpatches.Rectangle((0.6, 3.7), 8.6, 0.55, facecolor="#E2E8F0", edgecolor="#94A3B8"))
    ax.add_patch(mpatches.Rectangle((0.6, 3.7), 1.8, 0.55, facecolor="#FCA5A5", edgecolor=PALETTE["rose"], lw=1.5))
    ax.add_patch(mpatches.Rectangle((7.2, 3.7), 2.0, 0.55, facecolor="#FCA5A5", edgecolor=PALETTE["rose"], lw=1.5))
    ax.text(1.5, 3.97, "3' Half (Dead)", ha="center", va="center", fontsize=7.5, fontweight="bold", color="#7F1D1D")
    ax.text(8.2, 3.97, "5' Half (Dead)", ha="center", va="center", fontsize=7.5, fontweight="bold", color="#7F1D1D")

    ax.text(0.6, 2.65, "viralbio-reporting (Virtual 2x Stitch + Modulo-L Stop Deduplication):", fontsize=8.4, fontweight="bold", color=PALETTE["emerald"])
    ax.add_patch(mpatches.Rectangle((0.6, 1.65), 4.2, 0.55, facecolor="#E2E8F0", edgecolor="#94A3B8"))
    ax.add_patch(mpatches.Rectangle((4.8, 1.65), 4.2, 0.55, facecolor="#F1F5F9", edgecolor="#94A3B8", ls="--"))
    ax.add_patch(mpatches.Rectangle((3.4, 1.65), 2.8, 0.55, facecolor=PALETTE["emerald_light"], edgecolor=PALETTE["emerald"], lw=2.0))
    ax.text(4.8, 1.92, "Intact Full-Length ORF Across Seam", ha="center", va="center", fontsize=7.6, fontweight="bold", color="#065F46")
    ax.text(5.0, 0.78, "100% identical proteins regardless of where the circle is cut -> Valid 5-column NCBI .tbl", ha="center", fontsize=7.6, style="italic", color="#334155")

    # Panel D: Compositional Closure vs Consensus Network
    ax = axes[1, 1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("D. Percentage Distortion -> True Co-Occurrence Network", fontsize=11, fontweight="bold", color=PALETTE["purple"], loc="left")
    ax.add_patch(FancyBboxPatch((0.2, 0.3), 9.6, 5.2, boxstyle="round,pad=0.1,rounding_size=0.25", facecolor="#F8FAFC", edgecolor="#CBD5E1"))

    ax.text(2.4, 4.85, "Raw Count % Trap", ha="center", fontsize=8.8, fontweight="bold", color=PALETTE["rose"])
    ax.add_patch(mpatches.Rectangle((0.8, 2.3), 1.3, 2.1, facecolor="#93C5FD", edgecolor=PALETTE["blue"]))
    ax.add_patch(mpatches.Rectangle((0.8, 1.7), 1.3, 0.6, facecolor="#FDE68A", edgecolor=PALETTE["amber"]))
    ax.text(1.45, 1.30, "Sample 1", ha="center", fontsize=7.6, color="#334155")
    ax.add_patch(mpatches.Rectangle((2.7, 2.0), 1.3, 2.4, facecolor="#FCA5A5", edgecolor=PALETTE["rose"]))
    ax.add_patch(mpatches.Rectangle((2.7, 1.7), 1.3, 0.3, facecolor="#FDE68A", edgecolor=PALETTE["amber"]))
    ax.text(3.35, 1.30, "Sample 2 (Bloom!)", ha="center", fontsize=7.4, color="#334155")

    _draw_arrow(ax, 4.35, 2.95, 5.65, 2.95, PALETTE["purple"], "viro-networks", label_dy=0.35)

    nodes = [(6.55, 3.55, "Virus B"), (8.75, 3.55, "Virus C"), (7.65, 1.85, "Virus D")]
    ax.plot([6.55, 8.75], [3.55, 3.55], color=PALETTE["emerald"], lw=3.2, zorder=1)
    ax.plot([6.55, 7.65], [3.55, 1.85], color="#CBD5E1", lw=1.2, ls=":", zorder=1)
    ax.plot([8.75, 7.65], [3.55, 1.85], color="#CBD5E1", lw=1.2, ls=":", zorder=1)
    ax.text(7.65, 4.55, "Verified Co-Infection\n(3 Estimators + Bootstrap)", ha="center", fontsize=7.6, fontweight="bold", color=PALETTE["emerald"])
    for nx, ny, nlbl in nodes:
        ax.add_patch(Circle((nx, ny), 0.46, facecolor=PALETTE["purple_light"], edgecolor=PALETTE["purple"], lw=1.8, zorder=2))
        ax.text(nx, ny, nlbl, ha="center", va="center", fontsize=7.5, fontweight="bold", color=PALETTE["navy"], zorder=3)
    ax.text(5.0, 0.68, "CLR + SparCC / Graphical / Rank Consensus + Sloan Neutral Drift Baseline", ha="center", fontsize=7.6, style="italic", color="#334155")

    _save_both(fig, "fig2_five_problems_solved", preview_dir)


# ============================================================================
# FIGURE 3: 4-Panel Empirical Gate 3 Benchmark Scorecard (Real Data!)
# ============================================================================


def build_fig3_benchmark_performance(preview_dir: pathlib.Path | None) -> None:
    scorecard = json.loads((RESULTS_DIR / "gate3_master_scorecard.json").read_text(encoding="utf-8"))
    bm1_rows = scorecard["BM1_divergent_rescue"]["rows"]
    bm2 = scorecard["BM2_circular_curation"]
    bm3_rows = scorecard["BM3_variant_calling"]["rows"]
    bm4 = scorecard["BM4_seam_annotation"]
    bm5 = scorecard["BM5_ecological_networks"]

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.4))
    fig.subplots_adjust(hspace=0.38, wspace=0.25, top=0.90)
    fig.suptitle(
        "Empirical Gate 3 Benchmark Validation Across All 5 Virome-Suite Modules",
        fontsize=14,
        fontweight="bold",
        color=PALETTE["navy"],
        y=0.97,
    )

    # Panel A: BM-1
    ax = axes[0, 0]
    divs = [r["divergence_pct"] for r in bm1_rows]
    vd_y = [r["viro_drift_genome_recovery_pct"] for r in bm1_rows]
    nt_y = [r["nucleotide_bait_recovery_pct"] for r in bm1_rows]
    gl_y = [r["unguided_global_nt_search_pct"] for r in bm1_rows]
    ax.plot(divs, vd_y, "-o", color=PALETTE["blue"], lw=2.6, markersize=7, label="viro-drift (6-Frame Catalytic Bait)")
    ax.plot(divs, nt_y, "--s", color=PALETTE["amber"], lw=2.0, markersize=6, label="Nucleotide Seed Baseline (19-mer)")
    ax.plot(divs, gl_y, ":^", color=PALETTE["gray"], lw=1.8, markersize=6, label="Unguided Global NT Search")
    ax.set_title("A. BM-1: Divergent Rescue in 98% Host DNA", fontsize=11, fontweight="bold", color=PALETTE["navy"], loc="left")
    ax.set_xlabel("Sequence Divergence from Known Reference (%)", fontsize=9.5)
    ax.set_ylabel("Target Locus Recovery (%)", fontsize=9.5)
    ax.set_ylim(-5, 110)
    ax.set_xticks(divs)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="center left", fontsize=8.0, frameon=True, facecolor="white")

    # Panel B: BM-2 (ylim raised to 128 so legend never touches 100% bar labels)
    ax = axes[0, 1]
    metrics_b = ["Exact 1.0x Monomer\n(With Tail Indels)", "Linear TIR Veto\nSpecificity (n=50)", "Pos-1 Landmark &\nE_junction AUROC"]
    vgt_vals = [
        bm2["vgt_exact_monomer_recovery_pct"],
        100.0 - bm2["vgt_false_circularization_pct"],
        bm2["e_junction_auroc"] * 100.0,
    ]
    base_vals = [
        bm2["tail_subtraction_baseline_exact_pct"],
        100.0 - bm2["naive_no_tir_veto_false_circ_pct"],
        50.0,
    ]
    x_b = np.arange(len(metrics_b))
    w_b = 0.34
    bars1 = ax.bar(x_b - w_b / 2, vgt_vals, w_b, color=PALETTE["teal"], label="viral-genome-tools")
    bars2 = ax.bar(x_b + w_b / 2, base_vals, w_b, color="#FDA4AF", edgecolor=PALETTE["rose"], label="Naive DTR Trimmer Baseline")
    for b in bars1:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 2.2, f"{b.get_height():.0f}%", ha="center", fontsize=8.3, fontweight="bold", color=PALETTE["teal"])
    for b in bars2:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 2.2, f"{b.get_height():.0f}%", ha="center", fontsize=8.3, color=PALETTE["rose"])
    ax.set_title("B. BM-2: Monomer Trimming & Linear TIR Protection", fontsize=11, fontweight="bold", color=PALETTE["navy"], loc="left")
    ax.set_ylabel("Accuracy / Specificity (%)", fontsize=9.5)
    ax.set_xticks(x_b)
    ax.set_xticklabels(metrics_b, fontsize=8.5)
    ax.set_ylim(0, 130)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)
    ax.legend(loc="upper right", ncol=2, fontsize=7.8, frameon=True, facecolor="white")

    # Panel C: BM-3
    ax = axes[1, 0]
    d_labels = [f"{r['depth_x']:,}x" for r in bm3_rows]
    qp_fdr = [r["qp_empirical_fdr_pct"] for r in bm3_rows]
    sr_fdr = [r["single_rep_binomial_fdr_pct"] for r in bm3_rows]
    pl_fdr = [r["pooled_binomial_fdr_pct"] for r in bm3_rows]
    x_c = np.arange(len(d_labels))
    ax.plot(x_c, sr_fdr, "--s", color=PALETTE["rose"], lw=2.2, markersize=6.5, label="Single-Replicate Binomial (Jackpot-Prone)")
    ax.plot(x_c, pl_fdr, ":^", color=PALETTE["amber"], lw=2.0, markersize=6.5, label="Pooled-Reads Binomial (No Dispersion)")
    ax.plot(x_c, qp_fdr, "-o", color=PALETTE["emerald"], lw=2.8, markersize=7.5, label="quasispecies-profiler (0.0% Empirical FDR)")
    ax.axhline(1.0, color="#64748B", lw=1.0, ls="-.", label="Target FDR Threshold (1%)")
    ax.set_title("C. BM-3: False Discovery Rate Under PCR Jackpots", fontsize=11, fontweight="bold", color=PALETTE["navy"], loc="left")
    ax.set_xlabel("Sequencing Coverage Depth per Replicate", fontsize=9.5)
    ax.set_ylabel("Empirical False Discovery Rate (%)", fontsize=9.5)
    ax.set_xticks(x_c)
    ax.set_xticklabels(d_labels, fontsize=9.0)
    ax.set_ylim(-4, 72)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="center right", fontsize=8.0, frameon=True, facecolor="white")

    # Panel D: BM-4 & BM-5
    ax = axes[1, 1]
    categories = ["BM-4: Intact Proteome\nAcross 30 Circle Cuts", "BM-5: Co-Occurrence F1\n(x100 Score, Closure Stress)"]
    ours_d = [bm4["viralbio_reporting_exact_proteome_pct"], bm5["viro_networks_mean_f1"] * 100.0]
    base_d = [bm4["linear_1x_caller_exact_proteome_pct"], bm5["raw_rel_pearson_mean_f1"] * 100.0]
    x_d = np.arange(len(categories))
    w_d = 0.32
    b_o = ax.bar(x_d - w_d / 2, ours_d, w_d, color=PALETTE["purple"], label="Virome-Suite Module")
    b_b = ax.bar(x_d + w_d / 2, base_d, w_d, color="#CBD5E1", edgecolor="#64748B", label="Standard Linear / Pearson Baseline")
    for b in b_o:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 2.2, f"{b.get_height():.1f}%", ha="center", fontsize=8.8, fontweight="bold", color=PALETTE["purple"])
    for b in b_b:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 2.2, f"{b.get_height():.1f}%", ha="center", fontsize=8.8, color="#475569")
    ax.set_title("D. BM-4 & BM-5: Seam Annotation & Network Recovery", fontsize=11, fontweight="bold", color=PALETTE["navy"], loc="left")
    ax.set_ylabel("Performance Score (%)", fontsize=9.5)
    ax.set_xticks(x_d)
    ax.set_xticklabels(categories, fontsize=9.0)
    ax.set_ylim(0, 122)
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)
    ax.legend(loc="upper right", fontsize=8.2, frameon=True, facecolor="white")

    _save_both(fig, "fig3_benchmark_performance", preview_dir)


# ============================================================================
# FIGURE 4: Cohort Ecological Showcase (Sloan Neutral Curve + Consensus Graph)
# ============================================================================


def build_fig4_cohort_showcase(preview_dir: pathlib.Path | None) -> None:
    df = pd.read_csv(FIXTURES_DIR / "virome_cohort_counts.tsv", sep="\t", index_col=0)
    mat = df.to_numpy(dtype=float)
    taxa = list(df.columns)
    edges, _ = infer_consensus_network(mat, taxa, n_boot=100, stability_threshold=0.80, seed=42)
    sloan_rows, s_diag = fit_sloan_neutral_model(mat, taxa)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.5))
    fig.subplots_adjust(wspace=0.24, top=0.85, bottom=0.14)
    fig.suptitle(
        "Multi-Sample Cohort Output Example: Neutral Drift Curve & Consensus Co-Infection Graph (n = 24 Samples)",
        fontsize=13.2,
        fontweight="bold",
        color=PALETTE["navy"],
        y=0.96,
    )

    p_grid = np.logspace(-2.2, -0.3, 200)
    pred_curve = sloan_predicted_occurrence(p_grid, s_diag["Ne_m"], s_diag["detection_limit"])
    n_eff = int(s_diag["n_samples_used"])
    z = 1.96
    band_lo = np.maximum(0.0, pred_curve - z * np.sqrt(pred_curve * (1.0 - pred_curve) / n_eff))
    band_hi = np.minimum(1.0, pred_curve + z * np.sqrt(pred_curve * (1.0 - pred_curve) / n_eff))

    ax1.fill_between(np.log10(p_grid), band_lo, band_hi, color=PALETTE["blue_light"], alpha=0.65, label="95% Neutral Prediction Band")
    ax1.plot(np.log10(p_grid), pred_curve, color=PALETTE["blue"], lw=2.4, label=f"Sloan Neutral Fit (Nm = {s_diag['Ne_m']:.2f})")

    custom_offsets = {
        "Badnavirus_A": (0.0, -0.075, "center"),
        "Begomovirus_C": (-0.04, 0.04, "right"),
        "Badnavirus_B": (-0.04, -0.055, "right"),
        "Mastrevirus_D": (0.0, 0.05, "center"),
        "Polerovirus_SCYLV": (-0.045, -0.055, "right"),
        "Poacevirus_SCSMV": (0.045, -0.055, "left"),
    }
    for r in sloan_rows:
        is_neutral = r.partition_class == "neutral_drift_consistent"
        col = PALETTE["slate"] if is_neutral else PALETTE["rose"]
        mk = "o" if is_neutral else "D"
        lx = math.log10(r.mean_rel_abundance)
        ax1.scatter([lx], [r.observed_occurrence], color=col, s=75, marker=mk, zorder=4)
        short_lbl = r.taxon.replace("Poacevirus_", "").replace("Polerovirus_", "")
        dx, dy, ha_align = custom_offsets.get(r.taxon, (0.0, 0.045, "center"))
        ax1.text(lx + dx, r.observed_occurrence + dy, short_lbl, ha=ha_align, fontsize=8.0, fontweight="bold", color=col)

    ax1.set_title("A. Sloan Neutral Community Drift Partitioning", fontsize=11, fontweight="bold", color=PALETTE["navy"], loc="left")
    ax1.set_xlabel("Log10 Mean Relative Abundance across Cohort", fontsize=9.5)
    ax1.set_ylabel("Observed Sample Occurrence Frequency", fontsize=9.5)
    ax1.set_ylim(0.0, 1.12)
    ax1.set_xlim(-2.1, -0.3)
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="lower right", fontsize=8.2, frameon=True, facecolor="white")

    # Panel B: Clean Circular Consensus Network
    ax2.set_xlim(-1.50, 1.50)
    ax2.set_ylim(-1.48, 1.38)
    ax2.axis("off")
    ax2.set_title("B. 3-Estimator Bootstrap Consensus Network", fontsize=11, fontweight="bold", color=PALETTE["navy"], loc="left")

    n_t = len(taxa)
    coords = {}
    for i, t in enumerate(taxa):
        ang = (2.0 * math.pi * i / n_t) + (math.pi / 6.0)
        coords[t] = (1.02 * math.cos(ang), 0.92 * math.sin(ang) + 0.05)

    for i in range(n_t):
        for j in range(i + 1, n_t):
            x1, y1 = coords[taxa[i]]
            x2, y2 = coords[taxa[j]]
            ax2.plot([x1, x2], [y1, y2], color="#E2E8F0", lw=1.0, ls=":", zorder=1)

    for e in edges:
        x1, y1 = coords[e.taxon_a]
        x2, y2 = coords[e.taxon_b]
        if e.sign > 0:
            ax2.plot([x1, x2], [y1, y2], color=PALETTE["emerald"], lw=4.2, zorder=3)
            mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
            offset_x = -0.22 if mx < 0 else 0.22
            ax2.text(
                mx + offset_x,
                my,
                f"+r={e.sparcc_r:.2f}\n({int(e.bootstrap_stability * 100)}% boot)",
                ha="center",
                va="center",
                fontsize=7.6,
                fontweight="bold",
                color=PALETTE["emerald"],
                bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": PALETTE["emerald"], "lw": 1.0},
                zorder=5,
            )
        else:
            ax2.plot([x1, x2], [y1, y2], color=PALETTE["rose"], lw=2.0, ls="--", alpha=0.75, zorder=2)

    connected_pos = {e.taxon_a for e in edges if e.sign > 0} | {e.taxon_b for e in edges if e.sign > 0}
    for t in taxa:
        nx, ny = coords[t]
        active = t in connected_pos
        bg = PALETTE["emerald_light"] if active else PALETTE["light_bg"]
        ed = PALETTE["emerald"] if active else PALETTE["slate"]
        ax2.add_patch(Circle((nx, ny), 0.25, facecolor=bg, edgecolor=ed, lw=2.2, zorder=4))
        ax2.text(nx, ny, t.replace("_", "\n"), ha="center", va="center", fontsize=7.6, fontweight="bold", color=PALETTE["navy"], zorder=5)

    # Clean legend bar at bottom of Panel B
    pos_line = mpatches.Patch(facecolor=PALETTE["emerald"], label="Positive Co-Occurrence (boot >= 80%)")
    neg_line = mpatches.Patch(facecolor=PALETTE["rose"], alpha=0.75, label="Negative / Compositional Exclusion")
    ax2.legend(handles=[pos_line, neg_line], loc="lower center", ncol=2, fontsize=7.8, frameon=True, facecolor="white")

    _save_both(fig, "fig4_cohort_showcase", preview_dir)


def main() -> None:
    preview = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else None
    build_fig1_pipeline_overview(preview)
    build_fig2_five_problems_solved(preview)
    build_fig3_benchmark_performance(preview)
    build_fig4_cohort_showcase(preview)
    print("All 4 SVG figures generated cleanly in docs/assets/!")


if __name__ == "__main__":
    main()
