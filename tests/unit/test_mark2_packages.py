"""Comprehensive 46-Test Verification & Adversarial Edge-Case Suite for Completion Mark 2.

Covers Choice B AST/Wheel Isolation (`T0.1`-`T0.2`) and all 5 standalone packages plus
the unified `virome-suite` orchestrator (`T1.1`-`T5.8`), verifying every post-review
mathematical resolution (`P0-PKG-1` through `P1-VN-3`).
"""

from __future__ import annotations

import ast
import csv
import json
import math
import pathlib
import random
import sys
import unittest
from typing import Dict, List, Tuple
from unittest.mock import patch

import numpy as np
import pandas as pd
from Bio import SeqIO

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent
for _pkg in (
    "viro_drift",
    "viral_genome_tools",
    "quasispecies_profiler",
    "viralbio_reporting",
    "viro_networks",
):
    _pkg_path = str(ROOT_DIR / "packages" / _pkg)
    if _pkg_path not in sys.path:
        sys.path.insert(0, _pkg_path)
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

FIXTURES_DIR = ROOT_DIR / "tests" / "fixtures"

# Package 1: viro-drift
from viro_drift.cli import main as viro_drift_main
from viro_drift.host_partition import build_host_kmer_index, partition_host_reads, revcomp as vd_revcomp
from viro_drift.micro_assembler import assemble_viral_contigs
from viro_drift.recruitment import iterative_recruit_reads
from viro_drift.translation_bait import (
    ReducedAlphabetSeedIndex,
    bait_reads_with_seeds,
    mask_low_complexity_nt,
    shannon_entropy_nt,
    translate_6frames,
)

# Package 2: viral-genome-tools
from viral_genome_tools.cli import main as vgt_main
from viral_genome_tools.junction_scorer import (
    TopologyDecision,
    classify_topology,
    compute_coverage_uniformity,
    compute_e_junction_score,
    compute_host_flank_ratio,
    expected_spanning_pairs,
    expected_split_reads,
    profile_reads_across_seam,
)
from viral_genome_tools.landmark_rotator import rotate_to_landmark
from viral_genome_tools.monomer_curator import (
    detect_terminal_inverted_repeat,
    revcomp as vgt_revcomp,
    trim_concatemer_to_monomer,
)

# Package 3: quasispecies-profiler
from quasispecies_profiler.beta_binomial_caller import (
    call_variants_across_replicates,
    estimate_background_error_and_prior,
    log_beta_binomial_pmf,
)
from quasispecies_profiler.cli import _parse_replicate_tsv, main as qp_main
from quasispecies_profiler.diversity_metrics import (
    compute_nei_gojobori_dnds,
    compute_nucleotide_diversity_pi,
    compute_shannon_entropy,
)

# Package 4: viralbio-reporting
from viralbio_reporting.circular_annotator import (
    AnnotatedCDS,
    annotate_genome_orfs,
    reconstruct_cds_nt_from_intervals,
    revcomp as vbr_revcomp,
    validate_cds_feature,
)
from viralbio_reporting.cli import main as vbr_main
from viralbio_reporting.ncbi_packager import (
    SequenceProvenanceChain,
    build_ncbi_submission_package,
    format_ncbi_feature_table,
    validate_tbl_syntax,
)

# Package 5: viro-networks
from viro_networks.cli import main as vn_main
from viro_networks.clr_consensus import (
    CohortInsufficientError,
    clr_transform,
    compute_sparcc_matrix,
    infer_consensus_network,
)
from viro_networks.hub_perturbation import evaluate_hub_perturbation
from viro_networks.sloan_neutral import fit_sloan_neutral_model

# Master Orchestrator: virome-suite
from orchestrator import audit_inventory
from orchestrator.cli import main as orchestrator_main


AA_TO_CODON: Dict[str, str] = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "TTA",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}

ALT_CODON: Dict[str, str] = {
    "A": "GCC", "C": "TGC", "D": "GAC", "E": "GAG", "F": "TTC",
    "G": "GGC", "H": "CAC", "I": "ATC", "K": "AAG", "L": "CTG",
    "M": "ATG", "N": "AAC", "P": "CCC", "Q": "CAG", "R": "CGG",
    "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAC",
}


def _backtranslate_peptide(peptide: str, use_alt_synonymous: bool = False) -> str:
    table = ALT_CODON if use_alt_synonymous else AA_TO_CODON
    return "".join(table[aa] for aa in peptide)


def _load_fasta_map(path: pathlib.Path) -> Dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(path), "fasta")}


# ============================================================================
# GROUP 0: Choice B Packaging & Dependency Isolation (2 Tests: T0.1 - T0.2)
# ============================================================================


def test_T0_1_ast_static_import_isolation() -> None:
    """T0.1: Zero cross-imports from `orchestrator` or sibling packages inside `packages/*`."""
    all_namespaces = {
        "viro_drift",
        "viral_genome_tools",
        "quasispecies_profiler",
        "viralbio_reporting",
        "viro_networks",
    }
    for pkg_ns in sorted(all_namespaces):
        pkg_root = ROOT_DIR / "packages" / pkg_ns
        forbidden = (all_namespaces - {pkg_ns}) | {"orchestrator"}
        py_files = list(pkg_root.rglob("*.py"))
        assert len(py_files) >= 3, f"Expected modules in {pkg_root}"
        for py_file in py_files:
            tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        top = alias.name.split(".")[0]
                        assert top not in forbidden, f"Forbidden import {top} in {py_file}"
                elif isinstance(node, ast.ImportFrom) and node.module:
                    top = node.module.split(".")[0]
                    assert top not in forbidden, f"Forbidden import {top} in {py_file}"


def test_T0_2_wheel_namespace_ownership_and_versions() -> None:
    """T0.2: Each Choice B `pyproject.toml` owns strictly its own namespace at v0.1.0 (`P0-PKG-1`)."""
    expected = {
        ROOT_DIR / "packages" / "viro_drift" / "pyproject.toml": ("viro-drift", 'packages = ["viro_drift"]'),
        ROOT_DIR / "packages" / "viral_genome_tools" / "pyproject.toml": (
            "viral-genome-tools",
            'packages = ["viral_genome_tools"]',
        ),
        ROOT_DIR / "packages" / "quasispecies_profiler" / "pyproject.toml": (
            "quasispecies-profiler",
            'packages = ["quasispecies_profiler"]',
        ),
        ROOT_DIR / "packages" / "viralbio_reporting" / "pyproject.toml": (
            "viralbio-reporting",
            'packages = ["viralbio_reporting"]',
        ),
        ROOT_DIR / "packages" / "viro_networks" / "pyproject.toml": (
            "viro-networks",
            'packages = ["viro_networks"]',
        ),
        ROOT_DIR / "pyproject.toml": ("virome-suite", 'packages = ["orchestrator"]'),
    }
    for toml_path, (dist_name, pkg_decl) in expected.items():
        text = toml_path.read_text(encoding="utf-8")
        assert f'name = "{dist_name}"' in text
        assert 'version = "0.1.0"' in text
        assert pkg_decl in text


# ============================================================================
# GROUP 1: Package 1 `viro-drift` (8 Tests: T1.1 - T1.8)
# ============================================================================


def test_T1_1_badnavirus_rt_6frame_seed_baiting() -> None:
    """T1.1: 6-frame seed baiting achieves 100% recall on Badnavirus RT reads across frames/strands."""
    seeds = _load_fasta_map(FIXTURES_DIR / "core_viral_seeds.faa")
    rt_nt = _backtranslate_peptide(seeds["Badnavirus_RT_Core"], use_alt_synonymous=True)
    pairs = [
        ("r_f0", "AC" + rt_nt[:140] + "GTACGT", None),
        ("r_rc", vd_revcomp("TGCA" + rt_nt[15:160]), None),
    ]
    matched_ids, hit_map = bait_reads_with_seeds(pairs, seeds, mode="kmer_seed")
    assert matched_ids == {"r_f0", "r_rc"}
    assert all(v == "Badnavirus_RT_Core" for v in hit_map.values())


def test_T1_2_potyviridae_rdrp_6frame_seed_baiting() -> None:
    """T1.2: 6-frame seed baiting achieves 100% recall on Potyviridae/Poacevirus RdRp reads."""
    seeds = _load_fasta_map(FIXTURES_DIR / "core_viral_seeds.faa")
    rdrp_nt = _backtranslate_peptide(seeds["Potyviridae_RdRp_Core"], use_alt_synonymous=True)
    pairs = [
        ("rdrp_1", "G" + rdrp_nt[:145], None),
        ("rdrp_2", vd_revcomp(rdrp_nt[18:168]), None),
    ]
    matched_ids, hit_map = bait_reads_with_seeds(pairs, seeds, mode="kmer_seed")
    assert matched_ids == {"rdrp_1", "rdrp_2"}
    assert hit_map["rdrp_1"] == "Potyviridae_RdRp_Core"


def test_T1_3_geminiviridae_rep_6frame_seed_baiting() -> None:
    """T1.3: 6-frame seed baiting achieves 100% recall on Geminivirus RCR Rep reads."""
    seeds = _load_fasta_map(FIXTURES_DIR / "core_viral_seeds.faa")
    rep_nt = _backtranslate_peptide(seeds["Geminiviridae_Rep_Core"], use_alt_synonymous=True)
    pairs = [
        ("rep_1", "AT" + rep_nt[:144], None),
        ("rep_2", vd_revcomp("C" + rep_nt[12:156]), None),
    ]
    matched_ids, hit_map = bait_reads_with_seeds(pairs, seeds, mode="kmer_seed")
    assert matched_ids == {"rep_1", "rep_2"}
    assert hit_map["rep_1"] == "Geminiviridae_Rep_Core"


def test_T1_4_divergent_read_rescue_in_95pct_host_background(tmp_path: pathlib.Path) -> None:
    """T1.4: Divergent viral locus (~25% nt divergence in seed) hidden in 95% host reads is rescued & assembled."""
    rng = random.Random(2026)
    host_genome = "".join(rng.choices("ACGT", k=4000))
    seeds = _load_fasta_map(FIXTURES_DIR / "core_viral_seeds.faa")

    # Mutate seed conservatively + synonymously (~25% nt divergence vs canonical codons)
    raw_pep = list(seeds["Badnavirus_RT_Core"])
    # Conservative swaps within reduced groups: I->V, K->R, D->E
    for idx, aa in enumerate(raw_pep):
        if idx % 5 == 0 and aa == "K":
            raw_pep[idx] = "R"
        elif idx % 7 == 0 and aa == "I":
            raw_pep[idx] = "V"
    divergent_seed_nt = _backtranslate_peptide("".join(raw_pep), use_alt_synonymous=True)
    flank_left = "".join(rng.choices("ACGT", k=220))
    flank_right = "".join(rng.choices("ACGT", k=220))
    viral_contig_truth = flank_left + divergent_seed_nt + flank_right

    r1_records: List[Tuple[str, str]] = []
    r2_records: List[Tuple[str, str]] = []

    # 190 host read pairs (~95%)
    for i in range(190):
        pos = rng.randint(0, len(host_genome) - 260)
        r1_records.append((f"host_{i}", host_genome[pos : pos + 110]))
        r2_records.append((f"host_{i}", vd_revcomp(host_genome[pos + 140 : pos + 250])))

    # 16 tiling viral read pairs across the 623 bp viral locus (including reverse-strand reads)
    step = 28
    v_idx = 0
    for start in range(0, len(viral_contig_truth) - 180, step):
        frag1 = viral_contig_truth[start : start + 110]
        frag2 = vd_revcomp(viral_contig_truth[start + 70 : start + 180])
        # Repeat 2x depth so de Bruijn min_edge_cov=2 is solidly supported
        for rep in range(2):
            r1_records.append((f"viral_{v_idx}_{rep}", frag1))
            r2_records.append((f"viral_{v_idx}_{rep}", frag2))
        v_idx += 1

    host_fa = tmp_path / "host.fasta"
    host_fa.write_text(f">Host_Chr1\n{host_genome}\n", encoding="utf-8")
    r1_fa = tmp_path / "r1.fasta"
    r2_fa = tmp_path / "r2.fasta"
    r1_fa.write_text("".join(f">{rid}/1\n{s}\n" for rid, s in r1_records), encoding="utf-8")
    r2_fa.write_text("".join(f">{rid}/2\n{s}\n" for rid, s in r2_records), encoding="utf-8")

    out_dir = tmp_path / "drift_out"
    rc = viro_drift_main(
        [
            "run",
            "--r1", str(r1_fa),
            "--r2", str(r2_fa),
            "--seeds", str(FIXTURES_DIR / "core_viral_seeds.faa"),
            "--host", str(host_fa),
            "-o", str(out_dir),
        ]
    )
    assert rc == 0
    summary = json.loads((out_dir / "viro_drift_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "CANDIDATES_ASSEMBLED"
    assert summary["assembled_contigs"] >= 1
    assembled = _load_fasta_map(out_dir / "assembled_contigs.fasta")
    best_len = max(len(s) for s in assembled.values())
    assert best_len >= 400


def test_T1_5_adversarial_low_complexity_polyA_and_stop_X_exclusion() -> None:
    """T1.5: Low-entropy homopolymers (`H < 1.2 bits`), `'X'`/`'*'`, and non-collinear hits yield 0 matches (`P0-VD-1`, `P0-VD-2`)."""
    seeds = _load_fasta_map(FIXTURES_DIR / "core_viral_seeds.faa")
    assert shannon_entropy_nt("A" * 24) == 0.0
    masked = mask_low_complexity_nt("A" * 60 + "ATATATATATATATATATATATAT" + "C" * 60)
    assert masked.count("N") >= 100

    # Construct adversarial reads: poly-A, dinucleotide repeat, and non-collinear swapped 5-mers
    rt_pep = seeds["Badnavirus_RT_Core"]
    non_collinear_nt = (
        _backtranslate_peptide(rt_pep[40:45])  # downstream 5-mer placed FIRST
        + "TAA"                                # in-frame stop codon interrupting peptide
        + _backtranslate_peptide(rt_pep[0:5])  # upstream 5-mer placed SECOND
    )
    adversarial_pairs = [
        ("poly_a", "A" * 150, "T" * 150),
        ("at_repeat", "ATAT" * 38, "TATA" * 38),
        ("nnnn_masked", "N" * 150, None),
        ("non_collinear_stop", non_collinear_nt + ("ACGT" * 25), None),
    ]
    matched_ids, _ = bait_reads_with_seeds(adversarial_pairs, seeds, mode="kmer_seed")
    assert len(matched_ids) == 0

    # Also verify seed index never indexes '*' or 'X'
    idx = ReducedAlphabetSeedIndex({"decoy": "XXXXX*MMMMM"})
    assert idx.match_read("NNNNNNNNNNNNNNNNNNNNNNNN") is None
    assert len(translate_6frames("ATGCGT")) == 6


def test_T1_6_host_retrotransposon_decoy_and_mate_partitioning() -> None:
    """T1.6: Host k-mer index partitions concordant/ambiguous/candidate pairs and blocks host decoys."""
    rng = random.Random(99)
    host_seq = "".join(rng.choices("ACGT", k=1200))
    viral_seq = "".join(rng.choices("TGCA", k=600))
    host_kmers = build_host_kmer_index([host_seq], k=21)

    chimeric_read = host_seq[100:160] + viral_seq[50:110]
    pairs = [
        ("concordant_host", host_seq[0:100], vd_revcomp(host_seq[200:300])),
        ("discordant_pair", host_seq[300:400], viral_seq[100:200]),
        ("chimeric_single", chimeric_read, None),
        ("pure_viral_seed", viral_seq[0:100], vd_revcomp(viral_seq[80:180])),
        ("pure_viral_ext", viral_seq[70:170], vd_revcomp(viral_seq[150:250])),
        ("host_decoy_in_pool", viral_seq[70:100] + host_seq[500:570], None),
    ]
    part = partition_host_reads(pairs, host_kmers, k=21)
    assert len(part.host_concordant) == 1
    assert {p[0] for p in part.host_spanning_ambiguous} >= {"discordant_pair", "chimeric_single"}

    rec = iterative_recruit_reads(
        all_candidate_pairs=[p for p in pairs if p[0].startswith(("pure_viral", "host_decoy"))],
        initial_seed_ids={"pure_viral_seed"},
        host_kmer_mask=host_kmers,
    )
    recruited_ids = {p[0] for p in rec.recruited_pairs}
    assert "pure_viral_ext" in recruited_ids
    assert "host_decoy_in_pool" not in recruited_ids
    assert rec.filtered_host_decoys_count >= 1


def test_T1_7_uninfected_negative_control_cli(tmp_path: pathlib.Path) -> None:
    """T1.7: Uninfected negative control exits 0 with `NO_VIRAL_CANDIDATES`; missing external binary exits 3."""
    rng = random.Random(7)
    neg_fq = tmp_path / "neg.fastq"
    lines = []
    for i in range(25):
        seq = "".join(rng.choices("ACGT", k=100))
        lines.append(f"@read_{i}\n{seq}\n+\n{'I' * 100}\n")
    neg_fq.write_text("".join(lines), encoding="utf-8")

    out_dir = tmp_path / "neg_out"
    rc = viro_drift_main(
        ["run", "--r1", str(neg_fq), "--seeds", str(FIXTURES_DIR / "core_viral_seeds.faa"), "-o", str(out_dir)]
    )
    assert rc == 0
    assert (out_dir / "assembled_contigs.fasta").read_text(encoding="utf-8") == ""
    summary = json.loads((out_dir / "viro_drift_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "NO_VIRAL_CANDIDATES"

    with patch("shutil.which", return_value=None):
        rc_missing_bin = viro_drift_main(
            [
                "run",
                "--r1", str(neg_fq),
                "--seeds", str(FIXTURES_DIR / "core_viral_seeds.faa"),
                "--mode", "diamond",
                "-o", str(tmp_path / "diamond_out"),
            ]
        )
        assert rc_missing_bin == 3

    # Verify external diamond wrapper path when mocked binary emits hit table
    dm_work = tmp_path / "dm_mock"
    dm_work.mkdir(parents=True, exist_ok=True)
    (dm_work / "diamond_hits.tsv").write_text("r_dm/1\tBadnavirus_RT_Core\t1e-12\n", encoding="utf-8")
    with patch("shutil.which", return_value="diamond"), patch("subprocess.run") as mock_run:
        rec_dm, map_dm = bait_reads_with_seeds(
            [("r_dm", "ACGT" * 30, "TGCA" * 30)],
            seeds={"Badnavirus_RT_Core": "MDDILVYSPGSKAIK"},
            mode="diamond",
            seeds_path=FIXTURES_DIR / "core_viral_seeds.faa",
            work_dir=dm_work,
        )
        assert rec_dm == {"r_dm"}
        assert map_dm["r_dm"] == "Badnavirus_RT_Core"
        assert mock_run.call_count == 2


def test_T1_8_two_phase_cycle_safe_debruijn_circular_unitig() -> None:
    """T1.8: Pure closed circular de Bruijn cycle terminates without infinite loop (`P0-VD-3`)."""
    rng = random.Random(31415)
    circle = "".join(rng.choices("ACGT", k=360))
    doubled = circle + circle
    # Tile every 10 bp around the entire circle twice so every k-mer has in_degree == out_degree == 1
    circular_pairs = []
    for rep in range(2):
        for i in range(0, len(circle), 10):
            circular_pairs.append((f"circ_{rep}_{i}", doubled[i : i + 90], None))

    contigs = assemble_viral_contigs(circular_pairs, k_values=(31,), min_contig_len=200)
    assert len(contigs) == 1
    assert contigs[0].is_circular_cycle is True
    assert len(contigs[0].sequence) >= len(circle)


# ============================================================================
# GROUP 2: Package 2 `viral-genome-tools` (10 Tests: T2.1 - T2.10)
# ============================================================================


def test_T2_1_badnavirus_2_25x_concatemer_curation() -> None:
    """T2.1: Badnavirus 2.25x concatemer trimmed to exact 2,400 bp (+/-0 bp) and rotated to `TGGTATCAGAGC`."""
    concat = _load_fasta_map(FIXTURES_DIR / "concatemer_test_panel.fasta")["Badnavirus_Concatemer_2.25x"]
    ref = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")["Badnavirus_Ref_Monomer"]
    trim = trim_concatemer_to_monomer(concat)
    assert trim.was_trimmed is True
    assert trim.monomer_length == 2400
    rot = rotate_to_landmark(trim.monomer_seq)
    assert rot.rotated_seq.startswith("TGGTATCAGAGC")
    assert rot.rotated_seq == ref


def test_T2_2_begomovirus_1_80x_concatemer_curation() -> None:
    """T2.2: Begomovirus 1.80x concatemer trimmed to exact 2,700 bp (+/-0 bp) and rotated to `TAATATTAC`."""
    concat = _load_fasta_map(FIXTURES_DIR / "concatemer_test_panel.fasta")["Begomovirus_Concatemer_1.80x"]
    ref = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")["Begomovirus_Ref_Monomer"]
    trim = trim_concatemer_to_monomer(concat)
    assert trim.was_trimmed is True
    assert trim.monomer_length == 2700
    rot = rotate_to_landmark(trim.monomer_seq)
    assert rot.rotated_seq.startswith("TAATATTAC")
    assert rot.rotated_seq == ref


def test_T2_3_circovirus_2_00x_concatemer_curation() -> None:
    """T2.3: Circovirus 2.00x concatemer trimmed to exact 1,800 bp (+/-0 bp) and rotated to `TAGTATTAC`."""
    concat = _load_fasta_map(FIXTURES_DIR / "concatemer_test_panel.fasta")["Circovirus_Concatemer_2.00x"]
    ref = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")["Circovirus_Ref_Monomer"]
    trim = trim_concatemer_to_monomer(concat)
    assert trim.was_trimmed is True
    assert trim.monomer_length == 1800
    rot = rotate_to_landmark(trim.monomer_seq)
    assert rot.rotated_seq.startswith("TAGTATTAC")
    assert rot.rotated_seq == ref


def test_T2_4_negative_strand_circular_concatemer_curation() -> None:
    """T2.4: Negative-strand (`-`) circular concatemer is flipped to `+`, trimmed to 1.0x, and rotated to Pos 1."""
    concat = _load_fasta_map(FIXTURES_DIR / "concatemer_test_panel.fasta")["Begomovirus_Concatemer_1.80x"]
    ref = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")["Begomovirus_Ref_Monomer"]
    rc_concat = vgt_revcomp(concat)
    trim = trim_concatemer_to_monomer(rc_concat)
    assert trim.monomer_length == 2700
    rot = rotate_to_landmark(trim.monomer_seq)
    assert rot.strand_flipped is True
    assert rot.rotated_seq.startswith("TAATATTAC")
    assert rot.rotated_seq == ref


def test_T2_5_indel_and_mismatch_in_dtr_overlap_tail() -> None:
    """T2.5: Revolution-2 anchor recovers exact 2,400 bp monomer despite ~1.5% error + 2 bp indel in DTR tail (`P0-VGT-1`)."""
    ref = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")["Badnavirus_Ref_Monomer"]
    rotated = ref[600:] + ref[:600]
    tail = list(rotated[:750])
    # Introduce 2 bp deletion inside overlap tail + 8 point substitutions (~1.3% divergence)
    del tail[250:252]
    for pos in (60, 130, 210, 330, 410, 520, 610, 700):
        tail[pos] = {"A": "G", "G": "A", "C": "T", "T": "C"}[tail[pos]]
    noisy_concat = rotated + "".join(tail)
    assert len(noisy_concat) == 2400 + 748

    trim = trim_concatemer_to_monomer(noisy_concat)
    assert trim.was_trimmed is True
    assert trim.monomer_length == 2400
    rot = rotate_to_landmark(trim.monomer_seq)
    assert rot.rotated_seq == ref


def test_T2_6_adversarial_linear_tir_veto(tmp_path: pathlib.Path) -> None:
    """T2.6: Linear genome with Terminal Inverted Repeats triggers TIR Veto and locks `LINEAR_CANDIDATE`."""
    tir_seq = _load_fasta_map(FIXTURES_DIR / "proven_linear_controls.fasta")["Linear_TIR_Negative_Control"]
    assert detect_terminal_inverted_repeat(tir_seq) is True
    trim = trim_concatemer_to_monomer(tir_seq)
    assert trim.has_terminal_inverted_repeat is True
    assert trim.was_trimmed is False
    assert trim.monomer_length == len(tir_seq)

    in_fa = tmp_path / "tir.fasta"
    in_fa.write_text(f">Linear_TIR_Control\n{tir_seq}\n", encoding="utf-8")
    out_dir = tmp_path / "tir_out"
    assert vgt_main(["curate", "--assembly", str(in_fa), "-o", str(out_dir)]) == 0
    summary = json.loads((out_dir / "curation_summary.json").read_text(encoding="utf-8"))
    assert summary["records"][0]["topology_decision"] == "linear_candidate"


def test_T2_7_unknown_circular_virus_orf_fallback_and_1mm_landmark() -> None:
    """T2.7: Unknown circular genome rotates to largest forward ORF (`ATG`); 1-mm motif triggers Pass-2 fallback (`P0-VGT-2`)."""
    rng = random.Random(101)
    # Build synthetic sequence free of known landmarks with one prominent 300 bp forward ORF
    filler = ("ACGT" * 200)[:700].replace("ATG", "ACG")
    orf_seq = "ATG" + ("GCTGAT" * 55) + "TAA"
    synth = filler[:300] + orf_seq + filler[300:]
    # Ensure no accidental registry landmark
    for motif in ("TGGTATCAGAGC", "TAATATTAC", "TAGTATTAC", "CAGTATTAC"):
        synth = synth.replace(motif, "CCCCCCCCC")

    rot_orf = rotate_to_landmark(synth)
    assert rot_orf.orientation_mode == "canonical_orf_orientation"
    assert rot_orf.landmark_score == 0.0
    assert rot_orf.rotated_seq.startswith("ATG")

    # 1-mismatch Pass-2 fallback on 12-bp Badnavirus PBS (TGGTATCAGAGC -> TGGTATCAGATC)
    pbs_1mm = synth[:200] + "TGGTATCAGATC" + synth[200:]
    rot_1mm = rotate_to_landmark(pbs_1mm)
    assert rot_1mm.orientation_mode == "biological_landmark_1mm"
    assert rot_1mm.landmark_score == 0.5
    assert rot_1mm.rotated_seq.startswith("TGGTATCAGATC")


def test_T2_8_analytical_read_geometry_e_junction_precision(tmp_path: pathlib.Path) -> None:
    """T2.8: `compute_e_junction_score` matches analytical closed-form equation to <1e-6 precision (`P1-VGT-3`)."""
    d_med, d_ori, n_sp, n_span, lm = 40.0, 38.0, 30, 8, 1.0
    lr, mu_ins, k_min = 150, 350.0, 20

    e_sp = d_med * (1.0 - (2.0 * k_min) / lr)
    e_span = (d_med / (2.0 * lr)) * max(0.0, mu_ins - 2.0 * lr)
    s_split = min(1.0, float(n_sp) / max(1.0, e_sp))
    s_span = min(1.0, float(n_span) / max(1.0, e_span))
    s_cov = math.exp(-abs(math.log((d_ori + 1.0) / (d_med + 1.0))))
    z = -2.20 + 2.10 * s_split + 1.60 * s_span + 1.20 * s_cov + 0.90 * lm
    expected_score = 100.0 / (1.0 + math.exp(-z))

    score, comps = compute_e_junction_score(
        n_split=n_sp,
        n_spanning_pairs=n_span,
        d_origin=d_ori,
        d_median=d_med,
        landmark_score=lm,
        read_length=lr,
        mean_insert_size=mu_ins,
        k_min_anchor=k_min,
    )
    assert abs(score - expected_score) < 1e-9
    assert abs(comps["expected_split_reads"] - expected_split_reads(d_med, lr, k_min)) < 1e-9
    assert abs(comps["expected_spanning_pairs"] - expected_spanning_pairs(d_med, lr, mu_ins)) < 1e-9

    # Also verify `profile_reads_across_seam` + CLI with seam-crossing paired reads and host check
    ref = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")["Circovirus_Ref_Monomer"]
    doubled = ref + ref
    r1_list, r2_list = [], []
    for rep in range(3):
        for pos in range(0, len(ref), 30):
            r1_list.append((f"p_{rep}_{pos}", doubled[pos : pos + 120]))
            r2_list.append((f"p_{rep}_{pos}", vgt_revcomp(doubled[pos + 180 : pos + 300])))
    r1_fa = tmp_path / "circ_r1.fasta"
    r2_fa = tmp_path / "circ_r2.fasta"
    host_fa = tmp_path / "circ_host.fasta"
    asm_fa = tmp_path / "circ_asm.fasta"
    r1_fa.write_text("".join(f">{i}\n{s}\n" for i, s in r1_list), encoding="utf-8")
    r2_fa.write_text("".join(f">{i}\n{s}\n" for i, s in r2_list), encoding="utf-8")
    host_fa.write_text(">HostRef\n" + ("GATTACA" * 80) + "\n", encoding="utf-8")
    asm_fa.write_text(f">Circovirus_Test\n{ref + ref[:400]}\n", encoding="utf-8")
    out_dir = tmp_path / "circ_curated"
    assert vgt_main(
        ["curate", "--assembly", str(asm_fa), "--r1", str(r1_fa), "--r2", str(r2_fa), "--host", str(host_fa), "-o", str(out_dir)]
    ) == 0
    rep_json = json.loads((out_dir / "curation_summary.json").read_text(encoding="utf-8"))
    assert rep_json["records"][0]["topology_decision"] == "circular_candidate"


def test_T2_9_low_depth_guard_abstention() -> None:
    """T2.9: Coverage below D_min = 10x sets S_cov = 0.0 and routes safely to `UNRESOLVED_TOPOLOGY`."""
    assert compute_coverage_uniformity(d_origin=9.0, d_median=9.0, d_min=10.0) == 0.0
    dec = classify_topology(e_junction_score=92.0, d_median=8.5, p_flank=0.0, has_tir=False)
    assert dec == TopologyDecision.UNRESOLVED_TOPOLOGY


def test_T2_10_host_flank_none_and_eve_classification() -> None:
    """T2.10: `P_flank` returns `None` on 0 junction reads, routes high-flank low-E to EVE, high-flank high-E to conflicted."""
    assert compute_host_flank_ratio(0, 0) is None
    assert classify_topology(20.0, 35.0, p_flank=0.40, has_tir=False) == TopologyDecision.INTEGRATED_CANDIDATE
    assert classify_topology(82.0, 35.0, p_flank=0.30, has_tir=False) == TopologyDecision.MIXED_OR_CONFLICTED
    assert classify_topology(18.0, 35.0, p_flank=0.02, has_tir=False) == TopologyDecision.LINEAR_CANDIDATE
    assert classify_topology(18.0, 35.0, p_flank=None, has_tir=False) == TopologyDecision.UNRESOLVED_TOPOLOGY


# ============================================================================
# GROUP 3: Package 3 `quasispecies-profiler` (10 Tests: T3.1 - T3.10)
# ============================================================================


def _fixture_variant_calls():
    sites = _parse_replicate_tsv(FIXTURES_DIR / "replicate_variants_fixture.tsv")
    calls = call_variants_across_replicates(sites, replicate_type="technical", alpha_fdr=0.01, n_boot=2000, seed=42)
    return {c.pos: c for c in calls}


def test_T3_1_spiked_isnv_4_5pct_called() -> None:
    """T3.1: True 4.5% iSNV (pos 150) is called with q < 0.01 and accurate VAF."""
    by_pos = _fixture_variant_calls()
    c = by_pos[150]
    assert c.is_called is True
    assert c.q_value < 0.01
    assert abs(c.vaf_hat - 0.045) < 0.002


def test_T3_2_spiked_isnv_2_0pct_called() -> None:
    """T3.2: True 2.0% iSNV (pos 320) is called with q < 0.01 and accurate VAF."""
    by_pos = _fixture_variant_calls()
    c = by_pos[320]
    assert c.is_called is True
    assert c.q_value < 0.01
    assert abs(c.vaf_hat - 0.020) < 0.002


def test_T3_3_spiked_isnv_1_0pct_called() -> None:
    """T3.3: True 1.0% iSNV (pos 580) is called with q < 0.01 and accurate VAF."""
    by_pos = _fixture_variant_calls()
    c = by_pos[580]
    assert c.is_called is True
    assert c.q_value < 0.01
    assert abs(c.vaf_hat - 0.010) < 0.002


def test_T3_4_spiked_isnv_0_6pct_called_and_overall_rmse() -> None:
    """T3.4: True 0.6% iSNV (pos 910) is called with q < 0.01 and overall VAF RMSE < 0.002 across all 4 true iSNVs."""
    by_pos = _fixture_variant_calls()
    c = by_pos[910]
    assert c.is_called is True
    assert c.q_value < 0.01
    expected = {150: 0.045, 320: 0.020, 580: 0.010, 910: 0.006}
    rmse = math.sqrt(sum((by_pos[p].vaf_hat - v) ** 2 for p, v in expected.items()) / len(expected))
    assert rmse < 0.002


def test_T3_5_background_pool_excludes_pcr_jackpots() -> None:
    """T3.5: Max-replicate background rule `<0.002` excludes single-replicate PCR jackpots from M_bg (`P0-QP-1`)."""
    sites = _parse_replicate_tsv(FIXTURES_DIR / "replicate_variants_fixture.tsv")
    ctx_err, global_err, _, sigma_rho = estimate_background_error_and_prior(sites, max_bg_vaf=0.002)
    # Only pos 1700 (~0.0003) and pos 1950 (~0.000233) have max replicate VAF < 0.002;
    # if jackpot [180/10000, 2/10000, 1/10000] had leaked in, global_err would exceed 0.0015!
    assert global_err <= 0.0005
    assert sigma_rho >= 3.5
    assert all(v <= 0.0005 for v in ctx_err.values())


def test_T3_6_pcr_jackpots_rejected_100pct(tmp_path: pathlib.Path) -> None:
    """T3.6: Single-replicate PCR jackpots (pos 1120, 1450) are 100% rejected (`REJECTED_PCR_JACKPOT`, `P0-QP-2`)."""
    by_pos = _fixture_variant_calls()
    for jpos in (1120, 1450):
        assert by_pos[jpos].status == "REJECTED_PCR_JACKPOT"
        assert by_pos[jpos].is_called is False
        assert by_pos[jpos].q_value == 1.0
    for bg_pos in (1700, 1950):
        assert by_pos[bg_pos].is_called is False

    out_dir = tmp_path / "qp_cli"
    rc = qp_main(
        [
            "call",
            "--manifest", str(FIXTURES_DIR / "replicate_variants_fixture.tsv"),
            "--ref", str(FIXTURES_DIR / "proven_circular_monomers.fasta"),
            "-o", str(out_dir),
        ]
    )
    assert rc == 0
    summary = json.loads((out_dir / "quasispecies_summary.json").read_text(encoding="utf-8"))
    assert summary["called_isnv_count"] == 4
    assert summary["rejected_pcr_jackpots"] == 2


def test_T3_7_ultra_deep_100000x_numerical_stability() -> None:
    """T3.7: `betaln` log-likelihood and caller produce zero NaN/inf at n = 100,000x ultra-deep coverage."""
    ll = log_beta_binomial_pmf(np.array([1500.0, 1480.0, 1520.0]), np.array([1e5, 1e5, 1e5]), 0.015, 1e-4)
    assert np.all(np.isfinite(ll))

    sites = [
        {"pos": 10, "ref": "A", "alt": "G", "k": [20, 18, 22], "n": [100000, 100000, 100000]},
        {"pos": 20, "ref": "C", "alt": "T", "k": [1500, 1490, 1510], "n": [100000, 100000, 100000]},
    ]
    res = call_variants_across_replicates(sites, n_boot=500, seed=11)
    assert all(math.isfinite(r.lrt_stat) and math.isfinite(r.vaf_hat) for r in res)
    assert res[1].is_called is True


def test_T3_8_physical_lod_gates_low_and_mid_depth() -> None:
    """T3.8: Depth <100x flags `UNASSESSED_LOW_DEPTH`; 0.5% VAF at 200x depth flags `BELOW_LOD_AT_DEPTH`."""
    sites = [
        {"pos": 1, "ref": "A", "alt": "G", "k": [2, 2, 2], "n": [60, 60, 60]},
        {"pos": 2, "ref": "C", "alt": "T", "k": [1, 1, 1], "n": [200, 200, 200]},
        {"pos": 3, "ref": "G", "alt": "A", "k": [0, 0, 0], "n": [1000, 1000, 1000]},
    ]
    res = {r.pos: r for r in call_variants_across_replicates(sites, n_boot=200)}
    assert res[1].status == "UNASSESSED_LOW_DEPTH"
    assert res[2].status == "BELOW_LOD_AT_DEPTH"
    assert res[3].status == "MONOMORPHIC_ZERO"


def test_T3_9_profile_likelihood_95_ci_scaling_and_bracketing() -> None:
    """T3.9: Sign-bracketed 95% profile likelihood CI contains true VAF and narrows from 1,000x to 10,000x (`P1-QP-4`)."""
    sites = [
        {"pos": 10, "ref": "A", "alt": "G", "k": [1, 1, 1], "n": [10000, 10000, 10000]},
        {"pos": 100, "ref": "C", "alt": "T", "k": [20, 20, 20], "n": [1000, 1000, 1000]},
        {"pos": 200, "ref": "C", "alt": "T", "k": [200, 200, 200], "n": [10000, 10000, 10000]},
    ]
    res = {r.pos: r for r in call_variants_across_replicates(sites, n_boot=400)}
    assert res[100].ci_95_low < 0.020 < res[100].ci_95_high
    assert res[200].ci_95_low < 0.020 < res[200].ci_95_high
    width_1k = res[100].ci_95_high - res[100].ci_95_low
    width_10k = res[200].ci_95_high - res[200].ci_95_low
    assert width_10k < width_1k


def test_T3_10_diversity_metrics_pi_shannon_and_nei_gojobori_dnds() -> None:
    """T3.10: Analytical check of Shannon H, nucleotide diversity pi (n<2 guard), and Nei-Gojobori dN/dS (dS=0 guard)."""
    h_val = compute_shannon_entropy([0.10], genome_length=10)
    expected_h = (-(0.10 * math.log2(0.10) + 0.90 * math.log2(0.90))) / 10.0
    assert abs(h_val - expected_h) < 1e-9

    pi_val = compute_nucleotide_diversity_pi([0.10, 0.50], [100.0, 1.0], genome_length=10)
    expected_pi = ((100.0 / 99.0) * 2.0 * 0.10 * 0.90) / 10.0  # second locus has n=1 -> safely skipped!
    assert abs(pi_val - expected_pi) < 1e-9

    cds = "ATG" + ("GCT" * 20) + "TAA"
    # Pure nonsynonymous variant (GCT=Ala -> GAT=Asp at 0-based index 4 -> dS = 0 -> dN/dS is None, zero ZeroDivisionError)
    dnds_zero_ds = compute_nei_gojobori_dnds(cds, [(4, "A", 0.20)])
    assert dnds_zero_ds["dN_dS"] is None and dnds_zero_ds["d_N"] > 0.0

    # Both synonymous (index 5: GCT->GCC) and nonsynonymous (index 4: GCT->GAT)
    dnds_both = compute_nei_gojobori_dnds(cds, [(4, "A", 0.15), (5, "C", 0.30)])
    assert dnds_both["dN_dS"] is not None and dnds_both["dN_dS"] > 0.0


# ============================================================================
# GROUP 4: Package 4 `viralbio-reporting` (8 Tests: T4.1 - T4.8)
# ============================================================================


def test_T4_1_circular_rotation_invariance_across_10_cuts_with_internal_atg() -> None:
    """T4.1: Modulo-L stop deduplication recovers 100% identical ORFs across 10 rotations slicing past internal `ATG` (`P0-VBR-1`)."""
    # Construct a 1,800 bp circular genome with:
    # - Gene A (+ strand, 120 aa) containing an internal Methionine `ATG` at aa 35
    # - Gene B (- strand, 95 aa)
    gene_a = "ATG" + ("GCTGAC" * 17) + "ATG" + ("CGTTCT" * 41) + "TAA"  # 357 bp (118 aa)
    gene_b_plus = "ATG" + ("CCAGAA" * 46) + "TAG"                       # 282 bp (93 aa)
    spacer_1 = "TAATAG" * 90                                            # stop-dense noncoding spacer (540 bp)
    spacer_2 = ("TGATAA" * 104)[:621]                                   # stop-dense noncoding spacer (621 bp)
    genome = gene_a + spacer_1 + vbr_revcomp(gene_b_plus) + spacer_2
    assert len(genome) == 1800

    baseline_orfs = annotate_genome_orfs(genome, topology="circular", min_aa=60)
    baseline_proteins = sorted((f.strand, f.aa_sequence) for f in baseline_orfs)
    assert len(baseline_proteins) >= 2

    # Cut at 10 positions around the circle, explicitly including cuts between Gene A's canonical ATG (0) and internal ATG (105)
    cut_offsets = [0, 36, 75, 105, 210, 480, 777, 1020, 1350, 1695]
    for cut in cut_offsets:
        rotated = genome[cut:] + genome[:cut]
        rot_orfs = annotate_genome_orfs(rotated, topology="circular", min_aa=60)
        rot_proteins = sorted((f.strand, f.aa_sequence) for f in rot_orfs)
        assert rot_proteins == baseline_proteins, f"Rotation variance at cut={cut}!"
        for feat in rot_orfs:
            ok, msg = validate_cds_feature(rotated, feat)
            assert ok, f"Invalid CDS at cut={cut}: {msg}"


def test_T4_2_single_revolution_deduplication_no_2x_duplicates() -> None:
    """T4.2: Annotation on circular monomers emits zero duplicate ORFs from Revolution 2."""
    monomers = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")
    for name, seq in monomers.items():
        orfs = annotate_genome_orfs(seq, topology="circular", min_aa=50)
        # Every (strand, stop_key_mod_L) must appear at most once
        stop_endpoints = [(f.strand, f.intervals_1based[-1][1]) for f in orfs]
        assert len(stop_endpoints) == len(set(stop_endpoints)), f"Duplicate stop endpoint in {name}"


def test_T4_3_negative_strand_origin_spanning_orf_projection() -> None:
    """T4.3: Negative-strand (`-`) origin-spanning ORF projects to valid descending compound intervals (`P0-VBR-2`)."""
    plus_cds = "ATG" + ("GCTTTC" * 50) + "TGA"  # 306 bp (101 aa)
    rc_cds = vbr_revcomp(plus_cds)
    spacer = "TAATGA" * 150  # 900 bp stop-dense spacer
    # Place half of rc_cds at the 3' end of genome and half at 5' start so the (-) ORF crosses the origin!
    split = 150
    circular_genome = rc_cds[split:] + spacer + rc_cds[:split]
    L = len(circular_genome)

    orfs = annotate_genome_orfs(circular_genome, topology="circular", min_aa=80)
    minus_wrapping = [f for f in orfs if f.strand == "-" and f.spans_origin]
    assert len(minus_wrapping) == 1
    feat = minus_wrapping[0]
    assert len(feat.intervals_1based) == 2
    (a1, b1), (a2, b2) = feat.intervals_1based
    assert a1 > b1 == 1
    assert a2 == L and a2 > b2
    ok, msg = validate_cds_feature(circular_genome, feat)
    assert ok, msg
    assert reconstruct_cds_nt_from_intervals(circular_genome, feat) == plus_cds


def test_T4_4_stop_free_circular_frame_capped() -> None:
    """T4.4: Stop-free circular reading frame (>1.0x revolution) is strictly capped/filtered (`nt_length <= L`)."""
    stop_free_circle = "ATG" + ("GCT" * 299)  # 900 bp with zero in-frame stop codons across infinite revolutions
    orfs = annotate_genome_orfs(stop_free_circle, topology="circular", min_aa=50)
    assert all(f.nt_length <= len(stop_free_circle) for f in orfs)
    assert len(orfs) == 0  # Zero genuine stop-terminated single-revolution ORFs


def test_T4_5_partial_linear_cds_codon_start_frames() -> None:
    """T4.5: 5'-partial linear CDS anchors interval at `<1` with `codon_start = f + 1` and passes validator (`P1-VBR-3`)."""
    body = "GCTGAT" * 45 + "TAA"  # 273 bp coding body ending at stop codon, lacking upstream ATG
    for frame_offset, prefix in ((0, ""), (1, "C"), (2, "CC")):
        lin_seq = prefix + body + ("TAATAG" * 20)
        orfs = annotate_genome_orfs(lin_seq, topology="linear", min_aa=60, include_partial_linear=True)
        partials = [f for f in orfs if f.partial_5prime and f.strand == "+"]
        assert len(partials) == 1
        p = partials[0]
        assert p.intervals_1based[0][0] == 1
        assert p.codon_start == frame_offset + 1
        ok, msg = validate_cds_feature(lin_seq, p)
        assert ok, msg
        tbl = format_ncbi_feature_table("LinTest", [p])
        assert tbl.splitlines()[1].startswith("<1\t")
        assert f"\t\t\tcodon_start\t{frame_offset + 1}" in tbl


def test_T4_6_ncbi_tbl_generation_on_proven_circular_monomers(tmp_path: pathlib.Path) -> None:
    """T4.6: `viralbio-reporting annotate` generates valid `.tbl` and `.fsa` packages on proven monomers."""
    out_circ = tmp_path / "annot_circ"
    rc1 = vbr_main(
        [
            "annotate",
            "--fasta", str(FIXTURES_DIR / "proven_circular_monomers.fasta"),
            "--topology", "circular",
            "--min-aa", "40",
            "-o", str(out_circ),
        ]
    )
    assert rc1 == 0
    summary = json.loads((out_circ / "annotation_summary.json").read_text(encoding="utf-8"))
    assert len(summary["records"]) == 3
    assert all(r["tbl_syntax_passed"] is True for r in summary["records"])


def test_T4_7_internal_tbl_syntax_validator_catches_malformed_tables() -> None:
    """T4.7: `validate_tbl_syntax` accepts valid compound `.tbl` files and rejects malformed headers/qualifiers."""
    bad_empty, errs_empty = validate_tbl_syntax("")
    assert bad_empty is False and len(errs_empty) >= 1

    malformed = ">Feature TestSeq\njoin(1..100)\t100\tCDS\nabc\t0\tCDS\n\tbad_qual\n"
    ok_mal, errs_mal = validate_tbl_syntax(malformed)
    assert ok_mal is False
    assert len(errs_mal) >= 3


def test_T4_8_sha256_provenance_chain_integrity(tmp_path: pathlib.Path) -> None:
    """T4.8: `SequenceProvenanceChain` verifies intact SHA-256 chain and detects tampering."""
    chain = SequenceProvenanceChain(tmp_path / "virus_prov.json")
    chain.append("raw", "a" * 64, "b" * 64, {"note": "trimmed DTR"})
    chain.append("rotate", "b" * 64, "c" * 64, {"note": "rotated to Pos 1"})
    assert chain.verify_integrity() is True

    # Tamper with previous_hash link
    chain.records[1].previous_hash = "0" * 64
    assert chain.verify_integrity() is False


# ============================================================================
# GROUP 5: Package 5 `viro-networks` & Orchestrator (8 Tests: T5.1 - T5.8)
# ============================================================================


def test_T5_1_planted_cooccurrence_recovery_on_cohort_fixture() -> None:
    """T5.1: Positive-floored SparCC + graphical model recovers planted `Poacevirus_SCSMV` <-> `Polerovirus_SCYLV` edge (`P0-VN-2`)."""
    df = pd.read_csv(FIXTURES_DIR / "virome_cohort_counts.tsv", sep="\t", index_col=0)
    mat = df.to_numpy(dtype=float)
    taxa = list(df.columns)

    sparcc = compute_sparcc_matrix(mat, pseudocount=0.5)
    assert np.all(np.isfinite(sparcc))

    edges, diag = infer_consensus_network(mat, taxa, n_boot=100, stability_threshold=0.80, seed=42)
    pairs = {tuple(sorted((e.taxon_a, e.taxon_b))): e for e in edges}
    target = ("Poacevirus_SCSMV", "Polerovirus_SCYLV")
    assert target in pairs
    assert pairs[target].sign == 1
    assert pairs[target].bootstrap_stability >= 0.80
    assert diag["consensus_edge_count"] >= 1


def test_T5_2_pseudocount_sensitivity_concordance() -> None:
    """T5.2: CLR transformation sums to 0 per sample and reports pseudocount sensitivity across {0.1, 0.5, 1.0}."""
    df = pd.read_csv(FIXTURES_DIR / "virome_cohort_counts.tsv", sep="\t", index_col=0)
    mat = df.to_numpy(dtype=float)
    clr = clr_transform(mat, pseudocount=0.5)
    assert np.allclose(np.sum(clr, axis=1), 0.0, atol=1e-10)

    edges, diag = infer_consensus_network(mat, list(df.columns), pseudocount_grid=(0.1, 0.5, 1.0), n_boot=60)
    assert 0.0 <= diag["pseudocount_extremes_jaccard"] <= 1.0
    assert all(0.0 <= e.pseudocount_concordance <= 1.0 for e in edges)


def test_T5_3_sloan_neutral_model_and_zero_variance_r2_guard() -> None:
    """T5.3: Sloan neutral model converges cleanly on cohort fixture AND constant-occupancy edge case (`P0-VN-1`)."""
    df = pd.read_csv(FIXTURES_DIR / "virome_cohort_counts.tsv", sep="\t", index_col=0)
    rows, diag = fit_sloan_neutral_model(df.to_numpy(dtype=float), list(df.columns))
    assert len(rows) == 6
    assert math.isfinite(diag["Ne_m"]) and diag["Ne_m"] > 0.0
    assert math.isfinite(diag["r_squared"])
    assert all(0.0 <= r.wilson_95_low <= r.wilson_95_high <= 1.0 for r in rows)

    # Adversarial constant-occupancy matrix (all taxa present in 100% of samples -> SS_tot == 0)
    const_mat = np.full((12, 4), 100.0) + np.arange(48, dtype=float).reshape(12, 4)
    _, const_diag = fit_sloan_neutral_model(const_mat, ["V1", "V2", "V3", "V4"])
    assert math.isfinite(const_diag["r_squared"]) and 0.0 <= const_diag["r_squared"] <= 1.0


def test_T5_4_adversarial_ultra_small_cohort_abstains_exit_4(tmp_path: pathlib.Path) -> None:
    """T5.4: Ultra-small cohort (N=3 < 5) raises `CohortInsufficientError` and CLI returns exit code 4."""
    tiny_tsv = tmp_path / "tiny_counts.tsv"
    tiny_tsv.write_text(
        "sample_id\tV1\tV2\tV3\nS1\t10\t0\t5\nS2\t0\t20\t5\nS3\t15\t10\t0\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "tiny_out"
    rc = vn_main(["infer", "--counts", str(tiny_tsv), "-o", str(out_dir)])
    assert rc == 4
    summary = json.loads((out_dir / "network_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "ABSTAIN_INSUFFICIENT_SAMPLES"


def test_T5_5_pure_python_hub_perturbation_without_networkx() -> None:
    """T5.5: Pure-Python Brandes betweenness and eigenvalue natural connectivity match analytical star graph (`P1-VN-3`)."""
    taxa = ["Hub", "Leaf1", "Leaf2", "Leaf3", "Leaf4"]
    edges = [("Hub", "Leaf1"), ("Hub", "Leaf2"), ("Hub", "Leaf3"), ("Hub", "Leaf4")]
    records = evaluate_hub_perturbation(taxa, edges)
    top = records[0]
    assert top.taxon == "Hub"
    assert top.degree == 4
    assert abs(top.degree_centrality - 1.0) < 1e-9
    # On a 5-node star, Hub lies on all 6 shortest paths between the 4 leaves -> normalized betweenness == 1.0
    assert abs(top.betweenness_centrality - 1.0) < 1e-9
    assert top.natural_connectivity_drop_pct > 50.0
    assert abs(top.lcc_fraction_after_removal - 0.25) < 1e-9


def test_T5_6_all_five_standalone_clis_version_and_error_codes(tmp_path: pathlib.Path) -> None:
    """T5.6: All 5 standalone Choice B CLIs handle `--version`, missing input (exit 2), and full cohort run (exit 0)."""
    for cli_fn in (viro_drift_main, vgt_main, qp_main, vbr_main, vn_main):
        assert cli_fn(["--version"]) == 0

    missing = str(tmp_path / "does_not_exist.file")
    out = str(tmp_path / "err_out")
    assert viro_drift_main(["run", "--r1", missing, "--seeds", missing, "-o", out]) == 2
    assert vgt_main(["curate", "--assembly", missing, "-o", out]) == 2
    assert qp_main(["call", "--manifest", missing, "-o", out]) == 2
    assert vbr_main(["annotate", "--fasta", missing, "-o", out]) == 2
    assert vn_main(["infer", "--counts", missing, "-o", out]) == 2

    # Full `viro-networks infer` CLI run on the 24-sample cohort fixture
    vn_out = tmp_path / "vn_full"
    assert vn_main(["infer", "--counts", str(FIXTURES_DIR / "virome_cohort_counts.tsv"), "-o", str(vn_out)]) == 0
    assert (vn_out / "consensus_edges.tsv").exists()
    assert (vn_out / "sloan_neutral_fit.tsv").exists()
    assert (vn_out / "topological_hub_candidates.tsv").exists()


def test_T5_7_master_orchestrator_end_to_end_pipeline(tmp_path: pathlib.Path) -> None:
    """T5.7: Master `virome-suite run` orchestrates curation, annotation, quasispecies, and networks end-to-end."""
    assert orchestrator_main(["--version"]) == 0
    assert orchestrator_main(["check-env", "--json"]) == 0
    assert orchestrator_main(["run", "--r1", str(tmp_path / "missing.fa"), "-o", str(tmp_path / "x")]) == 2

    e2e_out = tmp_path / "e2e_suite"
    rc = orchestrator_main(
        [
            "run",
            "--r1", str(FIXTURES_DIR / "concatemer_test_panel.fasta"),
            "--replicates", str(FIXTURES_DIR / "replicate_variants_fixture.tsv"),
            "--counts", str(FIXTURES_DIR / "virome_cohort_counts.tsv"),
            "-o", str(e2e_out),
        ]
    )
    assert rc == 0
    summary = json.loads((e2e_out / "pipeline_summary.json").read_text(encoding="utf-8"))
    assert summary["stages"]["viral_genome_tools_exit"] == 0
    assert summary["stages"]["annotated_genomes"] == 3
    assert summary["stages"]["quasispecies_profiler_exit"] == 0
    assert summary["stages"]["viro_networks_exit"] == 0
    assert (e2e_out / "pipeline_provenance.json").exists()

    # Also exercise raw-reads + `--seeds` + `--host` branch in `orchestrator.runner`
    seeds = _load_fasta_map(FIXTURES_DIR / "core_viral_seeds.faa")
    rt_nt = _backtranslate_peptide(seeds["Badnavirus_RT_Core"])
    locus = ("ACGTGCAT" * 30) + rt_nt + ("TGCATGCA" * 30)
    r1_s = tmp_path / "seed_r1.fasta"
    r2_s = tmp_path / "seed_r2.fasta"
    host_s = tmp_path / "seed_host.fasta"
    r1_s.write_text(
        "".join(f">s_{i}\n{locus[p : p + 120]}\n" for i, p in enumerate(range(0, len(locus) - 180, 20)) for _ in range(2)),
        encoding="utf-8",
    )
    r2_s.write_text(
        "".join(
            f">s_{i}\n{vd_revcomp(locus[p + 60 : p + 180])}\n"
            for i, p in enumerate(range(0, len(locus) - 180, 20))
            for _ in range(2)
        ),
        encoding="utf-8",
    )
    host_s.write_text(">H1\n" + ("CCGGAATT" * 80) + "\n", encoding="utf-8")
    e2e_reads_out = tmp_path / "e2e_reads"
    assert orchestrator_main(
        [
            "run",
            "--r1", str(r1_s),
            "--r2", str(r2_s),
            "--seeds", str(FIXTURES_DIR / "core_viral_seeds.faa"),
            "--host", str(host_s),
            "-o", str(e2e_reads_out),
        ]
    ) == 0


def test_T5_8_orchestrator_audit_inventory_classifier_and_runner(tmp_path: pathlib.Path) -> None:
    """T5.8: `orchestrator.audit_inventory` accurately classifies legacy script categories and writes manifest TSV."""
    scbv = tmp_path / "r-pkg_tst" / "scripts" / "SCBV_all"
    virome = tmp_path / "r-pkg_tst" / "scripts" / "VIROME_all"
    workflow = tmp_path / "r-pkg_tst" / "workflow"
    (scbv / "Systembiology_databases" / "interproscan").mkdir(parents=True)
    (scbv / "08_secondary_solution" / "plan_d" / "finishing").mkdir(parents=True)
    (scbv / "backup_20260428").mkdir(parents=True)
    virome.mkdir(parents=True)
    workflow.mkdir(parents=True)

    sample_files = {
        scbv / "Systembiology_databases" / "interproscan" / "helper.py": "import os\n",
        scbv / "08_secondary_solution" / "plan_d" / "finishing" / "run_spades.sh": "spades.py -k 21\n",
        scbv / "backup_20260428" / "old.sh": "bowtie2\n",
        scbv / "08D_01_seed_mining.py": "# diamond TGGTATCAGAGC /mnt/c\n",
        scbv / "04_canonical_monomer_curator.py": "# TGGTATCAGAGC\n",
        scbv / "05_variant_calling_quasispecies.py": "# samtools lofreq\n",
        scbv / "07_generate_ncbi_submission.py": "# table2asn\n",
        virome / "sb_07_consensus_network.py": "# sparcc\n",
        virome / "11_totivirus_extraction.sh": "# fungal\n",
        virome / "00_setup_workstation_env.sh": "# setup\n",
        virome / "misc_monitor.sh": "echo status\n",
        workflow / "dup_copy.sh": "echo status\n",
    }
    for fpath, content in sample_files.items():
        fpath.write_text(content, encoding="utf-8")

    audit_inventory.main(root=tmp_path)
    out_manifest = tmp_path / "scripts_manifest.tsv"
    assert out_manifest.exists()
    with out_manifest.open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 11
    pkgs = {r["target_package"] for r in rows}
    assert {
        "viro_drift",
        "viral_genome_tools",
        "quasispecies_profiler",
        "viralbio_reporting",
        "viro_networks",
        "orchestrator",
    }.issubset(pkgs)


def test_T5_9_multi_sample_multi_replicate_cohort_runner(tmp_path: pathlib.Path) -> None:
    """T5.9: `virome-suite run-cohort` processes multi-sample + multi-replicate sample sheets end-to-end."""
    # Missing or invalid samplesheet returns exit code 2
    assert orchestrator_main(["run-cohort", "--samplesheet", str(tmp_path / "missing.tsv"), "-o", str(tmp_path / "o")]) == 2
    bad_sheet = tmp_path / "bad_sheet.tsv"
    bad_sheet.write_text("wrong_col\nval\n", encoding="utf-8")
    assert orchestrator_main(["run-cohort", "--samplesheet", str(bad_sheet), "-o", str(tmp_path / "o")]) == 2

    monomers = _load_fasta_map(FIXTURES_DIR / "proven_circular_monomers.fasta")
    badna = monomers["Badnavirus_Ref_Monomer"]
    begomo = monomers["Begomovirus_Ref_Monomer"]
    circo = monomers["Circovirus_Ref_Monomer"]

    # Create a 6-sample cohort (S01 has 3 technical replicates with a true 5% variant at index 120; S02..S06 have 1-2 replicates)
    data_dir = tmp_path / "cohort_inputs"
    data_dir.mkdir(parents=True, exist_ok=True)
    sheet_lines = ["sample_id\treplicate_id\tr1_path\tr2_path\tgroup"]

    # Sample S01: 3 replicates tiling a 450 bp slice of Badnavirus with a consistent 20% A/G polymorphism at offset 90
    s01_ref = badna[:450]
    mut_base = "G" if s01_ref[90] != "G" else "A"
    s01_alt = s01_ref[:90] + mut_base + s01_ref[91:]
    for rep_idx in range(1, 4):
        rep_fa = data_dir / f"S01_Rep{rep_idx}.fasta"
        seqs = [f">Badna_contig\n{badna + badna[:200]}\n"]
        # Add 25 wild-type tiling reads + 6 mutant reads across the region so depth >= 25 and alt ~19%
        for c_i in range(25):
            seqs.append(f">wt_{rep_idx}_{c_i}\n{s01_ref[30:180]}\n")
        for m_i in range(6):
            seqs.append(f">mut_{rep_idx}_{m_i}\n{s01_alt[30:180]}\n")
        rep_fa.write_text("".join(seqs), encoding="utf-8")
        sheet_lines.append(f"Sample_01\tRep{rep_idx}\t{rep_fa.name}\tNA\tSymptomatic")

    # Samples S02..S06: distinct combinations of Badnavirus, Begomovirus, and Circovirus across 5 more samples (total N=6 >= 5)
    sample_specs = [
        ("Sample_02", [badna, begomo], "Symptomatic"),
        ("Sample_03", [begomo, circo], "Asymptomatic"),
        ("Sample_04", [badna, begomo], "Symptomatic"),
        ("Sample_05", [circo], "Asymptomatic"),
        ("Sample_06", [badna, begomo, circo], "Symptomatic"),
    ]
    for sid, virus_list, grp in sample_specs:
        fa_p = data_dir / f"{sid}_Rep1.fasta"
        fa_p.write_text(
            "".join(f">{sid}_v{idx}\n{vseq + vseq[:150]}\n" for idx, vseq in enumerate(virus_list, start=1)),
            encoding="utf-8",
        )
        sheet_lines.append(f"{sid}\tRep1\t{fa_p.name}\t\t{grp}")

    sheet_tsv = data_dir / "study_samplesheet.tsv"
    sheet_tsv.write_text("\n".join(sheet_lines) + "\n", encoding="utf-8")

    cohort_out = tmp_path / "cohort_results"
    rc = orchestrator_main(
        [
            "run-cohort",
            "--samplesheet", str(sheet_tsv),
            "--bootstraps", "40",
            "--stability", "0.60",
            "-o", str(cohort_out),
        ]
    )
    assert rc == 0
    csum = json.loads((cohort_out / "cohort_summary.json").read_text(encoding="utf-8"))
    assert csum["n_samples"] == 6
    assert csum["total_replicate_libraries"] == 8
    assert csum["catalog_viral_taxa"] >= 3
    assert csum["viro_networks_status"] == "SUCCESS"
    assert csum["samples"]["Sample_01"]["quasispecies_profiler_exit"] == 0
    assert (cohort_out / "04_cohort_matrices" / "cohort_abundance_counts.tsv").exists()
    assert (cohort_out / "04_cohort_matrices" / "cohort_metadata.tsv").exists()
    assert (cohort_out / "05_cohort_ecological_networks" / "sloan_neutral_fit.tsv").exists()


if __name__ == "__main__":
    unittest.main()
