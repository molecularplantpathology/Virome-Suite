"""Unit & Behavioral Acceptance Tests for Completion Mark 1 (Gate 1 Foundation).

Executable directly via standard-library `python -m unittest discover -s tests/unit -v`
or `pytest`.
"""

from __future__ import annotations

import csv
import pathlib
import sys
import tempfile
import unittest

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from orchestrator.branching import (
    build_contig_evidence_record,
    classify_topology_route,
)
from orchestrator.cli import main as orchestrator_main
from orchestrator.preflight import run_preflight_check
from orchestrator.schemas import (
    MULTI_FAMILY_LANDMARK_REGISTRY,
    AppendOnlyProvenanceLog,
    TopologyState,
    compute_coverage_uniformity,
    compute_host_flank_ratio,
    expected_spanning_pairs,
    expected_split_reads,
)


class TestMark1Foundation(unittest.TestCase):
    """Behavioral verification suite for Completion Mark 1 (Gate 1)."""

    def test_scripts_manifest_inventory_complete(self) -> None:
        """Deliverable D1: scripts_manifest.tsv audits all 263 scripts."""
        manifest = ROOT_DIR / "scripts_manifest.tsv"
        self.assertTrue(manifest.exists(), "scripts_manifest.tsv must exist")
        with manifest.open("r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 263, f"Expected 263 scripts in manifest, found {len(rows)}")

        required_cols = {
            "script_id",
            "source_collection",
            "relative_path",
            "sha256",
            "language",
            "status",
            "license_provenance",
            "target_package",
            "external_tools",
            "generalization_action",
        }
        self.assertTrue(required_cols.issubset(set(rows[0].keys())))

    def test_multi_virus_fixtures_present(self) -> None:
        """Deliverable D2: Multi-family proven viral test fixtures exist."""
        fixtures_dir = ROOT_DIR / "tests" / "fixtures"
        expected_files = [
            "proven_circular_monomers.fasta",
            "concatemer_test_panel.fasta",
            "proven_linear_controls.fasta",
            "core_viral_seeds.faa",
            "replicate_variants_fixture.tsv",
            "virome_cohort_counts.tsv",
        ]
        for fname in expected_files:
            fpath = fixtures_dir / fname
            self.assertTrue(fpath.exists() and fpath.stat().st_size > 0, f"Missing fixture: {fname}")

    def test_benchmark_protocol_locked(self) -> None:
        """Deliverable D3: Pre-registered benchmark protocol YAML is valid and complete."""
        proto_path = ROOT_DIR / "benchmarks" / "benchmark_protocol.yaml"
        self.assertTrue(proto_path.exists())
        text = proto_path.read_text(encoding="utf-8")
        self.assertIn('protocol_version: "1.0.0"', text)
        self.assertIn("BM-1_divergent_viral_rescue:", text)
        self.assertIn("BM-6_end_to_end_workflow_integration:", text)

    def test_read_geometry_scaling_non_binary(self) -> None:
        """Round 4/5 mathematical fix: expected split and spanning counts scale continuously."""
        exp_split = expected_split_reads(d_median=30.0, read_length=150, k_min_anchor=20)
        self.assertAlmostEqual(exp_split, 22.0, places=6)

        exp_span = expected_spanning_pairs(d_median=30.0, read_length=150, mean_insert_size=350.0)
        self.assertAlmostEqual(exp_span, 5.0, places=6)

        # Overlapping reads (mu_ins <= 2 * L_r) never produce negative expectations
        self.assertEqual(
            expected_spanning_pairs(d_median=50.0, read_length=150, mean_insert_size=250.0),
            0.0,
        )

    def test_coverage_uniformity_and_host_flank_guards(self) -> None:
        """Verify low-depth guard (D_median < 10x -> 0.0) and unassessed P_flank -> None."""
        self.assertEqual(compute_coverage_uniformity(d_origin=0.0, d_median=0.0, d_min=10.0), 0.0)
        self.assertEqual(compute_coverage_uniformity(d_origin=8.0, d_median=8.0, d_min=10.0), 0.0)
        self.assertAlmostEqual(
            compute_coverage_uniformity(d_origin=50.0, d_median=50.0, d_min=10.0),
            1.0,
            places=9,
        )

        # Zero informative reads must return None (Unassessed), never 0.0
        self.assertIsNone(compute_host_flank_ratio(n_split=0, n_host_viral_chimera=0))
        self.assertAlmostEqual(
            compute_host_flank_ratio(n_split=18, n_host_viral_chimera=2) or 0.0,
            0.10,
            places=9,
        )

    def test_topology_routing_precedence(self) -> None:
        """Verify branching router distinguishes circular, linear, integrated, conflicted, and abstain."""
        rec_circ = build_contig_evidence_record(
            contig_id="Badnavirus_Candidate",
            length_bp=2400,
            d_median=40.0,
            d_origin=38.0,
            n_split_reads=30,
            n_spanning_pairs=8,
            n_host_viral_chimeras=0,
            landmark_score=1.0,
            landmark_name="tRNA-Met_PBS",
            raw_sha256="aaa",
            curated_sha256="bbb",
        )
        self.assertEqual(rec_circ.topology_state, TopologyState.CIRCULAR_CANDIDATE)
        self.assertGreaterEqual(rec_circ.e_junction_score, 70.0)

        state_low = classify_topology_route(e_junction_score=85.0, d_median=6.0, p_flank=0.0)
        self.assertEqual(state_low, TopologyState.UNRESOLVED_TOPOLOGY)

        state_conflict = classify_topology_route(e_junction_score=80.0, d_median=35.0, p_flank=0.25)
        self.assertEqual(state_conflict, TopologyState.MIXED_OR_CONFLICTED)

    def test_multi_family_landmark_registry(self) -> None:
        """Verify landmark registry spans Badnavirus, Caulimovirus, Geminiviridae, Circoviridae, Nanoviridae."""
        families = {m.family for m in MULTI_FAMILY_LANDMARK_REGISTRY}
        self.assertTrue(
            {"Caulimoviridae", "Geminiviridae", "Circoviridae", "Nanoviridae"}.issubset(families)
        )

    def test_append_only_provenance_chain(self) -> None:
        """Verify hash-chained append-only provenance manifest."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = pathlib.Path(tmpdir) / "provenance.json"
            plog = AppendOnlyProvenanceLog(log_file)
            e1 = plog.append_event("dtr_trim", "viral-genome-tools", "1.0.0", "hash_in", "hash_mid")
            e2 = plog.append_event(
                "landmark_rotate", "viral-genome-tools", "1.0.0", "hash_mid", "hash_out"
            )
            self.assertEqual(e2.previous_event_hash, e1.event_hash())
            self.assertTrue(log_file.exists())

    def test_cli_preflight_and_audit_summary(self) -> None:
        """Deliverable D4: `virome-suite check-env` and `audit-summary` execute cleanly."""
        report = run_preflight_check(ROOT_DIR)
        self.assertTrue(report["core_python_ready"])
        self.assertTrue(report["workspace_artifacts"]["scripts_manifest_present"])
        self.assertEqual(orchestrator_main(["check-env"]), 0)
        self.assertEqual(
            orchestrator_main(["audit-summary", "--manifest", str(ROOT_DIR / "scripts_manifest.tsv")]),
            0,
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
