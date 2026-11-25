"""End-to-End Branching Execution Engine (`orchestrator.runner`)."""

from __future__ import annotations

import hashlib
import json
import pathlib
from typing import Any, Dict, Optional

from Bio import SeqIO

from orchestrator import __version__
from orchestrator.schemas import AppendOnlyProvenanceLog


def execute_end_to_end_pipeline(
    r1_path: pathlib.Path,
    output_dir: pathlib.Path,
    r2_path: Optional[pathlib.Path] = None,
    seeds_path: Optional[pathlib.Path] = None,
    host_path: Optional[pathlib.Path] = None,
    replicate_manifest: Optional[pathlib.Path] = None,
    cohort_counts: Optional[pathlib.Path] = None,
) -> Dict[str, Any]:
    """Orchestrate `viro_drift` -> morphology branching (`viral_genome_tools`) -> `viralbio_reporting` (+ optional variants & networks)."""
    from viro_drift.cli import main as drift_main
    from viral_genome_tools.cli import main as vgt_main
    from viralbio_reporting.cli import main as vbr_main

    output_dir.mkdir(parents=True, exist_ok=True)
    prov = AppendOnlyProvenanceLog(output_dir / "pipeline_provenance.json")
    r1_hash = hashlib.sha256(r1_path.read_bytes()).hexdigest()[:16]

    stage_status: Dict[str, Any] = {}

    # Stage 1: Divergent Discovery (`viro-drift`)
    drift_dir = output_dir / "01_viro_drift"
    if seeds_path is not None and seeds_path.exists():
        drift_argv = ["run", "--r1", str(r1_path), "--seeds", str(seeds_path), "-o", str(drift_dir)]
        if r2_path is not None:
            drift_argv.extend(["--r2", str(r2_path)])
        if host_path is not None:
            drift_argv.extend(["--host", str(host_path)])
        rc_drift = drift_main(drift_argv)
        candidate_fasta = drift_dir / "assembled_contigs.fasta"
        stage_status["viro_drift_exit"] = rc_drift
    else:
        # Treat r1_path directly as assembled contigs FASTA when no seeds are supplied
        candidate_fasta = r1_path
        stage_status["viro_drift_exit"] = "SKIPPED_DIRECT_ASSEMBLY_INPUT"

    cand_hash = (
        hashlib.sha256(candidate_fasta.read_bytes()).hexdigest()[:16]
        if candidate_fasta.exists()
        else "EMPTY"
    )
    prov.append_event("stage1_discovery", "viro-drift", __version__, r1_hash, cand_hash)

    # Stage 2: Morphology Branching & Circular Curation (`viral-genome-tools`)
    curation_dir = output_dir / "02_genome_curation"
    curated_records = []
    if candidate_fasta.exists() and candidate_fasta.stat().st_size > 0:
        vgt_argv = ["curate", "--assembly", str(candidate_fasta), "-o", str(curation_dir)]
        if seeds_path is not None and r1_path.exists():
            vgt_argv.extend(["--r1", str(r1_path)])
            if r2_path is not None:
                vgt_argv.extend(["--r2", str(r2_path)])
        if host_path is not None:
            vgt_argv.extend(["--host", str(host_path)])
        rc_vgt = vgt_main(vgt_argv)
        stage_status["viral_genome_tools_exit"] = rc_vgt
        cur_fa = curation_dir / "curated_genomes.fasta"
        if cur_fa.exists():
            curated_records = list(SeqIO.parse(str(cur_fa), "fasta"))
            cur_hash = hashlib.sha256(cur_fa.read_bytes()).hexdigest()[:16]
            prov.append_event("stage2_curation", "viral-genome-tools", __version__, cand_hash, cur_hash)

    # Stage 3: Functional Branching -> Seam-Aware Annotation (`viralbio-reporting`)
    annot_dir = output_dir / "03_annotation_genbank"
    annotated_count = 0
    for rec in curated_records:
        desc = rec.description
        topo = "linear" if "topology=linear_candidate" in desc else "circular"
        rec_fa = annot_dir / f"{rec.id}_input.fasta"
        rec_fa.parent.mkdir(parents=True, exist_ok=True)
        rec_fa.write_text(f">{rec.id}\n{str(rec.seq)}\n", encoding="utf-8")
        rc_vbr = vbr_main(
            ["annotate", "--fasta", str(rec_fa), "--topology", topo, "--min-aa", "40", "-o", str(annot_dir / rec.id)]
        )
        if rc_vbr == 0:
            annotated_count += 1
    stage_status["annotated_genomes"] = annotated_count

    # Stage 4 (Optional): Replicate Intra-Host Variant Profiling (`quasispecies-profiler`)
    if replicate_manifest is not None and replicate_manifest.exists():
        from quasispecies_profiler.cli import main as qp_main

        var_dir = output_dir / "04_quasispecies_variants"
        stage_status["quasispecies_profiler_exit"] = qp_main(
            ["call", "--manifest", str(replicate_manifest), "-o", str(var_dir)]
        )

    # Stage 5 (Optional): Cohort Co-Occurrence Network & Sloan Model (`viro-networks`)
    if cohort_counts is not None and cohort_counts.exists():
        from viro_networks.cli import main as vn_main

        net_dir = output_dir / "05_ecological_networks"
        stage_status["viro_networks_exit"] = vn_main(
            ["infer", "--counts", str(cohort_counts), "-o", str(net_dir)]
        )

    summary = {
        "orchestrator_version": __version__,
        "stages": stage_status,
        "provenance_log": str(output_dir / "pipeline_provenance.json"),
    }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
