"""End-to-End Single-Sample & Multi-Sample Cohort Execution Engine (`orchestrator.runner`)."""

from __future__ import annotations

import collections
import csv
import hashlib
import json
import pathlib
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO

from orchestrator import __version__
from orchestrator.schemas import AppendOnlyProvenanceLog


def _detect_seq_fmt(path: pathlib.Path) -> str:
    if path.suffix.lower() in (".fq", ".fastq"):
        return "fastq"
    with path.open("r", encoding="utf-8") as fh:
        first = fh.readline().strip()
    return "fastq" if first.startswith("@") else "fasta"


def _resolve_sheet_path(raw_val: str, base_dir: pathlib.Path) -> Optional[pathlib.Path]:
    val = (raw_val or "").strip()
    if not val or val.upper() in ("NA", "NONE", "NULL", "-"):
        return None
    cand = pathlib.Path(val)
    if cand.is_absolute() or cand.exists():
        return cand
    rel = base_dir / cand
    return rel


def _revcomp(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGTRYMKWSBDHVN", "TGCAYRKMWSVHDBN"))[::-1]


def _build_replicate_pileup_tsv(
    ref_seq: str,
    replicate_reads: List[Tuple[str, pathlib.Path, Optional[pathlib.Path]]],
    out_tsv: pathlib.Path,
    k_anchor: int = 19,
) -> bool:
    """Anchor replicate reads onto curated reference and emit multi-replicate allele count TSV."""
    ref = ref_seq.upper()
    g_len = len(ref)
    if g_len < k_anchor + 5 or len(replicate_reads) < 2:
        return False

    # Index unique k-mer anchors on reference
    anchor_pos: Dict[str, int] = {}
    ambig: set[str] = set()
    for i in range(0, g_len - k_anchor + 1):
        km = ref[i : i + k_anchor]
        if km in anchor_pos:
            ambig.add(km)
        else:
            anchor_pos[km] = i
    for km in ambig:
        anchor_pos.pop(km, None)

    n_reps = len(replicate_reads)
    # depth[rep_idx][pos] and alt_counts[rep_idx][pos][base]
    depths = [[0] * g_len for _ in range(n_reps)]
    alt_counts: List[List[Dict[str, int]]] = [
        [collections.Counter() for _ in range(g_len)] for _ in range(n_reps)
    ]

    for rep_idx, (_rep_id, r1_p, r2_p) in enumerate(replicate_reads):
        for rp in (r1_p, r2_p):
            if rp is None or not rp.exists():
                continue
            fmt = _detect_seq_fmt(rp)
            for rec in SeqIO.parse(str(rp), fmt):
                raw_s = str(rec.seq).upper()
                for cand in (raw_s, _revcomp(raw_s)):
                    start_ref: Optional[int] = None
                    for off in (0, max(0, len(cand) - k_anchor), len(cand) // 2):
                        if off + k_anchor <= len(cand):
                            hit = anchor_pos.get(cand[off : off + k_anchor])
                            if hit is not None:
                                cand_start = hit - off
                                if 0 <= cand_start and cand_start + len(cand) <= g_len:
                                    start_ref = cand_start
                                    break
                    if start_ref is None:
                        continue
                    mismatches = sum(
                        1 for a, b in zip(cand, ref[start_ref : start_ref + len(cand)]) if a != b
                    )
                    if mismatches <= max(2, int(0.08 * len(cand))):
                        for j, ch in enumerate(cand):
                            if ch in "ACGT":
                                p_idx = start_ref + j
                                depths[rep_idx][p_idx] += 1
                                if ch != ref[p_idx]:
                                    alt_counts[rep_idx][p_idx][ch] += 1
                        break

    rows: List[Dict[str, Any]] = []
    bg_added = 0
    for p_idx in range(g_len):
        ref_b = ref[p_idx]
        rep_n = [depths[r][p_idx] for r in range(n_reps)]
        if min(rep_n) < 20:
            continue
        # Identify dominant alternate base across replicates
        pooled_alts: collections.Counter[str] = collections.Counter()
        for r in range(n_reps):
            pooled_alts.update(alt_counts[r][p_idx])
        if pooled_alts:
            alt_b = pooled_alts.most_common(1)[0][0]
            rep_k = [alt_counts[r][p_idx].get(alt_b, 0) for r in range(n_reps)]
            row: Dict[str, Any] = {"pos": p_idx + 1, "ref": ref_b, "alt": alt_b}
            for r in range(n_reps):
                row[f"rep{r + 1}_k"] = rep_k[r]
                row[f"rep{r + 1}_n"] = rep_n[r]
            rows.append(row)
        elif bg_added < 20:
            # Include representative monomorphic background loci for Empirical Bayes error calibration
            alt_b = "G" if ref_b != "G" else "A"
            row = {"pos": p_idx + 1, "ref": ref_b, "alt": alt_b}
            for r in range(n_reps):
                row[f"rep{r + 1}_k"] = 0
                row[f"rep{r + 1}_n"] = rep_n[r]
            rows.append(row)
            bg_added += 1

    if not rows:
        return False

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["pos", "ref", "alt"] + [
        col for r in range(n_reps) for col in (f"rep{r + 1}_k", f"rep{r + 1}_n")
    ]
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return True


def execute_end_to_end_pipeline(
    r1_path: pathlib.Path,
    output_dir: pathlib.Path,
    r2_path: Optional[pathlib.Path] = None,
    seeds_path: Optional[pathlib.Path] = None,
    host_path: Optional[pathlib.Path] = None,
    replicate_manifest: Optional[pathlib.Path] = None,
    cohort_counts: Optional[pathlib.Path] = None,
) -> Dict[str, Any]:
    """Orchestrate single-sample `viro_drift` -> `viral_genome_tools` -> `viralbio_reporting` (+ optional variants & networks)."""
    from viral_genome_tools.cli import main as vgt_main
    from viralbio_reporting.cli import main as vbr_main
    from viro_drift.cli import main as drift_main

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


def execute_cohort_pipeline(
    samplesheet_path: pathlib.Path,
    output_dir: pathlib.Path,
    seeds_path: Optional[pathlib.Path] = None,
    host_path: Optional[pathlib.Path] = None,
    replicate_type: str = "technical",
    alpha_fdr: float = 0.01,
    bootstraps: int = 100,
    stability: float = 0.80,
) -> Dict[str, Any]:
    """Execute multi-sample, multi-replicate batch cohort orchestration across all 5 packages."""
    from quasispecies_profiler.cli import main as qp_main
    from viro_networks.cli import main as vn_main

    delim = "," if samplesheet_path.suffix.lower() == ".csv" else "\t"
    with samplesheet_path.open("r", encoding="utf-8") as fh:
        raw_rows = list(csv.DictReader(fh, delimiter=delim))

    if not raw_rows:
        raise ValueError(f"Cohort samplesheet contains zero rows: {samplesheet_path}")
    cols = set(raw_rows[0].keys())
    if not {"sample_id", "r1_path"}.issubset(cols):
        raise ValueError("Cohort samplesheet must include 'sample_id' and 'r1_path' columns.")

    base_dir = samplesheet_path.resolve().parent
    grouped: Dict[str, List[Dict[str, Any]]] = collections.OrderedDict()

    for idx, row in enumerate(raw_rows, start=1):
        sid = (row.get("sample_id") or "").strip()
        if not sid:
            raise ValueError(f"Row {idx}: missing sample_id")
        r1_p = _resolve_sheet_path(row.get("r1_path", ""), base_dir)
        if r1_p is None or not r1_p.exists():
            raise FileNotFoundError(f"Row {idx} ({sid}): r1_path not found: {row.get('r1_path')}")
        r2_p = _resolve_sheet_path(row.get("r2_path", ""), base_dir)
        if r2_p is not None and not r2_p.exists():
            raise FileNotFoundError(f"Row {idx} ({sid}): r2_path not found: {row.get('r2_path')}")
        rep_man = _resolve_sheet_path(row.get("replicate_manifest", ""), base_dir)
        rep_id = (row.get("replicate_id") or "").strip() or f"Rep{len(grouped.get(sid, [])) + 1}"
        group_label = (row.get("group") or row.get("condition") or "Unspecified").strip()
        grouped.setdefault(sid, []).append(
            {
                "sample_id": sid,
                "replicate_id": rep_id,
                "r1_path": r1_p,
                "r2_path": r2_p,
                "group": group_label,
                "replicate_manifest": rep_man,
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    prov = AppendOnlyProvenanceLog(output_dir / "cohort_provenance.json")
    sheet_hash = hashlib.sha256(samplesheet_path.read_bytes()).hexdigest()[:16]

    samples_root = output_dir / "samples"
    per_sample_summaries: Dict[str, Any] = {}
    catalog_seqs: Dict[str, str] = collections.OrderedDict()
    sample_read_files: Dict[str, List[pathlib.Path]] = {}
    metadata_rows: List[Dict[str, str]] = []

    for sid, rep_rows in grouped.items():
        s_dir = samples_root / sid
        s_dir.mkdir(parents=True, exist_ok=True)
        metadata_rows.append(
            {
                "sample_id": sid,
                "replicate_count": str(len(rep_rows)),
                "group": rep_rows[0]["group"],
            }
        )
        sample_read_files[sid] = [
            p for r in rep_rows for p in (r["r1_path"], r["r2_path"]) if p is not None
        ]

        # Pool replicate R1 / R2 files when multiple replicates exist for consensus assembly/curation
        if len(rep_rows) == 1:
            pooled_r1 = rep_rows[0]["r1_path"]
            pooled_r2 = rep_rows[0]["r2_path"]
        else:
            pool_dir = s_dir / "00_pooled_replicates"
            pool_dir.mkdir(parents=True, exist_ok=True)
            ext = rep_rows[0]["r1_path"].suffix or ".fasta"
            pooled_r1 = pool_dir / f"{sid}_pooled_R1{ext}"
            with pooled_r1.open("w", encoding="utf-8") as out1:
                for rr in rep_rows:
                    txt = rr["r1_path"].read_text(encoding="utf-8")
                    out1.write(txt if txt.endswith("\n") else txt + "\n")
            if all(rr["r2_path"] is not None for rr in rep_rows):
                pooled_r2 = pool_dir / f"{sid}_pooled_R2{ext}"
                with pooled_r2.open("w", encoding="utf-8") as out2:
                    for rr in rep_rows:
                        txt2 = rr["r2_path"].read_text(encoding="utf-8")
                        out2.write(txt2 if txt2.endswith("\n") else txt2 + "\n")
            else:
                pooled_r2 = None

        explicit_manifest = next(
            (rr["replicate_manifest"] for rr in rep_rows if rr["replicate_manifest"] is not None),
            None,
        )
        s_res = execute_end_to_end_pipeline(
            r1_path=pooled_r1,
            output_dir=s_dir,
            r2_path=pooled_r2,
            seeds_path=seeds_path,
            host_path=host_path,
            replicate_manifest=explicit_manifest,
        )

        # Register curated genomes into cohort-wide non-redundant reference catalog (deduplicate exact sequences)
        cur_fa = s_dir / "02_genome_curation" / "curated_genomes.fasta"
        sample_curated: List[Tuple[str, str]] = []
        if cur_fa.exists():
            for rec in SeqIO.parse(str(cur_fa), "fasta"):
                seq_u = str(rec.seq).upper()
                sample_curated.append((rec.id, seq_u))
                if seq_u not in catalog_seqs.values():
                    taxon_label = rec.id if rec.id not in catalog_seqs else f"{sid}_{rec.id}"
                    catalog_seqs[taxon_label] = seq_u

        # Automatic multi-replicate variant calling when M >= 2 replicates and no explicit manifest was supplied
        if explicit_manifest is None and len(rep_rows) >= 2 and sample_curated:
            auto_tsv = s_dir / "04_quasispecies_variants" / "auto_replicate_counts.tsv"
            rep_tuples = [
                (rr["replicate_id"], rr["r1_path"], rr["r2_path"]) for rr in rep_rows
            ]
            built = _build_replicate_pileup_tsv(sample_curated[0][1], rep_tuples, auto_tsv)
            if built:
                rc_qp = qp_main(
                    [
                        "call",
                        "--manifest", str(auto_tsv),
                        "--replicate-type", replicate_type,
                        "--alpha", str(alpha_fdr),
                        "-o", str(s_dir / "04_quasispecies_variants"),
                    ]
                )
                s_res["stages"]["quasispecies_profiler_exit"] = rc_qp

        per_sample_summaries[sid] = s_res["stages"]

    # Build cohort-wide non-redundant reference FASTA & cross-sample abundance count matrix
    matrices_dir = output_dir / "04_cohort_matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)
    meta_tsv = matrices_dir / "cohort_metadata.tsv"
    with meta_tsv.open("w", encoding="utf-8", newline="") as mfh:
        mw = csv.DictWriter(mfh, fieldnames=["sample_id", "replicate_count", "group"], delimiter="\t")
        mw.writeheader()
        mw.writerows(metadata_rows)

    counts_tsv = matrices_dir / "cohort_abundance_counts.tsv"
    catalog_fa = matrices_dir / "cohort_viral_catalog.fasta"
    taxa_list = list(catalog_seqs.keys())

    with catalog_fa.open("w", encoding="utf-8") as cfh:
        for t_name, t_seq in catalog_seqs.items():
            cfh.write(f">{t_name}\n{t_seq}\n")

    # Index 21-mers per catalog taxon for fast cross-sample read/contig abundance quantification
    taxon_kmers: Dict[str, set[str]] = {}
    for t_name, t_seq in catalog_seqs.items():
        ks = {t_seq[i : i + 21] for i in range(0, max(0, len(t_seq) - 20))}
        ks |= {_revcomp(x) for x in ks}
        taxon_kmers[t_name] = ks

    with counts_tsv.open("w", encoding="utf-8", newline="") as tfh:
        writer = csv.writer(tfh, delimiter="\t")
        writer.writerow(["sample_id"] + taxa_list)
        for sid in grouped:
            row_counts = [0] * len(taxa_list)
            for fpath in sample_read_files[sid]:
                fmt = _detect_seq_fmt(fpath)
                for rec in SeqIO.parse(str(fpath), fmt):
                    s = str(rec.seq).upper()
                    if len(s) < 21:
                        continue
                    probe_kmers = {s[0:21], s[len(s) // 2 : len(s) // 2 + 21], s[-21:]}
                    for t_idx, t_name in enumerate(taxa_list):
                        if probe_kmers & taxon_kmers[t_name]:
                            row_counts[t_idx] += 1
            writer.writerow([sid] + row_counts)

    counts_hash = hashlib.sha256(counts_tsv.read_bytes()).hexdigest()[:16]
    prov.append_event("cohort_abundance_aggregation", "virome-suite", __version__, sheet_hash, counts_hash)

    # Stage 5: Ecological Co-Occurrence Network & Sloan Neutral Model across Cohort (`viro-networks`)
    net_dir = output_dir / "05_cohort_ecological_networks"
    if len(grouped) >= 5 and len(taxa_list) >= 2:
        rc_vn = vn_main(
            [
                "infer",
                "--counts", str(counts_tsv),
                "--metadata", str(meta_tsv),
                "--bootstraps", str(bootstraps),
                "--stability", str(stability),
                "-o", str(net_dir),
            ]
        )
        net_status: Any = "SUCCESS" if rc_vn == 0 else f"EXIT_{rc_vn}"
        prov.append_event("cohort_network_inference", "viro-networks", __version__, counts_hash, str(net_status))
    else:
        net_dir.mkdir(parents=True, exist_ok=True)
        net_status = (
            "ABSTAIN_INSUFFICIENT_SAMPLES"
            if len(grouped) < 5
            else "SKIPPED_FEWER_THAN_2_TAXA"
        )
        (net_dir / "network_summary.json").write_text(
            json.dumps({"status": net_status, "n_samples": len(grouped), "n_taxa": len(taxa_list)}, indent=2),
            encoding="utf-8",
        )

    cohort_summary = {
        "orchestrator_version": __version__,
        "n_samples": len(grouped),
        "total_replicate_libraries": sum(len(v) for v in grouped.values()),
        "catalog_viral_taxa": len(taxa_list),
        "cohort_abundance_matrix": str(counts_tsv),
        "cohort_metadata": str(meta_tsv),
        "viro_networks_status": net_status,
        "samples": per_sample_summaries,
        "provenance_log": str(output_dir / "cohort_provenance.json"),
    }
    (output_dir / "cohort_summary.json").write_text(json.dumps(cohort_summary, indent=2), encoding="utf-8")
    return cohort_summary
