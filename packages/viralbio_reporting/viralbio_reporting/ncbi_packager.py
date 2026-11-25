"""NCBI 5-Column Feature Table Exporter, Provenance Chain & `table2asn` Harness (`viralbio_reporting.ncbi_packager`)."""

from __future__ import annotations

import dataclasses
import datetime
import hashlib
import json
import pathlib
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Sequence, Tuple

from viralbio_reporting.circular_annotator import AnnotatedCDS, validate_cds_feature


@dataclasses.dataclass
class ProvenanceRecord:
    """Hash-chained transformation provenance entry."""

    step: str
    input_sha256: str
    output_sha256: str
    details: Dict[str, Any]
    previous_hash: str
    timestamp_utc: str = dataclasses.field(
        default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat()
    )

    def record_hash(self) -> str:
        payload = json.dumps(dataclasses.asdict(self), sort_keys=True).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()


class SequenceProvenanceChain:
    """Self-contained append-only SHA-256 provenance logger."""

    def __init__(self, log_path: pathlib.Path | str) -> None:
        self.log_path = pathlib.Path(log_path)
        self.records: List[ProvenanceRecord] = []

    def append(
        self,
        step: str,
        input_sha256: str,
        output_sha256: str,
        details: Optional[Dict[str, Any]] = None,
    ) -> ProvenanceRecord:
        prev = self.records[-1].record_hash() if self.records else "GENESIS"
        rec = ProvenanceRecord(
            step=step,
            input_sha256=input_sha256,
            output_sha256=output_sha256,
            details=details or {},
            previous_hash=prev,
        )
        self.records.append(rec)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        self.log_path.write_text(
            json.dumps(
                {
                    "record_count": len(self.records),
                    "head_sha256": rec.record_hash(),
                    "records": [
                        {**dataclasses.asdict(r), "entry_sha256": r.record_hash()}
                        for r in self.records
                    ],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return rec

    def verify_integrity(self) -> bool:
        expected_prev = "GENESIS"
        for r in self.records:
            if r.previous_hash != expected_prev:
                return False
            expected_prev = r.record_hash()
        return True


def format_ncbi_feature_table(seq_id: str, features: Sequence[AnnotatedCDS]) -> str:
    """Format compliant NCBI 5-column `.tbl` feature table (`P0-VBR-2`, `P1-VBR-3`)."""
    lines: List[str] = [f">Feature {seq_id}"]
    for feat in features:
        n_iv = len(feat.intervals_1based)
        for idx, (st_1, en_1) in enumerate(feat.intervals_1based):
            s_str = f"<{st_1}" if (idx == 0 and feat.partial_5prime) else str(st_1)
            e_str = f">{en_1}" if (idx == n_iv - 1 and feat.partial_3prime) else str(en_1)
            if idx == 0:
                lines.append(f"{s_str}\t{e_str}\tCDS")
            else:
                lines.append(f"{s_str}\t{e_str}")
        lines.append(f"\t\t\tproduct\t{feat.product}")
        lines.append(f"\t\t\tcodon_start\t{feat.codon_start}")
        lines.append("\t\t\ttransl_table\t1")
        lines.append(f"\t\t\tnote\t{feat.feature_id}")
    return "\n".join(lines) + "\n"


def validate_tbl_syntax(tbl_text: str) -> Tuple[bool, List[str]]:
    """Strict internal syntax check for NCBI 5-column `.tbl` compliance."""
    errors: List[str] = []
    lines = tbl_text.splitlines()
    if not lines or not lines[0].startswith(">Feature "):
        return False, ["Missing '>Feature <SeqID>' header on line 1"]

    seen_cds = False
    for line_no, ln in enumerate(lines[1:], start=2):
        if not ln.strip():
            continue
        if ln.startswith(">Feature "):
            seen_cds = False
            continue
        if "join(" in ln or ".." in ln:
            errors.append(f"Line {line_no}: Illegal flatfile syntax ('join(' or '..') in .tbl")
        cols = ln.split("\t")
        if ln[0] != "\t":
            # Interval line
            if len(cols) not in (2, 3):
                errors.append(f"Line {line_no}: Expected 2 or 3 tab-separated columns on interval line")
            for coord_tok in cols[:2]:
                clean = coord_tok.lstrip("<>")
                if not clean.isdigit() or int(clean) < 1:
                    errors.append(f"Line {line_no}: Invalid 1-based coordinate '{coord_tok}'")
            if len(cols) == 3:
                seen_cds = True
            elif not seen_cds:
                errors.append(f"Line {line_no}: Continuation interval before feature key")
        else:
            # Qualifier line: must start with 3 leading tabs (`\t\t\tkey\tval`)
            if len(cols) < 4 or cols[0] != "" or cols[1] != "" or cols[2] != "":
                errors.append(f"Line {line_no}: Qualifier must begin with 3 leading tabs")

    return len(errors) == 0, errors


def build_ncbi_submission_package(
    seq_id: str,
    sequence: str,
    features: Sequence[AnnotatedCDS],
    output_dir: pathlib.Path,
    topology: str = "circular",
    organism: str = "Unclassified virus",
    host: str = "Plant host",
    country: str = "Unspecified",
    run_table2asn: bool = True,
) -> Dict[str, Any]:
    """Export `.fsa`, `.tbl`, `source.tsv`, provenance JSON, and run validation."""
    output_dir.mkdir(parents=True, exist_ok=True)
    clean_seq = sequence.upper().strip()

    # Validate every CDS feature biologically before writing
    bio_errors: List[str] = []
    for feat in features:
        ok, msg = validate_cds_feature(clean_seq, feat)
        if not ok:
            bio_errors.append(f"{feat.feature_id}: {msg}")

    fsa_path = output_dir / f"{seq_id}.fsa"
    tbl_path = output_dir / f"{seq_id}.tbl"
    src_path = output_dir / "source.tsv"
    prov_path = output_dir / "provenance_chain.json"

    topo_tag = "[topology=circular]" if topology.lower() == "circular" else "[topology=linear]"
    fsa_text = f">{seq_id} [organism={organism}] {topo_tag} [gcode=1]\n{clean_seq}\n"
    tbl_text = format_ncbi_feature_table(seq_id, features)
    src_text = f"Sequence_ID\tOrganism\tHost\tCountry\n{seq_id}\t{organism}\t{host}\t{country}\n"

    fsa_path.write_text(fsa_text, encoding="utf-8")
    tbl_path.write_text(tbl_text, encoding="utf-8")
    src_path.write_text(src_text, encoding="utf-8")

    syntax_ok, syntax_errors = validate_tbl_syntax(tbl_text)

    seq_hash = hashlib.sha256(clean_seq.encode("utf-8")).hexdigest()
    tbl_hash = hashlib.sha256(tbl_text.encode("utf-8")).hexdigest()
    chain = SequenceProvenanceChain(prov_path)
    chain.append("normalize_fasta", seq_hash, hashlib.sha256(fsa_text.encode("utf-8")).hexdigest(), {"topology": topology})
    chain.append("generate_ncbi_tbl", seq_hash, tbl_hash, {"feature_count": len(features)})

    t2a_status = "SKIPPED_BINARY_NOT_IN_PATH"
    t2a_fatal = 0
    t2a_exe = shutil.which("table2asn")
    if run_table2asn and t2a_exe is not None:
        proc = subprocess.run(
            [t2a_exe, "-i", str(fsa_path), "-f", str(tbl_path), "-V", "vb"],
            capture_output=True,
            text=True,
            check=False,
        )
        t2a_status = "EXECUTED" if proc.returncode == 0 else f"EXIT_{proc.returncode}"
        val_file = output_dir / f"{seq_id}.val"
        if val_file.exists():
            for ln in val_file.read_text(encoding="utf-8", errors="ignore").splitlines():
                if "ERROR:" in ln or "REJECT:" in ln:
                    t2a_fatal += 1

    return {
        "seq_id": seq_id,
        "length_bp": len(clean_seq),
        "topology": topology,
        "feature_count": len(features),
        "origin_spanning_count": sum(1 for f in features if f.spans_origin),
        "biological_validation_passed": len(bio_errors) == 0,
        "biological_errors": bio_errors,
        "tbl_syntax_passed": syntax_ok,
        "tbl_syntax_errors": syntax_errors,
        "table2asn_status": t2a_status,
        "table2asn_fatal_errors": t2a_fatal,
        "provenance_verified": chain.verify_integrity(),
        "sequence_sha256": seq_hash,
        "tbl_sha256": tbl_hash,
    }
