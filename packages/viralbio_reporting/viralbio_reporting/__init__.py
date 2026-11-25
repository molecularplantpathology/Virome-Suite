"""Package 4: viralbio_reporting — Seamless Circular Seam-Wrapping ORF Annotation & NCBI Deposition."""

from viralbio_reporting.circular_annotator import (
    AnnotatedCDS,
    annotate_genome_orfs,
    reconstruct_cds_nt_from_intervals,
    revcomp,
    validate_cds_feature,
)
from viralbio_reporting.ncbi_packager import (
    ProvenanceRecord,
    SequenceProvenanceChain,
    build_ncbi_submission_package,
    format_ncbi_feature_table,
    validate_tbl_syntax,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "AnnotatedCDS",
    "ProvenanceRecord",
    "SequenceProvenanceChain",
    "annotate_genome_orfs",
    "build_ncbi_submission_package",
    "format_ncbi_feature_table",
    "reconstruct_cds_nt_from_intervals",
    "revcomp",
    "validate_cds_feature",
    "validate_tbl_syntax",
]
