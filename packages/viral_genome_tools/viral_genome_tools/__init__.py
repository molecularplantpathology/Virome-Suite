"""Package 2: viral_genome_tools — Circular Viral Monomer Curation, Landmark Rotation & E_junction Scoring."""

from viral_genome_tools.junction_scorer import (
    JunctionReadEvidence,
    TopologyDecision,
    classify_topology,
    compute_coverage_uniformity,
    compute_e_junction_score,
    compute_host_flank_ratio,
    expected_spanning_pairs,
    expected_split_reads,
    profile_reads_across_seam,
)
from viral_genome_tools.landmark_rotator import (
    MULTI_FAMILY_LANDMARK_REGISTRY,
    LandmarkMotif,
    LandmarkRotationResult,
    rotate_to_landmark,
)
from viral_genome_tools.monomer_curator import (
    MonomerCurationResult,
    detect_terminal_inverted_repeat,
    revcomp,
    trim_concatemer_to_monomer,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "JunctionReadEvidence",
    "LandmarkMotif",
    "LandmarkRotationResult",
    "MULTI_FAMILY_LANDMARK_REGISTRY",
    "MonomerCurationResult",
    "TopologyDecision",
    "classify_topology",
    "compute_coverage_uniformity",
    "compute_e_junction_score",
    "compute_host_flank_ratio",
    "detect_terminal_inverted_repeat",
    "expected_spanning_pairs",
    "expected_split_reads",
    "profile_reads_across_seam",
    "revcomp",
    "rotate_to_landmark",
    "trim_concatemer_to_monomer",
]
