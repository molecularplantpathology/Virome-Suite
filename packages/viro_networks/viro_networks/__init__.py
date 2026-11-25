"""Package 5: viro_networks — Compositional Virome Co-Occurrence Networks & Sloan Neutral Model."""

from viro_networks.clr_consensus import (
    CohortInsufficientError,
    ConsensusEdge,
    clr_transform,
    compute_shrinkage_partial_corr,
    compute_sparcc_matrix,
    compute_spearman_clr,
    infer_consensus_network,
)
from viro_networks.hub_perturbation import HubVulnerabilityRecord, evaluate_hub_perturbation
from viro_networks.sloan_neutral import (
    SloanTaxonFit,
    fit_sloan_neutral_model,
    sloan_predicted_occurrence,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "CohortInsufficientError",
    "ConsensusEdge",
    "HubVulnerabilityRecord",
    "SloanTaxonFit",
    "clr_transform",
    "compute_shrinkage_partial_corr",
    "compute_sparcc_matrix",
    "compute_spearman_clr",
    "evaluate_hub_perturbation",
    "fit_sloan_neutral_model",
    "infer_consensus_network",
    "sloan_predicted_occurrence",
]
