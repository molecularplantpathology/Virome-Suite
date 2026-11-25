"""Package 3: quasispecies_profiler — Replicate Beta-Binomial Intra-Host Variant Caller."""

from quasispecies_profiler.beta_binomial_caller import (
    VariantCallRecord,
    call_variants_across_replicates,
    compute_profile_ci_95,
    estimate_background_error_and_prior,
    log_beta_binomial_pmf,
    penalized_log_likelihood,
)
from quasispecies_profiler.diversity_metrics import (
    compute_nei_gojobori_dnds,
    compute_nucleotide_diversity_pi,
    compute_shannon_entropy,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "VariantCallRecord",
    "call_variants_across_replicates",
    "compute_nei_gojobori_dnds",
    "compute_nucleotide_diversity_pi",
    "compute_profile_ci_95",
    "compute_shannon_entropy",
    "estimate_background_error_and_prior",
    "log_beta_binomial_pmf",
    "penalized_log_likelihood",
]
