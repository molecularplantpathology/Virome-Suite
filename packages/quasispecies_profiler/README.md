# quasispecies-profiler

[![PyPI](https://img.shields.io/pypi/v/quasispecies-profiler.svg)](https://pypi.org/project/quasispecies-profiler/)
[![CI & Coverage](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml/badge.svg)](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Replicate-Aware Hierarchical Empirical Bayes Beta-Binomial Intra-Host Variant Calling, Jackpot Filtering, and Population Diversity Profiling.**

`quasispecies-profiler` is a standalone package within the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

---

## Statistical Methodology

1. **Physical Limit-of-Detection (LOD) Pre-Gates (`quasispecies_profiler.beta_binomial_caller`):**
   - Mean replicate depth $\bar{n}_i < 100\times \implies$ `UNASSESSED_LOW_DEPTH`.
   - $100\times \le \bar{n}_i < 1,000\times$ with pooled VAF $< 3.0\% \implies$ `BELOW_LOD_AT_DEPTH`.
   - $\bar{n}_i \ge 1,000\times \implies$ evaluates candidate variants down to $0.5\%$ VAF ($0.005$).
2. **Maximum-Replicate Background Error Estimation:**
   Selects monomorphic reference sites strictly via $\max_{1 \le m \le M}(k_{jm}/n_{jm}) < 0.002$, preventing single-replicate amplification jackpots from contaminating the substitution-context background error rate $\epsilon_{\text{bg}}(c_i)$.
3. **Log-Space Penalized Beta-Binomial Likelihood & Technical Jackpot Veto:**
   Evaluates exact log-PMFs via `scipy.special.betaln` (numerically stable beyond $100,000\times$ coverage) with diffuse log-dispersion prior scale $\sigma_{\rho} \ge 3.5$. In `--replicate-type technical` mode ($M \ge 3$), loci active above background in $\le 1$ replicate (or exhibiting high overdispersion $\hat{\rho}_i > 0.005$ with inter-replicate $\mathrm{CV} > 1.0$) are vetoed as `REJECTED_PCR_JACKPOT`.
4. **Vectorized Null Parametric Bootstrap & Profile Likelihood CIs:**
   - Draws $B = 2,000$ vectorized Beta-Binomial null replicates per candidate locus (with asymptotic $\frac{1}{2}\chi^2_1$ extreme-tail resolution) and controls false discovery rate via Benjamini-Hochberg ($q \le \alpha$) across the evaluable candidate family.
   - Computes sign-bracketed 95% profile likelihood confidence intervals on VAF via `scipy.optimize.brentq`.
5. **Population Genetics Metrics (`quasispecies_profiler.diversity_metrics`):**
   Computes genome-normalized Shannon entropy $H$ (bits), unbiased nucleotide diversity $\pi$ ($n_i < 2$ guarded), and Nei-Gojobori synonymous/nonsynonymous selection ratios ($d_N, d_S, d_N/d_S$ with Jukes-Cantor saturation and $d_S = 0$ guards).

---

## Installation

```bash
pip install quasispecies-profiler
```

---

## Command-Line Usage (`quasispecies-profiler`)

### Input Replicate Manifest TSV Format (`replicate_counts.tsv`)

Supports any number of replicates $M \ge 2$ (`rep1_k/rep1_n`, `rep2_k/rep2_n`, `rep3_k/rep3_n`, ...):

```tsv
pos	ref	alt	rep1_k	rep1_n	rep2_k	rep2_n	rep3_k	rep3_n
150	A	G	46	1000	44	1000	45	1000
320	C	T	210	10000	195	10000	202	10000
1120	A	T	180	10000	2	10000	1	10000
1700	G	T	3	10000	2	10000	4	10000
```

### Run Variant Calling

```bash
quasispecies-profiler call \
  --manifest replicate_counts.tsv \
  --ref reference_monomer.fasta \
  --replicate-type technical \
  --alpha 0.01 \
  --output variants_out/
```

### CLI Options (`quasispecies-profiler call`)

| Argument | Required | Default | Description |
| :--- | :---: | :---: | :--- |
| `--manifest` | Yes | — | Tab-separated replicate count table (`pos`, `ref`, `alt`, `rep*_k`, `rep*_n`) |
| `--ref` | No | `None` | Optional reference FASTA used to normalize genome-wide diversity ($\pi, H$) by length $L$ |
| `--replicate-type` | No | `technical` | Replication model: `technical` (enables single-replicate jackpot veto) or `biological` |
| `--alpha` | No | `0.01` | Target Benjamini-Hochberg false discovery rate threshold |
| `-o`, `--output` | Yes | — | Output directory for `variants_called.tsv` and `quasispecies_summary.json` |

---

## Python API Usage

```python
from quasispecies_profiler.beta_binomial_caller import call_variants_across_replicates
from quasispecies_profiler.diversity_metrics import (
    compute_nucleotide_diversity_pi,
    compute_shannon_entropy,
    compute_nei_gojobori_dnds,
)

sites = [
    {"pos": 150, "ref": "A", "alt": "G", "k": [46, 44, 45], "n": [1000, 1000, 1000]},
    {"pos": 1120, "ref": "A", "alt": "T", "k": [180, 2, 1], "n": [10000, 10000, 10000]},
    {"pos": 1700, "ref": "G", "alt": "T", "k": [3, 2, 4], "n": [10000, 10000, 10000]},
]
calls = call_variants_across_replicates(sites, replicate_type="technical", alpha_fdr=0.01)
for c in calls:
    print(c.pos, c.status, c.is_called, f"VAF={c.vaf_hat:.4f} [{c.ci_95_low:.4f}, {c.ci_95_high:.4f}]", f"q={c.q_value:.3e}")
```

---

## Output Files & Column Schema

- **`variants_called.tsv`:** Per-site statistical table containing `pos`, `ref`, `alt`, `mean_depth`, `vaf_hat`, `rho_hat`, `ci_95_low`, `ci_95_high`, `lrt_stat`, `p_value`, `q_value`, `status` (`PASS_ISNV`, `REJECTED_PCR_JACKPOT`, `BACKGROUND_NOISE`, `BELOW_LOD_AT_DEPTH`, `UNASSESSED_LOW_DEPTH`, `MONOMORPHIC_ZERO`, `FILTERED_FDR`), and `is_called`.
- **`quasispecies_summary.json`:** Cohort/sample summary reporting `total_evaluated_loci`, `called_isnv_count`, `rejected_pcr_jackpots`, `shannon_entropy_H`, and `nucleotide_diversity_pi`.
