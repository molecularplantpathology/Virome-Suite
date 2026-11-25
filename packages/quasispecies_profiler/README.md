# quasispecies-profiler

**Replicate-Aware Hierarchical Empirical Bayes Beta-Binomial Intra-Host Variant Caller & Diversity Profiler**

Part of the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

## Key Features
- **PCR Amplification Jackpot Filtering:** Separates genuine intra-host viral variants (0.5%–5% VAF) from single-replicate RT-PCR overdispersion artifacts across technical and biological replicates.
- **Ultra-Deep Numerical Stability:** Log-space `betaln`/`gammaln` Beta-Binomial optimization stable beyond $100,000\times$ sequencing coverage.
- **Empirical Null Calibration & Profile CIs:** Vectorized null parametric bootstrap calibration with Benjamini-Hochberg FDR control, 95% Profile Likelihood Confidence Intervals, and depth-dependent Limits of Detection.
- **Population Genetics Metrics:** Genome-wide Shannon entropy ($H$), nucleotide diversity ($\pi$), and Nei-Gojobori $d_N/d_S$ selection ratios.

## Installation

```bash
pip install quasispecies-profiler
```

## Quick Start

```bash
quasispecies-profiler call --manifest replicates_counts.tsv -o out_variants/
```
