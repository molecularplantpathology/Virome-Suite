# viro-networks

**Compositional Virome Co-Occurrence Consensus Networks, Sloan Neutral Community Assembly & Hub Vulnerability**

Part of the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

## Key Features
- **Compositional CLR & Pseudocount Sensitivity Grid:** Centered Log-Ratio normalization evaluated across $\delta \in \{0.1, 0.5, 1.0\}$ to guard against sparse zero-inflation bias.
- **Three-Estimator Bootstrap Consensus Network:** Combines positive-floored SparCC, correlation-scaled shrinkage partial correlations, and rank associations across bootstrap resamples ($\ge 80\%$ stability).
- **Sloan Neutral Community Model Solver:** Fits neutral migration parameter $N_e m$ on untransformed relative abundances and occurrence frequencies with Wilson 95% intervals and Benjamini-Hochberg departure flags.
- **Pure-Python Topological Hub Perturbation:** Brandes betweenness centrality and eigenvalue natural connectivity drop upon node removal.

## Installation

```bash
pip install viro-networks
```

## Quick Start

```bash
viro-networks infer --counts virome_cohort_counts.tsv -o out_networks/
```
