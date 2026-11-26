# viro-networks

[![PyPI](https://img.shields.io/pypi/v/viro-networks.svg?color=blue)](https://pypi.org/project/viro-networks/)
[![CI & Coverage](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml/badge.svg)](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Compositional CLR Consensus Co-Occurrence Networks, Sloan Neutral Community Modeling, and Pure-Python Topological Hub Perturbation.**

`viro-networks` is a standalone package within the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

---

## Ecological & Statistical Methodology

1. **Sample-Size Safety Gate & CLR Normalization (`viro_networks.clr_consensus`):**
   - Enforces a minimum cohort size guard ($N_{\text{samples}} < 5 \implies$ `CohortInsufficientError` / CLI exit code `4` `ABSTAIN_INSUFFICIENT_SAMPLES`) to prevent underpowered false-positive networks.
   - Applies Centered Log-Ratio (CLR) transformation across a pseudocount sensitivity grid $\delta \in \{0.1, 0.5, 1.0\}$ and reports edge Jaccard concordance across pseudocount extremes.
2. **Three-Estimator Bootstrap Consensus Network (`viro_networks.clr_consensus`):**
   Retains an undirected signed edge $(i, j)$ if and only if **$\ge 2$ of 3 compositional estimators** agree in sign above threshold AND sample-level bootstrap stability across $B = 100$ resamples meets `stability_threshold` ($\ge 0.80$):
   - **Positive-Floored SparCC:** Log-ratio variation basis variance floored at $\omega_i^2 \ge \max(0.10\,\mathrm{Var}(\mathrm{CLR}_{\cdot i}), 10^{-6})$ ($|r^{\text{SparCC}}| \ge 0.35$).
   - **Shrinkage Graphical Partial Correlation:** Regularized CLR correlation precision inversion $\mathbf{\Theta} = ((1 - \alpha)\mathbf{R}_{\text{CLR}} + \alpha \mathbf{I})^{-1}$ with $\alpha = 0.15$ ($|\omega_{ij}| \ge 0.25$).
   - **CLR Rank Association:** Spearman rank correlation on CLR coordinates ($|\rho^{\text{rank}}_{ij}| \ge 0.35$).
3. **Domain-Clipped Sloan Neutral Community Model (`viro_networks.sloan_neutral`):**
   Fits Sloan's (2006) neutral immigration-drift parameter $N_e m \in [10^{-3}, 10^4]$ predicting taxon occupancy $\hat{f}_i = 1 - \mathrm{I}_d(2 N_e m p_i, 2 N_e m (1 - p_i))$ from mean relative abundance $p_i$, reporting zero-variance-guarded $R^2$, 95% Wilson score intervals, and Benjamini-Hochberg FDR-adjusted binomial departure $q$-values (`neutral_drift_consistent`, `candidate_above_neutral_partition`, `candidate_below_neutral_partition`).
4. **Topological Hub Perturbation (`viro_networks.hub_perturbation`):**
   Computes degree centrality, exact Brandes shortest-path betweenness centrality, largest connected component (LCC) retention, and eigenvalue natural connectivity drop $\Delta \bar{\lambda}\%$ upon *in silico* single-node removal using pure NumPy/Python (zero `networkx` dependency).

---

## Installation

```bash
pip install viro-networks
```

---

## Command-Line Usage (`viro-networks`)

### Input Cohort Count Matrix (`virome_cohort_counts.tsv`)

Rows represent samples ($N \ge 5$) and columns represent viral taxa/OTUs (raw mapped read counts or RPKM/CPM):

```tsv
sample_id	Badnavirus_A	Badnavirus_B	Poacevirus_SCSMV	Polerovirus_SCYLV	Begomovirus_C	Mastrevirus_D
SAMPLE_01	169	560	0	0	552	68
SAMPLE_02	1327	1053	3486	2921	665	112
SAMPLE_03	833	0	0	0	269	320
...
```

### Run Network & Neutral Model Inference

```bash
viro-networks infer \
  --counts virome_cohort_counts.tsv \
  --metadata cohort_metadata.tsv \
  --bootstraps 100 \
  --stability 0.80 \
  --output network_results/
```

### CLI Options (`viro-networks infer`)

| Argument | Required | Default | Description |
| :--- | :---: | :---: | :--- |
| `--counts` | Yes | — | Tab-separated sample $\times$ taxon abundance matrix TSV (first column = `sample_id`) |
| `--metadata` | No | `None` | Optional sample metadata TSV |
| `--bootstraps` | No | `100` | Number of sample-level bootstrap resamples for consensus stability |
| `--stability` | No | `0.80` | Minimum bootstrap selection frequency in $[0, 1]$ required to retain an edge |
| `-o`, `--output` | Yes | — | Output directory for network edges, Sloan neutral fit, and hub rankings |

**Exit Codes:** `0` = success, `2` = missing/malformed count table, `4` = safe abstention when cohort has $N < 5$ samples (`ABSTAIN_INSUFFICIENT_SAMPLES`).

---

## Python API Usage

```python
import pandas as pd
from viro_networks.clr_consensus import infer_consensus_network
from viro_networks.sloan_neutral import fit_sloan_neutral_model
from viro_networks.hub_perturbation import evaluate_hub_perturbation

df = pd.read_csv("virome_cohort_counts.tsv", sep="\t", index_col=0)
counts = df.to_numpy(dtype=float)
taxa = list(df.columns)

edges, net_diag = infer_consensus_network(counts, taxa, n_boot=100, stability_threshold=0.80)
sloan_rows, sloan_diag = fit_sloan_neutral_model(counts, taxa)
hubs = evaluate_hub_perturbation(taxa, [(e.taxon_a, e.taxon_b) for e in edges])
```

---

## Output Files & Column Schema

- **`consensus_edges.tsv`:** `taxon_a`, `taxon_b`, `sign` (`+1` or `-1`), `sparcc_r`, `partial_omega`, `spearman_rho`, `estimators_supporting`, `bootstrap_stability`, `pseudocount_concordance`.
- **`sloan_neutral_fit.tsv`:** `taxon`, `mean_rel_abundance`, `observed_occurrence`, `predicted_occurrence`, `wilson_95_low`, `wilson_95_high`, `p_value`, `q_value`, `partition_class`.
- **`topological_hub_candidates.tsv`:** `taxon`, `degree`, `degree_centrality`, `betweenness_centrality`, `natural_connectivity_drop_pct`, `lcc_fraction_after_removal`.
- **`network_summary.json`:** Cohort sample/taxon counts, pseudocount Jaccard concordance, $N_e m$ estimate, and $R^2$ goodness-of-fit.
