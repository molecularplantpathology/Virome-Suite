# Virome-Suite

[![CI/CD](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml/badge.svg)](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

**An Integrated, Modular Python Ecosystem for Divergent Viral Pathogen Discovery, Circular Genome Curation, Replicate-Aware Quasispecies Profiling, Seam-Aware GenBank Annotation, and Ecological Co-Occurrence Network Inference.**

Developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

---

## Modular Package Ecosystem

`Virome-Suite` is distributed as **5 standalone, independently installable domain packages** plus **1 unified branching orchestrator**:

| PyPI Package | CLI Command | Primary Scientific Capability |
| :--- | :--- | :--- |
| **`viro-drift`** | `viro-drift` | Evidence-preserving host read partitioning, 6-frame translated catalytic seed baiting (`RT`, `RdRp`, `Rep`), and cycle-safe local micro-assembly. |
| **`viral-genome-tools`** | `viral-genome-tools` | Terminal Inverted Repeat (TIR) veto, $1.1\times\text{--}3.0\times$ tandem concatemer monomerization, ICTV landmark rotation (`tRNA-Met` PBS / stem-loop), and read-geometry circular junction scoring ($E_{\text{junction}}$). |
| **`quasispecies-profiler`** | `quasispecies-profiler` | Hierarchical Empirical Bayes Beta-Binomial intra-host variant calling across replicates, null parametric bootstrap FDR calibration, 95% profile likelihood CIs, and diversity metrics ($\pi, H, d_N/d_S$). |
| **`viralbio-reporting`** | `viralbio-reporting` | 360° rotation-invariant $2\times$ circular wraparound ORF annotation, $(+)/(-)$ origin-spanning 5-column NCBI `.tbl` generation, SHA-256 provenance chains, and `table2asn` validation. |
| **`viro-networks`** | `viro-networks` | Centered Log-Ratio (CLR) normalization, 3-estimator bootstrap consensus networks (SparCC + shrinkage partial correlation + rank), Sloan Neutral Community Model fitting, and topological hub perturbation. |
| **`virome-suite`** | `virome-suite` | Unified branching orchestrator routing circular DNA and linear RNA viral genomes across all five modules with preflight environment verification (`virome-suite check-env`). |

---

## Installation

Install the complete ecosystem in one command:

```bash
pip install virome-suite
```

Or install any individual domain tool standalone:

```bash
pip install viro-drift
pip install viral-genome-tools
pip install quasispecies-profiler
pip install viralbio-reporting
pip install viro-networks
```

For optional external bioinformatics binaries (`bowtie2`, `spades.py`, `megahit`, `diamond`, `hmmer`, `samtools`, `lofreq`, `table2asn`), create the pinned Conda environment:

```bash
conda env create -f environment.yml
conda activate virome_suite
virome-suite check-env
```

---

## Quick Start

```bash
# 1. Verify Python scientific stack and external binaries
virome-suite check-env

# 2. Run the branching end-to-end pipeline
virome-suite run \
  --r1 sample_R1.fasta \
  --r2 sample_R2.fasta \
  --seeds core_viral_seeds.faa \
  --output results/
```

---

## Maintainer & Contact

* **Maintainer:** Dr. Shweta Singh (`shwetabac@gmail.com`)
* **Affiliation:** Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute
* **Repository:** [https://github.com/molecularplantpathology/Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite)
