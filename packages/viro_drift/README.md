# viro-drift

**Adaptive Divergent Viral Pathogen Discovery via 6-Frame Catalytic Seed Baiting & Local Micro-Assembly**

Part of the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

## Key Features
- **Evidence-Preserving Host Read Partitioning:** Separates host-concordant reads, ambiguous host-viral chimeric bridges (preserved for downstream integration checks), and unmapped viral candidates.
- **Low-Complexity & Homopolymer Shield:** Sliding-window nucleotide Shannon entropy masking ($H < 1.20\text{ bits}$) with strict stop/ambiguous residue exclusion.
- **6-Frame Collinear Catalytic Baiting:** Recruits divergent viral reads (15%–35% nucleotide divergence) using reduced-alphabet amino acid $k$-mers (`kmer_seed` self-contained engine) or external `DIAMOND blastx` / `HMMER hmmsearch`.
- **Two-Phase Cycle-Safe Micro-Assembly:** Strand-normalized iterative paired recruitment and local de Bruijn assembly supporting both linear unitigs and closed circular viral rings.

## Installation

```bash
pip install viro-drift
```

## Quick Start

```bash
viro-drift run --r1 unmapped_R1.fasta --r2 unmapped_R2.fasta --seeds viral_core_seeds.faa -o out_drift/
```
