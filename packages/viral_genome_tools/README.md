# viral-genome-tools

**Circular Viral Monomer Curation, Multi-Family Landmark Rotation & Physical Circular Junction Scoring**

Part of the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

## Key Features
- **Terminal Inverted Repeat (TIR) Veto:** Prevents false circularization of linear viral genomes bearing inverted terminal repeats.
- **Indel-Tolerant Concatemer Monomerization:** Resolves $1.1\times\text{--}3.0\times$ rolling-circle and pararetroviral tandem concatemers to exact $1.0\times$ monomers via Revolution-2 anchor banded alignment.
- **Multi-Family ICTV Landmark Canonicalization:** Two-pass priority matching across *Caulimoviridae* (`tRNA-Met` PBS), *Geminiviridae* (`TAATATTAC` stem-loop), *Circoviridae*, and *Nanoviridae*, with deterministic ORF-start fallback for novel circular viruses.
- **Read-Geometry Circular Junction Score ($E_{\text{junction}}$) & Orthogonal $P_{\text{flank}}$:** Continuous depth-scaled junction scoring with a $<10\times$ low-depth guard and explicit `unresolved_topology` abstention tier.

## Installation

```bash
pip install viral-genome-tools
```

## Quick Start

```bash
viral-genome-tools curate --assembly raw_contigs.fasta --r1 reads_R1.fasta --r2 reads_R2.fasta -o out_curation/
```
