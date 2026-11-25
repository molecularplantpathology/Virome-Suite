# viral-genome-tools

[![PyPI](https://img.shields.io/pypi/v/viral-genome-tools.svg)](https://pypi.org/project/viral-genome-tools/)
[![CI & Coverage](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml/badge.svg)](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Exact $1.0\times$ Circular Concatemer Monomerization, Multi-Family Origin Landmark Rotation, and Read-Geometry Circular Junction Scoring ($E_{\text{junction}}$).**

`viral-genome-tools` is a standalone package within the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

---

## Algorithmic Pipeline

1. **Terminal Inverted Repeat (TIR) Linear Veto (`viral_genome_tools.monomer_curator`):**
   Precedence Rule 0 screens $5'$ prefixes against reverse-complemented $3'$ suffixes ($40\text{--}500\text{ bp}$). Contigs carrying genuine Terminal Inverted Repeats ($\ge 90\%$ identity over $\ge 50\text{ bp}$, e.g. *Adenoviridae*, *Parvoviridae*, linear elements) are immediately locked to `linear_candidate` with zero bases trimmed.
2. **Revolution-2 Anchor Monomerization (`viral_genome_tools.monomer_curator`):**
   Resolves $1.1\times\text{--}3.0\times$ tandem assembly concatemers by anchoring index `0` of Revolution 1 to the 0-based start coordinate $p_2$ of Revolution 2 via banded alignment ($\ge 95\%$ identity, $\ge 40\text{ bp}$ direct terminal overlap). Defining monomer length $P^* = p_2$ recovers the exact biological monomer ($\pm 0\text{ bp}$) even when sequencing indels occur inside the downstream overlap tail.
3. **Two-Pass Priority Landmark Canonicalization (`viral_genome_tools.landmark_rotator`):**
   Standardizes strand orientation to `(+)` and rotates Coordinate 1 to conserved ICTV family landmarks across circular wraparound boundaries:
   - *Caulimoviridae* (*Badnavirus* / *Caulimovirus*): `TGGTATCAGAGC` (`tRNA-Met` primer-binding site)
   - *Geminiviridae* (*Begomovirus* / *Mastrevirus* / *Curtovirus*): `TAATATTAC` (stem-loop nonanucleotide)
   - *Circoviridae* / *Nanoviridae*: `TAGTATTAC` / `CAGTATTAC` (origin nonanucleotide)
   - Unknown circular templates lacking registered motifs fall back deterministically to the longest forward-strand `ATG` reading frame (`canonical_orf_orientation`).
4. **Read-Geometry Circular Junction Score & Orthogonal $P_{\text{flank}}$ (`viral_genome_tools.junction_scorer`):**
   Evaluates observed origin-crossing split reads and spanning mate pairs against analytical library expectations:
   $$\mathbb{E}[N_{\text{split}}] = D_{\text{median}}\left(1 - \frac{2 k_{\min}}{L_r}\right), \qquad \mathbb{E}[N_{\text{span}}] = \frac{D_{\text{median}}}{2 L_r}\max(0, \mu_{\text{ins}} - 2 L_r)$$
   Combines split support, spanning pair support, origin coverage uniformity ($D_{\text{median}} \ge 10\times$ guard), and landmark confidence into $E_{\text{junction}} \in [0, 100]$, paired with orthogonal host-chimera ratio $P_{\text{flank}}$ (`None` when unassessed) to classify contigs into `circular_candidate`, `linear_candidate`, `integrated_candidate`, `mixed_or_conflicted`, or `unresolved_topology`.

---

## Installation

```bash
pip install viral-genome-tools
```

---

## Command-Line Usage (`viral-genome-tools`)

```bash
viral-genome-tools curate \
  --assembly raw_assembled_contigs.fasta \
  --r1 sample_R1.fastq \
  --r2 sample_R2.fastq \
  --host host_reference.fasta \
  --insert-size 350 \
  --output curated_out/
```

### CLI Options (`viral-genome-tools curate`)

| Argument | Required | Default | Description |
| :--- | :---: | :---: | :--- |
| `--assembly` | Yes | — | Input assembled contigs FASTA |
| `--r1` | No | `None` | Forward reads (FASTA/FASTQ) for physical seam read-support profiling |
| `--r2` | No | `None` | Reverse paired-end reads (FASTA/FASTQ) for spanning pair profiling |
| `--host` | No | `None` | Optional host reference FASTA for orthogonal $P_{\text{flank}}$ chimera calculation |
| `--custom-motif` | No | `None` | Optional user-supplied Coordinate-1 nucleotide landmark motif |
| `--insert-size` | No | `350.0` | Mean paired-end library insert size ($\mu_{\text{ins}}$) in bp |
| `-o`, `--output` | Yes | — | Output directory for curated FASTA, TSV table, and JSON summary |

---

## Python API Usage

```python
from viral_genome_tools.monomer_curator import trim_concatemer_to_monomer
from viral_genome_tools.landmark_rotator import rotate_to_landmark
from viral_genome_tools.junction_scorer import compute_e_junction_score, classify_topology

trim_res = trim_concatemer_to_monomer(raw_contig_seq)
if not trim_res.has_terminal_inverted_repeat:
    rot_res = rotate_to_landmark(trim_res.monomer_seq)
    e_score, details = compute_e_junction_score(
        n_split=28, n_spanning_pairs=7, d_origin=39.0, d_median=40.0,
        landmark_score=rot_res.landmark_score, read_length=150, mean_insert_size=350.0,
    )
    decision = classify_topology(e_score, d_median=40.0, p_flank=0.0, has_tir=False)
    print(rot_res.matched_family, len(rot_res.rotated_seq), e_score, decision.value)
```

---

## Output Files & Schema

- **`curated_genomes.fasta`:** Canonical $1.0\times$ circular monomers (and unmodified linear contigs) with topology and $E_{\text{junction}}$ header metadata.
- **`curation_report.tsv` & `curation_summary.json`:** Tabular per-contig audit containing `contig_id`, `original_length`, `curated_length`, `was_trimmed`, `has_tir`, `matched_family`, `landmark_name`, `orientation_mode`, `d_median`, `e_junction_score`, `p_flank`, `topology_decision`, `raw_sha256`, and `curated_sha256`.
