# viro-drift

[![PyPI](https://img.shields.io/pypi/v/viro-drift.svg)](https://pypi.org/project/viro-drift/)
[![CI & Coverage](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml/badge.svg)](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Divergent Sequence Rescue via Evidence-Preserving Host Partitioning, 6-Frame Translated Catalytic Seed Baiting, and Two-Phase Cycle-Safe Micro-Assembly.**

`viro-drift` is a standalone, zero-internal-dependency package within the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

---

## Algorithmic Pipeline

1. **Evidence-Preserving Host Read Partitioning (`viro_drift.host_partition`):**
   Partitions paired-end or single-end reads against an optional host reference $k$-mer index ($k = 21$) into three mutually exclusive streams:
   - `host_concordant`: Pure host read pairs filtered prior to assembly.
   - `host_spanning_ambiguous`: Discordant pairs or chimeric reads carrying both host and non-host segments (preserved for downstream integration screening rather than silently discarded).
   - `viral_candidate`: Unmapped candidate read pairs forwarded to catalytic baiting.
2. **Low-Complexity Masking & 6-Frame Collinear Seed Baiting (`viro_drift.translation_bait`):**
   - Masks sliding 24-bp windows with nucleotide Shannon entropy $H(w) < 1.20\text{ bits}$ to `'N'`.
   - Translates all 6 reading frames (`+0, +1, +2, -0, -1, -2`) and splits peptides strictly on stop (`'*'`) and ambiguous (`'X'`) residues so masked regions never index or cross-match.
   - Supports a built-in pure-Python collinear reduced-alphabet amino acid matcher (`--mode kmer_seed`, 11 physicochemical equivalence groups, diagonal consistency $|(q_2 - q_1) - (p_2 - p_1)| \le 3\text{ aa}$) as well as external `--mode diamond` (`DIAMOND blastx`) and `--mode hmmer` (`hmmsearch`).
3. **Strand-Aware Canonical Recruitment (`viro_drift.recruitment`):**
   Iteratively rescues paired mates and unmapped reads sharing high-entropy ($H \ge 1.30\text{ bits}$) terminal canonical 21-mers $\min(\kappa, \overline{\kappa})$ while blocking host-masked $k$-mers and halting automatically on growth plateau ($<1\%$) or repeat branching ($>4.0$).
4. **Two-Phase Cycle-Safe De Bruijn Micro-Assembly (`viro_drift.micro_assembler`):**
   Builds strand-normalized local de Bruijn graphs across an odd-$k$ cascade (`21, 33, 55`), extracting linear/branching unitigs in Phase 1 and closed circular rings (`in_degree = out_degree = 1`) in Phase 2 with rotation-invariant reverse-complement deduplication.

---

## Installation

```bash
pip install viro-drift
```

---

## Command-Line Usage (`viro-drift`)

```bash
viro-drift run \
  --r1 sample_R1.fastq \
  --r2 sample_R2.fastq \
  --seeds core_viral_seeds.faa \
  --host host_genome.fasta \
  --mode kmer_seed \
  --max-rounds 5 \
  --output drift_results/
```

### CLI Options (`viro-drift run`)

| Argument | Required | Default | Description |
| :--- | :---: | :---: | :--- |
| `--r1` | Yes | — | Forward reads in FASTA or FASTQ format |
| `--r2` | No | `None` | Reverse paired-end reads in FASTA or FASTQ format |
| `--seeds` | Yes | — | Conserved protein domain seeds in FASTA (`.faa`) format (`RdRp`, `RT`, `Rep`, etc.) |
| `--host` | No | `None` | Optional host genome FASTA for 21-mer host partitioning and decoy masking |
| `--mode` | No | `kmer_seed` | Baiting backend: `kmer_seed` (pure Python), `diamond`, or `hmmer` |
| `--max-rounds` | No | `5` | Maximum iterative canonical $k$-mer read recruitment rounds |
| `-o`, `--output` | Yes | — | Output directory for assembled contigs and JSON run summary |

**Exit Codes:** `0` = success (including clean negative controls with `NO_VIRAL_CANDIDATES`), `2` = missing/invalid input files, `3` = requested external binary (`diamond`/`hmmsearch`) not found in `PATH`.

---

## Python API Usage

```python
from viro_drift.host_partition import build_host_kmer_index, partition_host_reads
from viro_drift.translation_bait import bait_reads_with_seeds
from viro_drift.recruitment import iterative_recruit_reads
from viro_drift.micro_assembler import assemble_viral_contigs

host_kmers = build_host_kmer_index([host_fasta_seq], k=21)
partition = partition_host_reads(read_pairs, host_kmers, k=21)

seed_ids, hit_domains = bait_reads_with_seeds(
    partition.viral_candidate,
    seeds={"Badnavirus_RT": "MDDILVYSPGSKAIKFAWKKLPVTVKIAQKYGVDPRKLLQYVDDILIFSPTSKELDQLLK"},
    mode="kmer_seed",
)

recruitment = iterative_recruit_reads(
    partition.viral_candidate,
    initial_seed_ids=seed_ids,
    host_kmer_mask=host_kmers,
    max_rounds=5,
)
contigs = assemble_viral_contigs(recruitment.recruited_pairs)
for c in contigs:
    print(c.contig_id, len(c.sequence), c.mean_kmer_depth, c.is_circular_cycle)
```

---

## Output Files

- **`assembled_contigs.fasta`:** Reconstructed candidate contigs annotated with `length`, `depth`, and `circular_cycle` header attributes (empty file if `NO_VIRAL_CANDIDATES`).
- **`viro_drift_summary.json`:** Structured JSON report containing `status`, `partition` read counts, `initial_seed_reads`, `matched_seed_domains`, `recruited_pairs`, `rounds_executed`, `stop_reason`, `filtered_host_decoys_count`, and `assembled_contigs`.
