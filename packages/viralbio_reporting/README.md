# viralbio-reporting

[![PyPI](https://img.shields.io/pypi/v/viralbio-reporting.svg)](https://pypi.org/project/viralbio-reporting/)
[![CI & Coverage](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml/badge.svg)](https://github.com/molecularplantpathology/Virome-Suite/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**360° Rotation-Invariant Circular Seam-Spanning ORF Annotation, NCBI 5-Column `.tbl`/`.fsa` Submission Packaging, and SHA-256 Provenance Chains.**

`viralbio-reporting` is a standalone package within the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

---

## Algorithmic Pipeline

1. **Modulo-$L$ Rotation-Invariant Circular ORF Discovery (`viralbio_reporting.circular_annotator`):**
   - Evaluates circular templates on a virtual $2\times$ sequence $S_{2\times} = S + S$ (and linear templates on $1\times$ $S$), restricting starts to $0 \le s_0 < L$ and bounding single-revolution coding lengths to $3 \le e_0 - s_0 \le L$ (including the terminal 3-bp stop codon).
   - Deduplicates candidate ORFs by their modulo-$L$ biological stop coordinate `(strand, (e_0 - 1) % L)`, retaining the longest single-revolution frame per stop codon. This mathematically guarantees **100% identical protein predictions across any 360° circular cut point**, even if the cut lands between a canonical start `ATG` and an internal Methionine `ATG`.
2. **Closed-Form 1-Based NCBI Compound Interval Projection:**
   - **Forward (`+`) Origin-Spanning ($e_0 > L$):** Projects to two tab-separated interval rows `[s_0 + 1, L]` and `[1, e_0 - L]`.
   - **Reverse (`-`) Origin-Spanning ($e_0 > L$ on $\mathrm{rc}(S + S)$):** Projects to descending 1-based rows `[L - s_0, 1]` and `[L, 2L - e_0 + 1]`.
   - **$5'$-Partial Linear CDS:** Anchors interval start at `<1` (or `<L` on `-` strand) and encodes reading-frame phase strictly via `codon_start` $\in \{1, 2, 3\}$ without double-shifting coordinates.
3. **Biological & Syntax Validation + SHA-256 Provenance (`viralbio_reporting.ncbi_packager`):**
   Reconstructs every CDS directly from its exported 1-based `.tbl` intervals to verify start/stop codon integrity and exact amino acid translation, runs internal 5-column `.tbl` syntax checks (plus headless NCBI `table2asn -V vb` whenever installed), and writes a tamper-evident SHA-256 hash chain (`provenance_chain.json`).

---

## Installation

```bash
pip install viralbio-reporting
```

---

## Command-Line Usage (`viralbio-reporting`)

```bash
viralbio-reporting annotate \
  --fasta curated_monomers.fasta \
  --topology circular \
  --min-aa 60 \
  --organism "Sugarcane bacilliform virus" \
  --host "Saccharum officinarum" \
  --country "India" \
  --output genbank_submission/
```

### CLI Options (`viralbio-reporting annotate`)

| Argument | Required | Default | Description |
| :--- | :---: | :---: | :--- |
| `--fasta` | Yes | — | Input curated circular or linear genome FASTA (single or multi-record) |
| `--topology` | No | `circular` | Sequence topology: `circular` (enables $2\times$ seam wraparound) or `linear` |
| `--min-aa` | No | `60` | Minimum protein length in amino acids |
| `--organism` | No | `"Unclassified virus"` | Organism source modifier written to `.fsa` and `source.tsv` |
| `--host` | No | `"Plant host"` | Host species qualifier written to `source.tsv` |
| `--country` | No | `"Unspecified"` | Geographic origin qualifier written to `source.tsv` |
| `-o`, `--output` | Yes | — | Output directory for GenBank submission files and validation summary |

---

## Python API Usage

```python
import pathlib
from viralbio_reporting.circular_annotator import annotate_genome_orfs, validate_cds_feature
from viralbio_reporting.ncbi_packager import build_ncbi_submission_package

orfs = annotate_genome_orfs(circular_seq, topology="circular", min_aa=60)
for feat in orfs:
    is_ok, msg = validate_cds_feature(circular_seq, feat)
    print(feat.feature_id, feat.strand, feat.intervals_1based, f"spans_origin={feat.spans_origin}", is_ok)

report = build_ncbi_submission_package(
    seq_id="Isolate_ICAR_01",
    sequence=circular_seq,
    features=orfs,
    output_dir=pathlib.Path("ncbi_pkg/"),
    topology="circular",
    organism="Sugarcane bacilliform virus",
    host="Saccharum officinarum",
    country="India",
)
```

---

## Output Files

- **`<SeqID>.fsa`:** NCBI GenBank submission FASTA with `[organism=...]`, `[topology=circular|linear]`, and `[gcode=1]` modifiers.
- **`<SeqID>.tbl`:** Strict 5-column tab-delimited NCBI feature table supporting origin-wrapping compound CDS intervals and partial linear qualifiers.
- **`source.tsv`:** Tab-delimited source modifier table (`Sequence_ID`, `Organism`, `Host`, `Country`).
- **`provenance_chain.json` & `annotation_summary.json`:** Cryptographic SHA-256 transformation chain and biological/syntax validation report.
