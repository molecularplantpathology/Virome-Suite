# viralbio-reporting

**Seamless Circular Wraparound ORF Annotation, SHA-256 Provenance Auditing & NCBI GenBank Package Builder**

Part of the [Virome-Suite](https://github.com/molecularplantpathology/Virome-Suite) ecosystem developed at the **Molecular Plant Pathology Lab, ICAR-Indian Sugarcane Research Institute**.

## Key Features
- **360° Rotation-Invariant Circular ORF Discovery:** $2\times$ virtual sequence wrap (`Seq + Seq`) with modulo-$L$ biological stop-codon deduplication, recovering identical intact proteins regardless of where Coordinate 1 is cut.
- **Exact $(+)$ and $(-)$ Strand Origin-Spanning Intervals:** Projects forward and reverse-complement seam-crossing ORFs into compliant multi-line 1-based NCBI 5-column `.tbl` intervals.
- **Type-Aware Biological Feature Validation:** Supports complete and $5'/3'$-partial linear CDS features with `codon_start` phase verification.
- **Hash-Chained Provenance & `table2asn` Harness:** Emits SHA-256 sequence provenance logs and ready-to-submit `.fsa`/`.tbl`/`source.tsv` files.

## Installation

```bash
pip install viralbio-reporting
```

## Quick Start

```bash
viralbio-reporting annotate --fasta curated_monomer.fasta --topology circular -o out_genbank/
```
