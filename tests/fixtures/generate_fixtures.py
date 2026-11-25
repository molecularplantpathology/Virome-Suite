"""Deterministic Generator for Multi-Virus Ground-Truth Validation Fixtures (`tests/fixtures/`).

Creates proven multi-family circular and linear viral reference sequences, tandem concatemers,
linear terminal-repeat negative controls, conserved catalytic protein baits, 3-replicate
variant count tables, and multi-sample cohort abundance matrices.
"""

from __future__ import annotations

import csv
import hashlib
import pathlib
import random


def _make_coding_seq(rng: random.Random, num_codons: int) -> str:
    sense_codons = [
        "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "AAT", "AAC",
        "GAT", "GAC", "TGT", "TGC", "CAA", "CAG", "GAA", "GAG",
        "GGT", "GGC", "CAT", "CAC", "ATT", "ATC", "CTT", "CTC",
        "TTA", "TTG", "AAA", "AAG", "ATG", "TTT", "TTC", "CCT",
        "CCC", "TCT", "TCC", "ACT", "ACC", "TGG", "TAT", "TAC",
        "GTT", "GTC", "GTA", "GTG",
    ]
    body = "".join(rng.choice(sense_codons) for _ in range(num_codons - 2))
    return "ATG" + body + "TAA"


def generate_all_fixtures(fixtures_dir: pathlib.Path) -> dict:
    fixtures_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261117)

    # 1. Proven Circular Virus Monomers (3 distinct families: Caulimoviridae, Geminiviridae, Circoviridae)
    badna_pbs = "TGGTATCAGAGC"  # 12 bp
    badna_orf1 = _make_coding_seq(rng, 220)  # 660 bp
    badna_orf2 = _make_coding_seq(rng, 260)  # 780 bp
    # Make a seam-crossing ORF3 when rotated, or placed near the 3' end wrapping into 5'
    badna_spacer = "".join(rng.choice("ACGT") for _ in range(249))
    badna_orf3 = _make_coding_seq(rng, 233)  # 699 bp
    badna_monomer = (badna_pbs + badna_spacer + badna_orf1 + badna_orf2 + badna_orf3)[:2400]

    gemini_ori = "TAATATTAC"  # 9 bp stem-loop nonanucleotide
    gemini_rep = _make_coding_seq(rng, 350)  # 1050 bp
    gemini_cp = _make_coding_seq(rng, 280)   # 840 bp
    gemini_spacer = "".join(rng.choice("ACGT") for _ in range(801))
    gemini_monomer = (gemini_ori + gemini_spacer + gemini_rep + gemini_cp)[:2700]

    circo_ori = "TAGTATTAC"  # 9 bp circovirus origin
    circo_rep = _make_coding_seq(rng, 300)  # 900 bp
    circo_cap = _make_coding_seq(rng, 220)  # 660 bp
    circo_spacer = "".join(rng.choice("ACGT") for _ in range(231))
    circo_monomer = (circo_ori + circo_spacer + circo_rep + circo_cap)[:1800]

    circ_fasta = fixtures_dir / "proven_circular_monomers.fasta"
    with circ_fasta.open("w", encoding="utf-8") as f:
        f.write(">Badnavirus_Ref_Monomer family=Caulimoviridae landmark=tRNA-Met_PBS length=2400\n")
        f.write(badna_monomer + "\n")
        f.write(">Begomovirus_Ref_Monomer family=Geminiviridae landmark=StemLoop_Nonanucleotide length=2700\n")
        f.write(gemini_monomer + "\n")
        f.write(">Circovirus_Ref_Monomer family=Circoviridae landmark=Origin_Nonanucleotide length=1800\n")
        f.write(circo_monomer + "\n")

    # 2. Stretched Concatemers (Rotated + duplicated 1.8x to 2.25x)
    def make_concatemer(seq: str, rot: int, mult: float) -> str:
        rotated = seq[rot:] + seq[:rot]
        target_len = int(round(len(seq) * mult))
        repeats = (target_len // len(seq)) + 2
        return (rotated * repeats)[:target_len]

    concat_fasta = fixtures_dir / "concatemer_test_panel.fasta"
    with concat_fasta.open("w", encoding="utf-8") as f:
        c_badna = make_concatemer(badna_monomer, rot=450, mult=2.25)
        f.write(f">Badnavirus_Concatemer_2.25x true_monomer_len=2400 concat_len={len(c_badna)}\n")
        f.write(c_badna + "\n")
        c_gemini = make_concatemer(gemini_monomer, rot=600, mult=1.80)
        f.write(f">Begomovirus_Concatemer_1.80x true_monomer_len=2700 concat_len={len(c_gemini)}\n")
        f.write(c_gemini + "\n")
        c_circo = make_concatemer(circo_monomer, rot=300, mult=2.00)
        f.write(f">Circovirus_Concatemer_2.00x true_monomer_len=1800 concat_len={len(c_circo)}\n")
        f.write(c_circo + "\n")

    # 3. Proven Linear Viruses & Terminal-Repeat Negative Controls
    poace_seq = "".join(rng.choice("ACGT") for _ in range(150)) + _make_coding_seq(rng, 900) + "".join(rng.choice("ACGT") for _ in range(150))
    polero_seq = "".join(rng.choice("ACGT") for _ in range(120)) + _make_coding_seq(rng, 720) + "".join(rng.choice("ACGT") for _ in range(120))
    comp = str.maketrans("ACGT", "TGCA")
    tir_arm = "".join(rng.choice("ACGT") for _ in range(180))
    tir_rc = tir_arm.translate(comp)[::-1]
    linear_tir = tir_arm + _make_coding_seq(rng, 700) + tir_rc

    lin_fasta = fixtures_dir / "proven_linear_controls.fasta"
    with lin_fasta.open("w", encoding="utf-8") as f:
        f.write(f">Poacevirus_Linear_Ref family=Potyviridae topology=linear length={len(poace_seq)}\n")
        f.write(poace_seq + "\n")
        f.write(f">Polerovirus_Linear_Ref family=Solemoviridae topology=linear length={len(polero_seq)}\n")
        f.write(polero_seq + "\n")
        f.write(f">Linear_TIR_Negative_Control topology=linear_with_TIR length={len(linear_tir)}\n")
        f.write(linear_tir + "\n")

    # 4. Conserved Catalytic Protein Seeds (FAA)
    seeds_faa = fixtures_dir / "core_viral_seeds.faa"
    with seeds_faa.open("w", encoding="utf-8") as f:
        f.write(">Badnavirus_RT_Core domain=Reverse_Transcriptase family=Caulimoviridae\n")
        f.write("MDDILVYSPGSKAIKFAWKKLPVTVKIAQKYGVDPRKLLQYVDDILIFSPTSKELDQLLK\n")
        f.write(">Potyviridae_RdRp_Core domain=RdRp_GDD family=Potyviridae\n")
        f.write("SGQPSTVVDNTLMVIIAMTYALRKLGVPEDILVGDDNLIYSLPERVREMVKTFGFKVT\n")
        f.write(">Geminiviridae_Rep_Core domain=RCR_Rep family=Geminiviridae\n")
        f.write("MPRLFKIYAKNYFLTYPNCSLSKEEALSQLKKLPTNIKFIKICRELHENGEPHLHVLIQF\n")

    # 5. Replicate Variant Calling Fixture (3 Replicates: True iSNVs vs PCR Jackpot Noise)
    var_tsv = fixtures_dir / "replicate_variants_fixture.tsv"
    var_rows = [
        # True spiked variants (consistent across all 3 replicates)
        {"pos": 150, "ref": "A", "alt": "G", "truth": "true_isnv", "expected_vaf": 0.045, "rep1_k": 46, "rep1_n": 1000, "rep2_k": 44, "rep2_n": 1000, "rep3_k": 45, "rep3_n": 1000},
        {"pos": 320, "ref": "C", "alt": "T", "truth": "true_isnv", "expected_vaf": 0.020, "rep1_k": 210, "rep1_n": 10000, "rep2_k": 195, "rep2_n": 10000, "rep3_k": 202, "rep3_n": 10000},
        {"pos": 580, "ref": "G", "alt": "A", "truth": "true_isnv", "expected_vaf": 0.010, "rep1_k": 105, "rep1_n": 10000, "rep2_k": 98, "rep2_n": 10000, "rep3_k": 101, "rep3_n": 10000},
        {"pos": 910, "ref": "T", "alt": "C", "truth": "true_isnv", "expected_vaf": 0.006, "rep1_k": 62, "rep1_n": 10000, "rep2_k": 58, "rep2_n": 10000, "rep3_k": 60, "rep3_n": 10000},
        # Single-replicate PCR jackpot noise (high in 1 rep, background in others)
        {"pos": 1120, "ref": "A", "alt": "T", "truth": "pcr_jackpot_artifact", "expected_vaf": 0.0, "rep1_k": 180, "rep1_n": 10000, "rep2_k": 2, "rep2_n": 10000, "rep3_k": 1, "rep3_n": 10000},
        {"pos": 1450, "ref": "C", "alt": "A", "truth": "pcr_jackpot_artifact", "expected_vaf": 0.0, "rep1_k": 1, "rep1_n": 10000, "rep2_k": 150, "rep2_n": 10000, "rep3_k": 3, "rep3_n": 10000},
        # Monomorphic background error sites
        {"pos": 1700, "ref": "G", "alt": "T", "truth": "monomorphic_bg", "expected_vaf": 0.0, "rep1_k": 3, "rep1_n": 10000, "rep2_k": 2, "rep2_n": 10000, "rep3_k": 4, "rep3_n": 10000},
        {"pos": 1950, "ref": "T", "alt": "G", "truth": "monomorphic_bg", "expected_vaf": 0.0, "rep1_k": 2, "rep1_n": 10000, "rep2_k": 3, "rep2_n": 10000, "rep3_k": 2, "rep3_n": 10000},
    ]
    with var_tsv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(var_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(var_rows)

    # 6. Multi-Sample Cohort Count Matrix for viro_networks (24 samples x 6 viral taxa)
    counts_tsv = fixtures_dir / "virome_cohort_counts.tsv"
    taxa = [
        "Badnavirus_A",
        "Badnavirus_B",
        "Poacevirus_SCSMV",
        "Polerovirus_SCYLV",
        "Begomovirus_C",
        "Mastrevirus_D",
    ]
    with counts_tsv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id"] + taxa)
        for s_idx in range(1, 25):
            # Planted co-occurrence between Poacevirus_SCSMV and Polerovirus_SCYLV
            base_co = rng.randint(800, 4000) if s_idx % 2 == 0 else 0
            row = [
                f"SAMPLE_{s_idx:02d}",
                rng.randint(100, 2500),
                rng.randint(50, 1200) if s_idx % 3 != 0 else 0,
                base_co + (rng.randint(50, 200) if base_co > 0 else 0),
                int(base_co * 0.85) + (rng.randint(40, 150) if base_co > 0 else 0),
                rng.randint(100, 800) if s_idx % 4 != 0 else 0,
                rng.randint(50, 600) if s_idx % 5 != 0 else 0,
            ]
            writer.writerow(row)

    summary = {}
    for p in sorted(fixtures_dir.glob("*")):
        if p.is_file() and p.name != "generate_fixtures.py":
            summary[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    return summary


if __name__ == "__main__":
    here = pathlib.Path(__file__).resolve().parent
    res = generate_all_fixtures(here)
    for k, v in res.items():
        print(f"Fixture {k:35s} SHA256={v}")
