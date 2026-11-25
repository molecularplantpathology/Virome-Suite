"""Package 1: viro_drift — Adaptive Divergent Viral Pathogen Discovery."""

from viro_drift.host_partition import PartitionedReads, build_host_kmer_index, partition_host_reads
from viro_drift.micro_assembler import AssembledContig, assemble_viral_contigs
from viro_drift.recruitment import RecruitmentResult, iterative_recruit_reads
from viro_drift.translation_bait import (
    ReducedAlphabetSeedIndex,
    bait_reads_with_seeds,
    mask_low_complexity_nt,
    shannon_entropy_nt,
    translate_6frames,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "AssembledContig",
    "PartitionedReads",
    "RecruitmentResult",
    "ReducedAlphabetSeedIndex",
    "assemble_viral_contigs",
    "bait_reads_with_seeds",
    "build_host_kmer_index",
    "iterative_recruit_reads",
    "mask_low_complexity_nt",
    "partition_host_reads",
    "shannon_entropy_nt",
    "translate_6frames",
]
