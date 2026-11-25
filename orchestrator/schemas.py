"""Shared Data Contracts, Read-Geometry Equations, and Provenance Schemas for Virome_Suite."""

from __future__ import annotations

import dataclasses
import datetime
import enum
import hashlib
import json
import math
import pathlib
from typing import Any, Dict, List, Optional


SCHEMA_VERSION = "1.0.0"


class TopologyState(str, enum.Enum):
    """Explicit genome morphology and topology classification states."""

    CIRCULAR_CANDIDATE = "circular_candidate"
    LINEAR_CANDIDATE = "linear_candidate"
    SEGMENTED_COMPONENT = "segmented_component"
    INTEGRATED_CANDIDATE = "integrated_candidate"
    MIXED_OR_CONFLICTED = "mixed_or_conflicted"
    UNRESOLVED_TOPOLOGY = "unresolved_topology"


@dataclasses.dataclass(frozen=True)
class LandmarkMotif:
    """Curated biological landmark for circular coordinate-1 canonicalization."""

    family: str
    genus: str
    landmark_name: str
    motif_seq: str
    description: str


# Multi-Family Circular Viral Landmark Registry (Generalizing beyond Badnaviruses)
MULTI_FAMILY_LANDMARK_REGISTRY: List[LandmarkMotif] = [
    LandmarkMotif(
        family="Caulimoviridae",
        genus="Badnavirus",
        landmark_name="tRNA-Met_PBS",
        motif_seq="TGGTATCAGAGC",
        description="Cytosolic plant tRNA-Met primer binding site (ICTV coordinate 1 anchor)",
    ),
    LandmarkMotif(
        family="Caulimoviridae",
        genus="Caulimovirus",
        landmark_name="tRNA-Met_PBS_Caulimo",
        motif_seq="TGGTATCAGAGCC",
        description="Caulimovirus 35S pregenomic RNA tRNA-Met primer binding site",
    ),
    LandmarkMotif(
        family="Geminiviridae",
        genus="Begomovirus",
        landmark_name="StemLoop_Nonanucleotide",
        motif_seq="TAATATTAC",
        description="Conserved rolling-circle origin stem-loop nonanucleotide (TAATATT/AC)",
    ),
    LandmarkMotif(
        family="Geminiviridae",
        genus="Mastrevirus",
        landmark_name="StemLoop_Nonanucleotide",
        motif_seq="TAATATTAC",
        description="Conserved rolling-circle origin stem-loop nonanucleotide (TAATATT/AC)",
    ),
    LandmarkMotif(
        family="Circoviridae",
        genus="Circovirus",
        landmark_name="Origin_Nonanucleotide",
        motif_seq="TAGTATTAC",
        description="Conserved circoviral rolling-circle replication origin nonanucleotide",
    ),
    LandmarkMotif(
        family="Nanoviridae",
        genus="Nanovirus",
        landmark_name="Origin_Nonanucleotide",
        motif_seq="TAGTATTAC",
        description="Conserved multipartite circular ssDNA origin nonanucleotide",
    ),
]


def expected_split_reads(
    d_median: float,
    read_length: int = 150,
    k_min_anchor: int = 20,
) -> float:
    """Compute physically expected split reads crossing a circular seam.

    E[N_split] = D_median * (1 - 2 * k_min / L_r)
    """
    if d_median <= 0 or read_length <= 2 * k_min_anchor:
        return 0.0
    return float(d_median) * (1.0 - (2.0 * float(k_min_anchor)) / float(read_length))


def expected_spanning_pairs(
    d_median: float,
    read_length: int = 150,
    mean_insert_size: float = 350.0,
) -> float:
    """Compute physically expected paired-end fragments bridging a circular seam.

    When per-base depth D_median counts aligned bases from both mates,
    fragment start density is D_median / (2 * L_r), and the eligible insert window
    is max(0, mu_ins - 2 * L_r).
    """
    if d_median <= 0 or read_length <= 0:
        return 0.0
    window = max(0.0, float(mean_insert_size) - 2.0 * float(read_length))
    return (float(d_median) / (2.0 * float(read_length))) * window


def compute_coverage_uniformity(
    d_origin: float,
    d_median: float,
    d_min: float = 10.0,
    epsilon: float = 1.0,
) -> float:
    """Compute bounded coverage uniformity across the circular origin with low-depth guard."""
    if d_median < d_min:
        return 0.0
    ratio = (float(d_origin) + epsilon) / (float(d_median) + epsilon)
    return math.exp(-abs(math.log(ratio)))


def compute_host_flank_ratio(
    n_split: int,
    n_host_viral_chimera: int,
) -> Optional[float]:
    """Compute orthogonal host-integration evidence ratio P_flank.

    Returns None (Unassessed / Insufficient Data) when zero informative reads exist,
    never confusing absence of evidence with evidence of absence.
    """
    total = int(n_split) + int(n_host_viral_chimera)
    if total <= 0:
        return None
    return float(n_host_viral_chimera) / float(total)


@dataclasses.dataclass
class ContigEvidenceRecord:
    """Versioned inter-package contract for a candidate viral sequence."""

    contig_id: str
    length_bp: int
    d_median: float
    d_origin: float
    n_split_reads: int
    n_spanning_pairs: int
    n_host_viral_chimeras: int
    landmark_score: float  # 0.0, 0.5, or 1.0
    landmark_name: str
    e_junction_score: float
    p_flank: Optional[float]
    topology_state: TopologyState
    raw_sha256: str
    curated_sha256: str
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        data = dataclasses.asdict(self)
        data["topology_state"] = self.topology_state.value
        return data


@dataclasses.dataclass
class ProvenanceEvent:
    """Single transformation record in the append-only SHA-256 provenance chain."""

    step_name: str
    tool_name: str
    tool_version: str
    input_sha256: str
    output_sha256: str
    parameters: Dict[str, Any]
    timestamp_utc: str = dataclasses.field(
        default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat()
    )
    previous_event_hash: str = "GENESIS"

    def event_hash(self) -> str:
        payload = json.dumps(dataclasses.asdict(self), sort_keys=True).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()


class AppendOnlyProvenanceLog:
    """Hash-chained provenance log recording sequence transformations."""

    def __init__(self, log_path: pathlib.Path) -> None:
        self.log_path = log_path
        self.events: List[ProvenanceEvent] = []

    def append_event(
        self,
        step_name: str,
        tool_name: str,
        tool_version: str,
        input_sha256: str,
        output_sha256: str,
        parameters: Optional[Dict[str, Any]] = None,
    ) -> ProvenanceEvent:
        prev_hash = self.events[-1].event_hash() if self.events else "GENESIS"
        ev = ProvenanceEvent(
            step_name=step_name,
            tool_name=tool_name,
            tool_version=tool_version,
            input_sha256=input_sha256,
            output_sha256=output_sha256,
            parameters=parameters or {},
            previous_event_hash=prev_hash,
        )
        self.events.append(ev)
        self._flush()
        return ev

    def _flush(self) -> None:
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": SCHEMA_VERSION,
            "event_count": len(self.events),
            "chain_head_sha256": self.events[-1].event_hash() if self.events else "EMPTY",
            "events": [
                {**dataclasses.asdict(e), "event_sha256": e.event_hash()} for e in self.events
            ],
        }
        self.log_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
