"""Preflight Environment & Dependency Health Checker for Virome_Suite (`virome-suite check-env`)."""

from __future__ import annotations

import importlib.metadata
import pathlib
import platform
import shutil
import sys
from typing import Any, Dict, List


REQUIRED_PYTHON_LIBRARIES = [
    ("numpy", "Core numerical arrays & matrix algebra"),
    ("scipy", "Beta-Binomial optimization & statistical distributions"),
    ("pandas", "Tabular manifest & abundance matrix processing"),
    ("Bio", "Biopython sequence & GenBank feature handling"),
]

EXTERNAL_BINARIES = [
    ("bowtie2", "viro_drift / host_partition", "Host read partitioning & alignment"),
    ("spades.py", "viro_drift", "Local multi-k-mer viral micro-assembly"),
    ("megahit", "viro_drift", "Metagenomic de Bruijn graph assembly"),
    ("diamond", "viro_drift", "6-frame translated protein seed baiting (blastx)"),
    ("hmmsearch", "viro_drift", "Profile HMM catalytic domain baiting"),
    ("samtools", "viral_genome_tools / quasispecies_profiler", "BAM sorting & depth profiling"),
    ("lofreq", "quasispecies_profiler", "Baseline Poisson-binomial variant comparator"),
    ("table2asn", "viralbio_reporting", "Official NCBI GenBank submission validator"),
]


def check_python_package(pkg_name: str) -> Dict[str, Any]:
    dist_lookup = {"Bio": "biopython", "yaml": "PyYAML"}.get(pkg_name, pkg_name)
    try:
        __import__(pkg_name)
        try:
            ver = importlib.metadata.version(dist_lookup)
        except importlib.metadata.PackageNotFoundError:
            ver = "installed"
        return {"name": dist_lookup, "installed": True, "version": ver}
    except ImportError:
        return {"name": dist_lookup, "installed": False, "version": "missing"}


def check_external_binary(bin_name: str, module: str, role: str) -> Dict[str, Any]:
    exe_path = shutil.which(bin_name)
    return {
        "binary": bin_name,
        "module": module,
        "role": role,
        "available": exe_path is not None,
        "path": exe_path or "not_in_PATH (install via `conda env create -f environment.yml`)",
    }


def run_preflight_check(root_dir: pathlib.Path | None = None) -> Dict[str, Any]:
    if root_dir is None:
        root_dir = pathlib.Path(__file__).resolve().parent.parent

    py_results = [check_python_package(pkg) for pkg, _ in REQUIRED_PYTHON_LIBRARIES]
    bin_results = [check_external_binary(b, m, r) for b, m, r in EXTERNAL_BINARIES]

    manifest_path = root_dir / "scripts_manifest.tsv"
    protocol_path = root_dir / "benchmarks" / "benchmark_protocol.yaml"
    fixtures_dir = root_dir / "tests" / "fixtures"

    return {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "python_dependencies": py_results,
        "external_binaries": bin_results,
        "workspace_artifacts": {
            "scripts_manifest_present": manifest_path.exists(),
            "benchmark_protocol_present": protocol_path.exists(),
            "test_fixtures_present": fixtures_dir.exists(),
        },
        "core_python_ready": all(item["installed"] for item in py_results[:3]),
    }


def format_preflight_report(report: Dict[str, Any]) -> str:
    lines: List[str] = []
    lines.append("=" * 88)
    lines.append(" VIROME_SUITE PREFLIGHT ENVIRONMENT & HEALTH CHECK (`virome-suite check-env`)")
    lines.append("=" * 88)
    lines.append(f" Python Runtime : {report['python_version']} ({report['platform']})")
    lines.append("-" * 88)
    lines.append(" [1] Python Scientific Stack:")
    for item in report["python_dependencies"]:
        status = "[OK]" if item["installed"] else "[MISSING]"
        lines.append(f"   {status:10s} {item['name']:15s} (version: {item['version']})")

    lines.append("-" * 88)
    lines.append(" [2] External Bioinformatics Executables (Managed via environment.yml):")
    for item in report["external_binaries"]:
        status = "[FOUND]" if item["available"] else "[OPTIONAL/CONDA]"
        lines.append(f"   {status:16s} {item['binary']:12s} -> {item['module']} ({item['role']})")

    lines.append("-" * 88)
    lines.append(" [3] Gate 1 / Mark 1 Foundation Artifacts:")
    art = report["workspace_artifacts"]
    lines.append(
        f"   {'[OK]' if art['scripts_manifest_present'] else '[MISSING]':10s} scripts_manifest.tsv (263-Script Audit Catalog)"
    )
    lines.append(
        f"   {'[OK]' if art['benchmark_protocol_present'] else '[MISSING]':10s} benchmarks/benchmark_protocol.yaml (Locked Test Rules)"
    )
    lines.append(
        f"   {'[OK]' if art['test_fixtures_present'] else '[MISSING]':10s} tests/fixtures/ (Multi-Virus Proven Test Genomes)"
    )
    lines.append("=" * 88)
    return "\n".join(lines)
