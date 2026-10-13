"""Probe explicitly selected browser converters without importing scientific callers."""
from __future__ import annotations

import importlib.metadata
import hashlib
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile

from .browser_common import require, sha256
from .browser_export import _binary_check


def probe():
    """Return converter identity/API/smoke evidence or raise for an invalid prefix.

Call this only for an explicitly selected conversion environment. A run with no
conversion prefix may still export UCSC/raw tracks with converter_missing status.
The probe creates small task-owned temporary files and never ingests into a server.
"""
    tools = {name: shutil.which(name) for name in ("bedGraphToBigWig", "clodius")}
    require(all(tools.values()), "Browser environment requires bedGraphToBigWig and clodius on PATH")
    try:
        versions = {name: importlib.metadata.version(name) for name in ("clodius", "pyBigWig")}
        import pyBigWig  # noqa: F401 -- require the decoded round-trip validator
    except (ImportError, importlib.metadata.PackageNotFoundError) as error:
        raise ValueError("Browser environment requires installed clodius and pyBigWig packages") from error
    contracts = {
        "bedfile": ("--chromsizes-filename", "--output-file", "--delimiter"),
        "bedpe": ("--chromsizes-filename", "--output-file", "--chr1-col", "--from1-col", "--to1-col", "--chr2-col", "--from2-col", "--to2-col"),
    }
    help_evidence = {}
    for mode, required in contracts.items():
        command = [tools["clodius"], "aggregate", mode, "--help"]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        text = result.stdout + result.stderr
        require(all(option in text for option in required), f"Clodius {mode} CLI lacks required exporter options")
        help_evidence[mode] = {"command": command, "required_options": list(required), "help_sha256": hashlib.sha256(text.encode()).hexdigest()}
    with tempfile.TemporaryDirectory(prefix="microc2compare-browser-probe-") as directory:
        root = Path(directory)
        sizes = root / "profile.chromsizes.tsv"
        sizes.write_text("probe_chr2\t80\nprobe_chr1\t100\n", encoding="utf-8")
        raw = root / "signal.bedGraph"
        raw.write_text("probe_chr1\t0\t10\t0\nprobe_chr1\t20\t30\t-2.5\nprobe_chr2\t10\t20\t1.25\n", encoding="utf-8")
        binary = root / "signal.bw"
        command = [tools["bedGraphToBigWig"], str(raw), str(sizes), str(binary)]
        subprocess.run(command, check=True, capture_output=True, text=True)
        require(binary.is_file() and binary.stat().st_size > 0, "Browser smoke conversion produced no bigWig")
        validation = _binary_check(binary, "bigwig", raw, {"probe_chr2": 80, "probe_chr1": 100})
        require("roundtrip_checked" in validation, "Browser smoke bigWig round-trip was not verified")
        smoke = {"status": "passed", "validation": validation, "input_rows": 3, "zero_and_negative_values_checked": True,
                 "output_sha256": sha256(binary), "output_bytes": binary.stat().st_size,
                 "chromosome_order": ["probe_chr2", "probe_chr1"], "coordinate_system_ingested": False, "render_tested": False}
    return {"status": "passed", "python": platform.python_version(), "packages": versions,
            "executables": {name: {"path": value, "sha256": sha256(value)} for name, value in tools.items()},
            "clodius_cli": help_evidence, "bigwig_smoke": smoke, "higlass_ingested": False, "higlass_render_tested": False}
