"""Fresh CLI registration/export with new samples, counts and reference identity."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest


REPO = Path(__file__).resolve().parents[1]
BASH = (Path("C:/Program Files/Git/bin/bash.exe") if os.name == "nt"
        else Path(shutil.which("bash") or "/unavailable/bash"))


def write_json(path, data):
    path.write_text(json.dumps(data), encoding="utf-8")


def write_table(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.skipif(not BASH.is_file(), reason="Bash unavailable")
def test_fresh_registered_two_versus_three_cli_and_export_action():
    # Git Bash must be able to write runtime/shell temporary files in the workspace.
    parent = REPO / ".devtools" / "browser-pipeline-tests"
    parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=parent) as directory:
        root = Path(directory)
        assembly = "synthetic_hg38_order_v9"
        sizes = root / "chromsizes.tsv"
        sizes.write_text("chr1\t1000000\nchr2\t500000\n")
        views = root / "views.bed"
        views.write_text("chr1\t0\t1000000\tchr1\nchr2\t0\t500000\tchr2\n")
        registry = root / "reference.json"
        write_json(registry, {"reference_id": assembly, "chromsizes": str(sizes), "chromsizes_sha256": digest(sizes),
                              "views": str(views), "views_sha256": digest(views),
                              "chromosome_classes": {"autosomes": ["chr1", "chr2"], "sex": []}})
        samples, membership, inputs = [], [], []
        for index in range(5):
            sid = f"NewObservation{index + 1}"
            condition = "reference_group" if index < 2 else "new_treatment"
            names = ["chr1", "chr2"] if index < 2 else ["chr2", "chr1"]
            bins = pd.DataFrame([(chrom, start, start + 100000) for chrom in names
                                 for start in range(0, 1000000 if chrom == "chr1" else 500000, 100000)],
                                columns=["chrom", "start", "end"])
            one, two = np.triu_indices(len(bins))
            pixels = pd.DataFrame({"bin1_id": one, "bin2_id": two, "count": np.ones(len(one), dtype="int64") * (index + 1)})
            matrix = root / f"{sid}.cool"
            cooler.create_cooler(str(matrix), bins, pixels, ordered=True, assembly=assembly)
            with h5py.File(matrix, "r+") as handle:
                weight = handle["bins"].create_dataset("weight", data=np.ones(len(bins)))
                weight.attrs["converged"] = True
            provenance = root / f"{sid}.provenance.json"
            write_json(provenance, {"reference_id": assembly, "count_units": "unique_contacts", "dedup_scope": "joint_library",
                                    "library_ids": [sid + "_library"], "source_ids": [sid + "_source"],
                                    "filters": {"mapq": 30, "min_distance_bp": 0, "chromosome_policy": "registered_autosomes"},
                                    "balancing": {"100000": {"weight_column": "weight", "convention": "multiplicative", "converged": True,
                                                               "recipe": {"kind": "known_synthetic_test_weights"}}}, "synthetic_fixture": True})
            samples.append({"sample_id": sid, "project": "fresh_project", "condition": condition, "assay": "microc",
                            "biological_replicate": str(index + 1 if index < 2 else index - 1), "subject_id": sid + "_subject",
                            "reference_id": assembly, "matrix_path": str(matrix), "provenance_path": str(provenance), "batch": "new_batch"})
            membership.append({"sample_id": sid, "library_id": sid + "_library", "source_id": sid + "_source", "source_kind": "independent_library",
                               "dedup_group_id": sid + "_dedup", "dedup_scope": "joint_library", "source_checksum": digest(matrix)})
            inputs += [matrix, provenance]
        samples_path, membership_path = root / "samples.tsv", root / "membership.tsv"
        write_table(samples_path, samples, list(samples[0]))
        write_table(membership_path, membership, list(membership[0]))
        before = {str(path): digest(path) for path in inputs}
        schema = json.loads((REPO / "schemas/config.schema.json").read_text())
        config = {key: "false" for key, value in schema["properties"].items() if key.startswith("RUN_") and value["type"] == "boolean"}
        config.update(PROJECT="fresh_project", ANALYSIS_SET_ID="new_cohort", RUN_ID="fresh_run", REFERENCE_ID=assembly,
                      ENVIRONMENT_MODE="development", PYTHON_BOOTSTRAP=Path(sys.executable).as_posix(), SAMPLESHEET=samples_path.as_posix(),
                      TECHNICAL_MEMBERSHIP=membership_path.as_posix(), COMPARISONS="", INDEPENDENT_FEATURES_REGISTRY="", REFERENCE_REGISTRY=registry.as_posix(),
                      OUTDIR=(root / "results").as_posix(), SCRATCH_ROOT=(root / "scratch").as_posix(), PS_RESOLUTION="100000",
                      TOTAL_CPUS="2", MAX_MEMORY_MB="4000", THREADS_PREFLIGHT="1", MEMORY_MB_PREFLIGHT="1000",
                      RUN_BROWSER_EXPORT="true", EXPORT_UCSC="true", EXPORT_HIGLASS="false", THREADS_BROWSER_EXPORT="1", MEMORY_MB_BROWSER_EXPORT="1000")
        project = root / "project.conf"
        project.write_text("\n".join(f'{key}="{value}"' for key, value in config.items()) + "\n", encoding="utf-8")
        runtime = root / "runtime"
        runtime.mkdir()
        environment = dict(os.environ, TMP=runtime.as_posix(), TEMP=runtime.as_posix(), TMPDIR=runtime.as_posix())
        command = [str(BASH), "--noprofile", "--norc", (REPO / "scripts/microc2compare.sh").as_posix(), "-c", project.as_posix()]
        result = subprocess.run(command + ["--until-stage", "browser_export"], cwd=REPO, env=environment,
                                capture_output=True, text=True, timeout=180)
        assert result.returncode == 0, result.stdout + result.stderr
        run = root / "results/fresh_project/new_cohort/fresh_run"
        audits = list((run / "browser_tracks").glob("*/export_audit.json"))
        assert len(audits) == 1
        audit = json.loads(audits[0].read_text())
        assert audit["assembly"] == assembly
        assert {condition: len(members) for condition, members in audit["condition_members"].items()} == {"reference_group": 2, "new_treatment": 3}
        assert audit["status"] == "passed" and not audit["contact_matrices_opened"]
        assert (run / "completion/run/browser_export.done.json").is_file()
        # An explicit export action is independently usable on the fresh run.
        action_out = root / "second_browser_bundle"
        result = subprocess.run(command + ["--export-browser", "--out", action_out.as_posix()], cwd=REPO, env=environment,
                                capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
        assert json.loads((action_out / "export_audit.json").read_text())["status"] == "passed"
        # Zero-call results have valid empty headerless native tracks. Exercise
        # the real stage wrapper/completion guard, rather than only the exporter.
        empty_calls = run / "03_features/loops/NewObservation1.100000.bedpe.tsv"
        empty_calls.parent.mkdir(parents=True, exist_ok=True)
        write_table(empty_calls, [], ["feature_id", "chrom1", "start1", "end1",
                                     "chrom2", "start2", "end2", "resolution"])
        empty_out = root / "empty_native_browser_bundle"
        # Refresh normal admission when changing the output option. The real
        # CLI then reaches scientific.py with current preflight prerequisites.
        with project.open("a", encoding="utf-8") as stream:
            stream.write(f'BROWSER_OUTPUT_DIR="{empty_out.as_posix()}"\n')
        result = subprocess.run(command + ["--resume", "--until-stage", "browser_export"],
                                cwd=REPO, env=environment, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "SUCCESS browser_export/run" in (run / "logs/browser_export.run.log").read_text()
        with (empty_out / "inventory.tsv").open(newline="") as stream:
            empty_tracks = list(csv.DictReader(stream, delimiter="\t"))
        assert empty_tracks and all(row["rows"] == "0" for row in empty_tracks)
        for row in empty_tracks:
            for key in ("raw_path", "bedpe_path"):
                if row[key]:
                    assert (empty_out / row[key]).stat().st_size == 0
        completion = json.loads((run / "completion/run/browser_export.done.json").read_text())
        assert completion["status"] == "success"
        assert json.loads((empty_out / "export_audit.json").read_text())["status"] == "passed"
        assert before == {str(path): digest(path) for path in inputs}
