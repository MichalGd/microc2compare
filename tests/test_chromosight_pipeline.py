"""Independent Chromosight discovery supplies its own central-stage universe."""
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from microc2compare.analyses import configure_selection
from microc2compare.contracts import atomic_json, enabled, file_hash, normalize, read_json, write_table
from microc2compare.state import complete, initialize, scopes


REPO = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("measurement", ["quantify", "apa", "chromosight", "annotation"])
def test_explicit_chromosight_discovery_supplies_selected_measurements(config, measurement):
    config["INDEPENDENT_FEATURES_REGISTRY"] = ""
    selected = configure_selection(config, "chromosight_discovery," + measurement)
    assert selected["RUN_CHROMOSIGHT_DISCOVERY"]
    assert selected["RUN_QC"]
    assert not selected["RUN_LOOPS"]
    assert enabled(selected, "quantify")


def test_explicit_both_callers_are_retained(config):
    selected = configure_selection(config, "loops,chromosight_discovery,quantify")
    assert selected["RUN_LOOPS"] and selected["RUN_CHROMOSIGHT_DISCOVERY"]


def test_fresh_chromosight_only_central_stages(tmp_path, config):
    cooler = pytest.importorskip("cooler")
    h5py = pytest.importorskip("h5py")
    pytest.importorskip("chromosight.utils.detection")
    from importlib.resources import files

    assembly, chrom, resolution = "independent_assembly_v7", "ctgNovel", 2000
    length, n = 320000, 160
    sizes, views, registry = tmp_path / "sizes.tsv", tmp_path / "views.bed", tmp_path / "reference.json"
    sizes.write_text(f"{chrom}\t{length}\n", encoding="utf-8")
    views.write_text(f"{chrom}\t0\t{length}\t{chrom}\n", encoding="utf-8")
    atomic_json(registry, {"reference_id": assembly, "chromsizes": str(sizes), "chromsizes_sha256": file_hash(sizes),
                          "views": str(views), "views_sha256": file_hash(views),
                          "chromosome_classes": {"autosomes": [chrom], "sex": []}})
    bins = cooler.binnify(pd.Series({chrom: length}), resolution)
    kernel = np.loadtxt(str(files("chromosight.kernels").joinpath("artificial_template_loops_type1.txt")))
    radius = kernel.shape[0] // 2
    x, y = np.triu_indices(n)
    base = np.full(len(x), 10000.)
    peak = np.zeros((n, n))
    # The published artificial template includes a baseline of one. Inject
    # only its excess to avoid manufacturing a rectangular edge pattern.
    peak[55-radius:56+radius, 95-radius:96+radius] = (kernel - 1) * 10000
    samples, membership, originals = [], [], [sizes, views, registry]
    for index, sid in enumerate(("nova-A", "nova.B", "nova_C")):
        matrix, provenance = tmp_path / (sid + ".cool"), tmp_path / (sid + ".json")
        counts = np.rint((base + peak[x, y]) * (1 + index * .07)).astype(np.int64)
        cooler.create_cooler(str(matrix), bins, pd.DataFrame({"bin1_id": x, "bin2_id": y, "count": counts}),
                             ordered=True, assembly=assembly)
        with h5py.File(matrix, "r+") as handle:
            weight = handle["bins"].create_dataset("weight", data=np.ones(n))
            weight.attrs["converged"] = True
        atomic_json(provenance, {"reference_id": assembly, "count_units": "unique_contacts", "dedup_scope": "joint_library",
                                "library_ids": [sid + "_library"], "source_ids": [sid + "_source"],
                                "filters": {"mapq": 30, "min_distance_bp": 0, "chromosome_policy": "registered_autosomes"},
                                "balancing": {str(resolution): {"weight_column": "weight", "convention": "multiplicative", "converged": True,
                                                                  "recipe": {"kind": "known_synthetic_test_weights"}}},
                                "synthetic_fixture": True})
        samples.append({"sample_id": sid, "project": "unrelated_project", "condition": "baseline" if index < 2 else "exposed",
                        "assay": "microc", "biological_replicate": str(index + 1 if index < 2 else 1),
                        "subject_id": sid + "_donor", "reference_id": assembly, "matrix_path": str(matrix), "provenance_path": str(provenance)})
        membership.append({"sample_id": sid, "library_id": sid + "_library", "source_id": sid + "_source", "source_kind": "independent_library",
                           "dedup_group_id": sid + "_dedup", "dedup_scope": "joint_library", "source_checksum": file_hash(matrix)})
        originals += [matrix, provenance]
    samples_path, membership_path = tmp_path / "samples.tsv", tmp_path / "membership.tsv"
    write_table(samples_path, samples, list(samples[0]))
    write_table(membership_path, membership, list(membership[0]))
    originals += [samples_path, membership_path]
    before = {str(p): file_hash(p) for p in originals}
    values = {key: str(value).lower() if isinstance(value, bool) else str(value) for key, value in config.items()}
    values.update(PROJECT="unrelated_project", ANALYSIS_SET_ID="independent_discovery", REFERENCE_ID=assembly,
                  SAMPLESHEET=str(samples_path), TECHNICAL_MEMBERSHIP=str(membership_path), REFERENCE_REGISTRY=str(registry),
                  COMPARISONS="", INDEPENDENT_FEATURES_REGISTRY="", ANNOTATION_REGISTRY="", ENVIRONMENT_MODE="development",
                  OUTDIR=str(tmp_path / "results"), SCRATCH_ROOT=str(tmp_path / "scratch"),
                  PS_RESOLUTION="2000", LOOP_RESOLUTION="2000", LOOP_DISCOVERY_RESOLUTIONS="2000", LOOP_QUANTIFICATION_RESOLUTIONS="2000",
                  LOOP_CALLER="mustache", LOOP_MIN_DISTANCE_BP="6000", LOOP_MAX_DISTANCE_BP="200000", MAX_DENSE_BINS="512",
                  CHROMOSIGHT_DISCOVERY_RESOLUTIONS="2000", CHROMOSIGHT_DISCOVERY_MIN_PEARSON="0.3",
                  TOTAL_CPUS="2", MAX_MEMORY_MB="4000", THREADS_CHROMOSIGHT="1", MEMORY_MB_CHROMOSIGHT="2000",
                  RUN_FINAL_REPORT="true", M2C_ANALYSES="chromosight_discovery,quantify")
    for kind in ("PREFLIGHT", "QC", "QUANTIFY", "REPORT"):
        values.update({f"THREADS_{kind}": "1", f"MEMORY_MB_{kind}": "2000"})
    c = normalize(values, tmp_path)
    assert c["RUN_CHROMOSIGHT_DISCOVERY"] and not c["RUN_LOOPS"]
    initialize(c, os.getpid(), "fresh-chromosight-test")
    root = Path(c["run_dir"])
    resolved = root / "00_metadata/config.json"
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1", TMP=str(tmp_path), TEMP=str(tmp_path))
    for stage in ("preflight", "register", "qc", "chromosight_discovery", "consensus", "quantify", "report"):
        assert enabled(c, stage)
        for sid in scopes(c, stage):
            result = subprocess.run([sys.executable, "-I", "-B", str(REPO / "scripts/scientific.py"), stage,
                                     "--config", str(resolved), "--scope", sid], cwd=REPO, env=environment,
                                    capture_output=True, text=True, timeout=600 if stage == "report" else 90)
            assert result.returncode == 0, stage + "/" + sid + "\n" + result.stdout + result.stderr
            assert complete(c, stage, sid), stage + "/" + sid
    assert not (root / "completion/nova-A/loops.done.json").exists()
    assert not list((root / "03_features/loops").glob("*.bedpe.tsv"))
    assert read_json(root / "00_metadata/environments.json")["loops"]["adapter_probe"]["versions"] == {"chromosight": "1.6.3"}
    for sample in samples:
        calls = pd.read_csv(root / f"03_features/loops/{sample['sample_id']}.2000.chromosight.tsv", sep="\t")
        assert len(calls) and calls.caller.eq("chromosight").all()
        closest = calls.iloc[np.argmin(abs(calls.start1 - 55*resolution) + abs(calls.start2 - 95*resolution))]
        assert abs(closest.start1 - 55*resolution) <= resolution and abs(closest.start2 - 95*resolution) <= resolution
    catalogue = pd.read_csv(root / "03_features/loops/catalogue.tsv", sep="\t")
    assert len(catalogue) and catalogue.supporting_callers.eq("chromosight").all()
    assert catalogue.biological_support.max() == 3
    counts = pd.read_csv(root / "04_quantification/discovered_loops/2000/counts.tsv", sep="\t")
    assert len(counts) and counts[[s["sample_id"] for s in samples]].notna().all().all()
    report = (root / "08_reports/index.html").read_text(encoding="utf-8")
    assert "Multi-resolution loop catalogue" in report and "nova-A.2000.chromosight.tsv" in report
    assert "Mustache 1.3.3 suppresses" not in report
    assert before == {str(p): file_hash(p) for p in originals}
    assert not list(root.rglob("*.cool")) and not list(root.rglob("*.mcool"))
