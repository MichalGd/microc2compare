"""Generic export acceptance: new cohorts, references and optional modules."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
from unittest.mock import patch
from types import SimpleNamespace

import pytest

from microc2compare.browser_export import export_browser, _binary_check


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def write_table(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def table(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def bundle(tmp_path, assembly="hg38", signals=True):
    run = tmp_path / "new-study" / "results" / "new-run"
    reference = tmp_path / "reference"
    reference.mkdir()
    sizes = reference / "sizes.tsv"
    sizes.write_text("chr1\t248956422\nchr2\t242193529\n", encoding="utf-8")
    registry = reference / "registry.json"
    write_json(registry, {"reference_id": assembly, "chromsizes": "sizes.tsv", "chromsizes_sha256": hashlib.sha256(sizes.read_bytes()).hexdigest()})
    samples = [{"sample_id": sid, "condition": condition, "reference_id": assembly}
               for sid, condition in (("Alpha1", "control"), ("Alpha2", "control"), ("Beta1", "stimulated"), ("Beta2", "stimulated"), ("Beta3", "stimulated"))]
    config = {"REFERENCE_ID": assembly, "REFERENCE_REGISTRY": str(registry), "run_dir": str(run), "RUN_ID": "new-run",
              "COMPARTMENT_RESOLUTION": 100, "INSULATION_RESOLUTION": 10, "INSULATION_WINDOWS_BP": "20,40",
              "INSULATION_PRIMARY_WINDOW_BP": 20, "ADA_RESOLUTION": 10, "MIN_PHASING_CORRELATION": .2}
    write_json(run / "00_metadata/config.json", config)
    write_json(run / "00_metadata/samples.json", {"samples": samples, "comparisons": []})
    axes = [{"chrom": "chr1", "start": start, "end": start + 10} for start in (0, 10, 20)]
    write_table(run / "00_metadata/masks/10.bins.tsv", axes, ("chrom", "start", "end"))
    write_table(run / "00_metadata/masks/100.bins.tsv", axes, ("chrom", "start", "end"))
    write_table(run / "00_metadata/masks/100.views.tsv", [{"chrom": "chr1", "start": 0, "end": 30, "name": "chr1"}], ("chrom", "start", "end", "name"))
    write_json(run / "00_metadata/masks/10.scope.json", {"resolution": 10, "analysis_chromosomes": ["chr1"], "samples": [
        {"sample_id": sample["sample_id"], "native_chromosomes": ["chr1", "chr2"] if i < 2 else ["chr2", "chr1"]} for i, sample in enumerate(samples)]})
    if signals:
        for index, sample in enumerate(samples):
            sid = sample["sample_id"]
            compartments = [dict(axis, E1=index + 1, E2=2, E3=-2, oriented_score=index + 1, selected_component="E1") for axis in axes]
            write_table(run / f"03_features/compartments/{sid}.tsv", compartments, (*axes[0], "E1", "E2", "E3", "oriented_score", "selected_component"))
            write_json(run / f"03_features/compartments/{sid}.phasing.json", [{"region": "chr1", "status": "oriented", "component": "E1", "correlation": .8}])
            insulation = [dict(axis, log2_insulation_score_20=-index, boundary_strength_20=.5 if j == 1 else "NA", is_boundary_20=j == 1,
                               log2_insulation_score_40=-index / 2, boundary_strength_40=.6 if j == 1 else "NA", is_boundary_40=j == 1) for j, axis in enumerate(axes)]
            write_table(run / f"03_features/insulation/{sid}.tsv", insulation, tuple(insulation[0]))
    return run, config, samples


def exported(tmp_path, config):
    out = tmp_path / "browser-output"
    outputs = export_browser(dict(config, BROWSER_OUTPUT_DIR=str(out)))
    assert out / "artifact_checksums.tsv" in outputs
    return out, table(out / "inventory.tsv")


def test_new_two_versus_three_cohort_and_hg38_without_loops(tmp_path):
    run, config, samples = bundle(tmp_path)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in run.rglob("*") if p.is_file()}
    out, inventory = exported(tmp_path, config)
    assert inventory and {row["assembly"] for row in inventory} == {"hg38"}
    assert not any(row["analysis"] == "loops" for row in inventory)
    means = {row["label"]: row for row in inventory if "selected_oriented.mean_complete" in row["path"]}
    assert set(means) == {"control", "stimulated"}
    for condition, value, denominator in (("control", 1.5, 2), ("stimulated", 4, 3)):
        text = (out / means[condition]["path"]).read_text(encoding="utf-8")
        assert f"all {denominator} observations required" in text
        assert float(text.splitlines()[1].split("\t")[3]) == value
    assert any("window_40" in row["path"] for row in inventory)
    assert all(before[p] == hashlib.sha256(p.read_bytes()).hexdigest() for p in before)
    audit = json.loads((out / "export_audit.json").read_text())
    assert audit["condition_members"]["stimulated"] == [s["sample_id"] for s in samples if s["condition"] == "stimulated"]
    assert not audit["contact_matrices_opened"]
    assert "basal" not in (out / "README.md").read_text()
    assert "mm39" not in (out / "README.md").read_text()


def add_loops(run, samples):
    native = {"feature_id": "native_1", "chrom1": "chr1", "start1": 10, "end1": 20, "chrom2": "chr1", "start2": 110, "end2": 120, "resolution": 10}
    write_table(run / f"03_features/loops/{samples[0]['sample_id']}.10.bedpe.tsv", [native], tuple(native))
    catalogue = dict(native, feature_id="candidate_X", start1=20, end1=30, representative_resolution=10, supporting_resolutions="10;20")
    write_table(run / "03_features/loops/catalogue.tsv", [catalogue], tuple(catalogue))
    membership = [dict(native, source=samples[i]["sample_id"], candidate_id="candidate_X") for i in (0, 0, 2, 3)]
    membership[1]["resolution"] = 20
    write_table(run / "03_features/loops/source_membership.tsv", membership, tuple(membership[0]))


def test_native_coordinates_and_condition_support_without_annotation(tmp_path):
    run, config, samples = bundle(tmp_path, assembly="custom_v7", signals=False)
    add_loops(run, samples)
    out, inventory = exported(tmp_path, config)
    native = next(row for row in inventory if row["analysis"] == "loops" and "native_10" in row["path"])
    representative = next(row for row in inventory if row["analysis"] == "loops" and "all_loops" in row["path"])
    assert (out / native["bedpe_path"]).read_text().split("\t")[1] == "10"
    assert (out / representative["bedpe_path"]).read_text().split("\t")[1] == "20"
    summary = table(out / "loops/shared/detection_support.tsv")[0]
    assert summary["support:control"] == "1" and summary["support:stimulated"] == "2"
    control = next(row for row in inventory if row["analysis"] == "loops" and row["label"] == "control" and "support_ge2" in row["path"])
    assert control["rows"] == "0"
    assert "unannotated" in (out / representative["path"]).read_text()


def test_tads_and_exact_ada_keep_saved_ids_and_arbitrary_denominators(tmp_path):
    run, config, samples = bundle(tmp_path, signals=False)
    conditions = {condition: [s["sample_id"] for s in samples if s["condition"] == condition] for condition in ("control", "stimulated")}
    fields = ("feature_id", "chrom", "start", "end", "window_bp", "mean_log2_insulation", "left_boundary_strength", "right_boundary_strength", "support", "n_samples", "support_fraction", "sample_ids")
    write_json(run / "03_features/tads/parameters.json", {"conditions": {k: k for k in conditions}})
    for window in (20, 40):
        for sample in samples:
            row = dict(zip(fields, ("saved_tad", "chr1", 0, 30, window, -.5, .8, .7, 1, 1, 1, sample["sample_id"])))
            write_table(run / f"03_features/tads/samples/{sample['sample_id']}.{window}.tsv", [row], fields)
            write_table(run / f"03_features/tads/samples/{sample['sample_id']}.{window}.exclusions.tsv", [], ("chrom", "start", "end", "window_bp", "reason"))
        for condition, members in conditions.items():
            row = dict(zip(fields, ("saved_condition_tad", "chr1", 0, 30, window, -.5, .8, .7, 1, len(members), 1 / len(members), members[0])))
            write_table(run / f"03_features/tads/conditions/{condition}.{window}.tsv", [row], fields)
            write_table(run / f"03_features/tads/conditions/{condition}.{window}.unsupported.tsv", [], fields)
    domain = {"feature_id": "saved_ADA_exact_id", "chrom": "chr1", "start": 0, "end": 30, "window_bp": 20}
    write_table(run / "03_features/boundary_derived_domains/consensus.tsv", [domain], tuple(domain))
    strengths = [{"feature_id": domain["feature_id"], "sample_id": s["sample_id"], "mean_intra_domain_oe": 0, "oe_sum": 0, "eligible_pixel_count": 2, "inference_status": "descriptive"} for s in samples]
    write_table(run / "04_quantification/domain_scores/strength.tsv", strengths, tuple(strengths[0]))
    out, inventory = exported(tmp_path, config)
    ada = next(row for row in inventory if row["analysis"] == "ada" and row["scope"] == "condition" and row["label"] == "stimulated" and "complete_case_mean_oe" in row["path"])
    assert (out / ada["raw_path"]).read_text().rstrip().endswith("\t0")
    coverage = next(row for row in inventory if row["analysis"] == "ada" and row["scope"] == "condition" and row["label"] == "stimulated" and "n_finite" in row["path"])
    assert (out / coverage["raw_path"]).read_text().rstrip().endswith("\t3")
    assert "saved_ADA_exact_id" in (out / next(r for r in inventory if r["analysis"] == "ada" and r["format"] == "bed")["raw_path"]).read_text()


def test_native_axis_profiles_and_missing_converters_are_explicit(tmp_path):
    run, config, samples = bundle(tmp_path)
    with patch("microc2compare.browser_export.shutil.which", return_value=None):
        out, inventory = exported(tmp_path, dict(config, EXPORT_HIGLASS=True))
    profiles = json.loads((out / "reference/coordinate_profiles.json").read_text())["profiles"]
    assert len(profiles) == 3
    assert {tuple(p["chromosomes"][i][0] for i in range(len(p["chromosomes"]))) for p in profiles} == {("chr1",), ("chr1", "chr2"), ("chr2", "chr1")}
    assert all(not Path(p["path"]).is_absolute() and (out / p["path"]).is_file() for p in profiles)
    assert {r["higlass_status"] for r in inventory} <= {"converter_missing: bedGraphToBigWig", "converter_missing: clodius", "empty; conversion not required"}
    assert not json.loads((out / "export_audit.json").read_text())["higlass_ingested"]


def test_converter_receipts_validate_binary_and_never_claim_ingestion(tmp_path):
    run, config, samples = bundle(tmp_path)
    calls = []
    def convert(command, **kwargs):
        calls.append(command)
        if command[0] == "fake-bigwig":
            output = Path(command[-1])
            output.write_bytes(b"\x26\xfc\x8f\x88" + b"fixture")
        else:
            import sqlite3
            output = Path(command[command.index("--output-file") + 1])
            with sqlite3.connect(output) as database:
                database.execute("CREATE TABLE intervals (start INTEGER, end INTEGER)")
        # UCSC converters can succeed silently; the stage still requires an
        # explicit nonempty conversion log carrying the command and status.
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
    def utility(name):
        return "fake-bigwig" if name == "bedGraphToBigWig" else "fake-clodius"
    with patch("microc2compare.browser_export.shutil.which", side_effect=utility), patch("microc2compare.browser_export.subprocess.run", side_effect=convert), patch.dict("sys.modules", {"pyBigWig": None}):
        out, inventory = exported(tmp_path, dict(config, EXPORT_HIGLASS=True))
    converted = [r for r in inventory if int(r["rows"]) > 0]
    assert all(r["higlass_status"] == "converted; not ingested or render-tested" for r in converted)
    receipts = json.loads((out / "higlass/conversion_receipts.json").read_text())
    assert receipts and all(not r["ingested"] and not r["render_tested"] for r in receipts)
    assert all((out / (r["path"] + ".conversion.log")).read_text().startswith("Exit code: 0\nCommand argv: ") for r in receipts)
    assert {r["filetype"] for r in receipts} == {"bigwig", "beddb"}
    recipes = json.loads((out / "higlass/ingestion_commands.json").read_text())["commands"]
    assert recipes and all(command[0] == "higlass-manage" for command in recipes)
    assert all(command[0] != "higlass-manage" for command in calls)


def test_invalid_converter_output_has_failed_export_audit(tmp_path):
    run, config, samples = bundle(tmp_path)
    def convert(command, **kwargs):
        Path(command[-1]).write_bytes(b"not-a-bigwig")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
    out = tmp_path / "failed-binary"
    with patch("microc2compare.browser_export.shutil.which", return_value="fake-bigwig"), patch("microc2compare.browser_export.subprocess.run", side_effect=convert):
        with pytest.raises(ValueError, match="not bigWig"):
            export_browser(dict(config, EXPORT_HIGLASS=True, BROWSER_OUTPUT_DIR=str(out)))
    assert not (out / "export_audit.json").is_file()
    assert json.loads((out / "export_failure.json").read_text())["status"] == "failed"


def test_optional_stripe_seed_overlap_preserved(tmp_path):
    run, config, samples = bundle(tmp_path, signals=False)
    row = {"stripe_id": "stripe_native", "sample_id": "Alpha1", "condition": "control", "reference_id": "hg38", "resolution": 10,
           "chrom1": "chr1", "start1": 10, "end1": 20, "chrom2": "chr1", "start2": 10, "end2": 80,
           "direction": "downstream", "matrix_location": "LT", "seed_bp": 10}
    write_table(run / "03_features/stripes/Alpha1.10.tsv", [row], tuple(row))
    summary = {"sample_id": "Alpha1", "condition": "control", "reference_id": "hg38", "resolution": 10,
               "stripe_count": 1, "biological_replicates_in_condition": 2}
    write_table(run / "03_features/stripes/summary.tsv", [summary], tuple(summary))
    write_table(run / "03_features/stripes/matches.tsv", [], ("resolution", "sample1", "sample2", "stripe1", "stripe2", "direction"))
    write_json(run / "03_features/stripes/summary.json", {"condition_members": {condition: [s["sample_id"] for s in samples if s["condition"] == condition] for condition in ("control", "stimulated")}})
    out, inventory = exported(tmp_path, config)
    geometry = next(r for r in inventory if r["format"] == "bedpe")
    assert (out / geometry["bedpe_path"]).read_text().split("\t")[:6] == ["chr1", "10", "20", "chr1", "10", "80"]
    assert not any(r["format"] == "interact" for r in inventory)
    assert len(list((out / "source_tables/stripes/support").glob("*"))) == 3


def add_condition_stripes(run, samples, with_support=True):
    calls = {}
    for sid, values in {
        "Alpha1": [("shared_id", 10, 20, 80, "12.3400")],
        "Alpha2": [],
        "Beta1": [("shared_id", 10, 20, 80, "23.4500")],
        "Beta2": [("stripe_b1", 10, 30, 90, "34.5600"), ("stripe_b2", 10, 40, 100, "45.6700")],
        "Beta3": [("stripe_c", 20, 30, 100, "56.7800")],
    }.items():
        condition = next(s["condition"] for s in samples if s["sample_id"] == sid)
        rows = [dict(stripe_id=fid, sample_id=sid, condition=condition, reference_id="hg38", resolution=10,
                     chrom1="chr1", start1=start, end1=end, chrom2="chr1", start2=start, end2=body_end,
                     direction="downstream", matrix_location="LT", seed_bp=start, relative_change=change,
                     original_descriptor="retained without averaging") for fid, start, end, body_end, change in values]
        fields = tuple(rows[0]) if rows else ("stripe_id", "sample_id", "condition", "reference_id", "resolution",
                  "chrom1", "start1", "end1", "chrom2", "start2", "end2", "direction", "matrix_location", "seed_bp", "relative_change", "original_descriptor")
        write_table(run / f"03_features/stripes/{sid}.10.tsv", rows, fields)
        write_table(run / f"03_features/stripes/{sid}.20.tsv", [], fields)
        calls[sid] = rows
    if not with_support:
        return calls
    summaries = [dict(sample_id=s["sample_id"], condition=s["condition"], reference_id="hg38", resolution=res,
                      stripe_count=len(calls[s["sample_id"]]) if res == 10 else 0,
                      biological_replicates_in_condition=sum(x["condition"] == s["condition"] for x in samples))
                 for res in (10, 20) for s in samples]
    write_table(run / "03_features/stripes/summary.tsv", summaries, tuple(summaries[0]))
    matches = [dict(resolution=10, sample1=a, stripe1=fa, sample2=b, stripe2=fb,
                    same_condition=(a.startswith("Beta") == b.startswith("Beta")), direction="downstream")
               for a, fa, b, fb in [
                   ("Alpha1", "shared_id", "Beta1", "shared_id"),
                   ("Alpha1", "shared_id", "Beta2", "stripe_b1"),
                   ("Alpha1", "shared_id", "Beta2", "stripe_b2"),
                   ("Beta1", "shared_id", "Beta2", "stripe_b1"),
                   ("Beta1", "shared_id", "Beta2", "stripe_b2"),
                   ("Beta2", "stripe_b1", "Beta3", "stripe_c"),
                   ("Beta2", "stripe_b2", "Beta3", "stripe_c"),
                   # Duplicate and reversed links still count the peer once.
                   ("Beta1", "shared_id", "Beta2", "stripe_b1"),
                   ("Beta2", "stripe_b1", "Beta1", "shared_id")]]
    write_table(run / "03_features/stripes/matches.tsv", matches, tuple(matches[0]))
    write_json(run / "03_features/stripes/summary.json",
               {"condition_members": {condition: [s["sample_id"] for s in samples if s["condition"] == condition]
                                      for condition in ("control", "stimulated")}})
    return calls


def test_condition_stripes_preserve_native_calls_and_direct_distinct_observation_support(tmp_path):
    run, config, samples = bundle(tmp_path, signals=False)
    calls = add_condition_stripes(run, samples)
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in run.rglob("*") if p.is_file()}
    with patch("microc2compare.browser_export.shutil.which", return_value=None):
        out, inventory = exported(tmp_path, dict(config, EXPORT_HIGLASS=True))
    rows = table(out / "stripes/conditions/stimulated/direct_support_10.tsv")
    by_id = {(row["sample_id"], row["stripe_id"]): row for row in rows}
    assert len(rows) == 4  # Two calls from Beta2 retain two different geometries.
    a = by_id[("Beta1", "shared_id")]
    assert a["n_direct_support_observations"] == "2" and a["n_direct_support_other_observations"] == "1"
    assert a["n_condition_observations"] == "3" and float(a["direct_support_fraction_all"]) == pytest.approx(2/3)
    assert a["direct_support_sample_ids"] == "Beta1;Beta2"  # No cross-condition or transitive Beta3 support.
    assert json.loads(a["direct_matched_call_ids"]) == [["Beta2", "stripe_b1"], ["Beta2", "stripe_b2"]]
    assert by_id[("Beta2", "stripe_b1")]["n_direct_support_observations"] == "3"
    assert by_id[("Beta3", "stripe_c")]["n_direct_support_observations"] == "2"
    control = table(out / "stripes/conditions/control/direct_support_10.tsv")[0]
    assert control["n_direct_support_observations"] == "1" and control["n_condition_observations"] == "2"
    assert control["n_condition_observations_with_saved_tables"] == "2"  # Empty Alpha2 table is retained.
    assert control["browser_id"] != a["browser_id"]  # Original ID shared across observations remains distinct.
    geometry = next(r for r in inventory if r["analysis"] == "stripes" and r["scope"] == "condition"
                    and r["label"] == "stimulated" and r["resolution"] == "10" and r["format"] == "bedpe")
    assert geometry["higlass_status"].startswith("converter_missing") and geometry["higlass_filetype"] == "bed2ddb"
    exported_geometry = [line.split("\t") for line in (out / geometry["bedpe_path"]).read_text().splitlines()]
    for native in calls["Beta1"] + calls["Beta2"] + calls["Beta3"]:
        saved = by_id[(native["sample_id"], native["stripe_id"])]
        assert all(saved[field] == str(value) for field, value in native.items())
        expected = [str(native[field]) for field in ("chrom1", "start1", "end1", "chrom2", "start2", "end2")]
        assert next(row[:6] for row in exported_geometry if row[6] == saved["browser_id"]) == expected
    assert all(before[p] == hashlib.sha256(p.read_bytes()).hexdigest() for p in before)
    empty = [r for r in inventory if r["analysis"] == "stripes" and r["scope"] == "condition" and r["resolution"] == "20"]
    assert len(empty) == 6 and all(r["rows"] == "0" for r in empty)
    assert table(out / "stripes/conditions/stimulated/direct_support_20.tsv") == []
    assert json.loads((out / "stripes/conditions/stimulated/direct_support_20.json").read_text())["native_call_count"] == 0


def test_condition_stripes_without_support_products_remain_unverified(tmp_path):
    run, config, samples = bundle(tmp_path, signals=False)
    add_condition_stripes(run, samples, with_support=False)
    out, inventory = exported(tmp_path, config)
    rows = table(out / "stripes/conditions/stimulated/direct_support_10.tsv")
    assert rows and all(r["direct_pairwise_support_status"] == "unverified_support_unavailable" for r in rows)
    assert all(r["n_direct_support_observations"] == r["direct_support_fraction_all"] == "NA" for r in rows)
    assert all(r["n_condition_observations"] == "3" for r in rows)
    assert any(r["analysis"] == "stripes" and r["scope"] == "condition" and r["format"] == "bedpe" for r in inventory)
    assert any(r["feature"] == "Stripe condition direct support" for r in table(out / "missing_features.tsv"))


def test_condition_stripes_validated_empty_matches_keep_source_only_support(tmp_path):
    run, config, samples = bundle(tmp_path, signals=False)
    add_condition_stripes(run, samples)
    # Distinct observations have disjoint anchors, so no pairwise match is valid.
    for index, sample in enumerate(samples):
        path = run / f"03_features/stripes/{sample['sample_id']}.10.tsv"
        rows = table(path)
        for row in rows:
            for field in ("start1", "end1", "start2", "end2", "seed_bp"):
                row[field] = str(int(row[field]) + index * 1000)
        with path.open(encoding="utf-8") as stream:
            fields = next(csv.reader(stream, delimiter="\t"))
        write_table(path, rows, fields)
    path = run / "03_features/stripes/matches.tsv"
    with path.open(encoding="utf-8") as stream:
        fields = next(csv.reader(stream, delimiter="\t"))
    write_table(path, [], fields)
    out, _ = exported(tmp_path, config)
    for condition, denominator in (("control", "2"), ("stimulated", "3")):
        rows = table(out / f"stripes/conditions/{condition}/direct_support_10.tsv")
        assert rows and all(r["n_direct_support_observations"] == "1" for r in rows)
        assert all(r["n_direct_support_other_observations"] == "0" for r in rows)
        assert all(r["n_condition_observations"] == denominator for r in rows)
        assert all(r["direct_pairwise_support_status"] == "saved_direct_pairwise_links" for r in rows)
        assert all(json.loads(r["direct_matched_call_ids"]) == [] for r in rows)


@pytest.mark.parametrize("corruption", ["direction", "same_condition", "self", "anchor_overlap", "chromosome"])
def test_condition_stripes_reject_inconsistent_direct_support_links(tmp_path, corruption):
    run, config, samples = bundle(tmp_path, signals=False)
    add_condition_stripes(run, samples)
    path = run / "03_features/stripes/matches.tsv"
    matches = table(path)
    if corruption == "direction":
        matches[0]["direction"] = "upstream"
    elif corruption == "same_condition":
        matches[0]["same_condition"] = "True"
    elif corruption == "self":
        matches[0].update(sample2=matches[0]["sample1"], stripe2=matches[0]["stripe1"])
    else:
        matches[0].update(sample2="Beta3", stripe2="stripe_c")
        if corruption == "chromosome":
            call_path = run / "03_features/stripes/Beta3.10.tsv"
            calls = table(call_path)
            calls[0].update(chrom1="chr2", chrom2="chr2")
            write_table(call_path, calls, tuple(calls[0]))
            for resolution in (10, 100):
                grid_path = run / f"00_metadata/masks/{resolution}.bins.tsv"
                grid = table(grid_path) + [{"chrom": "chr2", "start": "20", "end": "30"}]
                write_table(grid_path, grid, ("chrom", "start", "end"))
    write_table(path, matches, tuple(matches[0]))
    with pytest.raises(ValueError, match="Stripe support link"):
        exported(tmp_path, config)
    assert not (tmp_path / "browser-output/export_audit.json").exists()


def test_bigwig_decoded_roundtrip_checks_float32_values_and_coordinates(tmp_path):
    raw = tmp_path / "signal.bedGraph"
    raw.write_text("chr1\t0\t10\t0.1\nchr1\t10\t20\t0.1\nchr1\t20\t30\t0\n")
    binary = tmp_path / "signal.bw"
    binary.write_bytes(b"\x26\xfc\x8f\x88" + b"fixture")
    values = [(0, 20, 0.10000000149011612), (20, 30, 0.0)]
    reader = SimpleNamespace(isBigWig=lambda: True, chroms=lambda: {"chr1": 30},
                             intervals=lambda chrom: values, close=lambda: None)
    with patch.dict("sys.modules", {"pyBigWig": SimpleNamespace(open=lambda path: reader)}):
        assert "roundtrip_checked" in _binary_check(binary, "bigwig", raw, {"chr1": 30})
        values[1] = (20, 30, 1.0)
        with pytest.raises(ValueError, match="differ from source"):
            _binary_check(binary, "bigwig", raw, {"chr1": 30})


@pytest.mark.parametrize("corruption", ["reference", "sample_reference", "grid", "support", "partial_signals"])
def test_corrupt_inputs_fail_without_passed_audit(tmp_path, corruption):
    run, config, samples = bundle(tmp_path)
    if corruption == "reference":
        (tmp_path / "reference/sizes.tsv").write_text("chr1\t1\n")
    elif corruption == "sample_reference":
        samples[0]["reference_id"] = "wrong"
        write_json(run / "00_metadata/samples.json", {"samples": samples})
    elif corruption == "grid":
        write_table(run / "00_metadata/masks/10.bins.tsv", [{"chrom": "chr1", "start": 1, "end": 0}], ("chrom", "start", "end"))
    elif corruption == "support":
        add_loops(run, samples)
        path = run / "03_features/loops/source_membership.tsv"
        rows = table(path)
        rows[0]["source"] = "unknown"
        write_table(path, rows, tuple(rows[0]))
    else:
        (run / "03_features/compartments/Alpha1.tsv").unlink()
    out = tmp_path / "browser-output"
    with pytest.raises(ValueError):
        export_browser(dict(config, BROWSER_OUTPUT_DIR=str(out)))
    assert not (out / "export_audit.json").exists()


def test_default_destination_and_refuse_overwrite(tmp_path):
    run, config, samples = bundle(tmp_path, signals=False)
    outputs = export_browser(config)
    audit = next(p for p in outputs if p.name == "export_audit.json")
    assert audit.parent.parent == run / "browser_tracks"
    with pytest.raises(ValueError, match="must be new"):
        export_browser(dict(config, BROWSER_OUTPUT_DIR=str(audit.parent)))


def test_standalone_cli_uses_new_run_metadata(tmp_path):
    run, config, samples = bundle(tmp_path, signals=False)
    script = Path(__file__).resolve().parents[1] / "scripts/export_browser.py"
    import sys
    result = subprocess.run([sys.executable, "-I", "-B", str(script), "--run", str(run), "--out", str(tmp_path / "cli-bundle")], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["status"] == "passed"


@pytest.mark.parametrize("original_available", [False, True], ids=["copied-run", "original-available"])
@pytest.mark.parametrize("changed", [False, True], ids=["same-reference", "reference-drift"])
def test_relocated_reference_binds_to_original_frozen_evidence(tmp_path, original_available, changed):
    run, config, samples = bundle(tmp_path, signals=False)
    original_registry = Path(config["REFERENCE_REGISTRY"])
    original_sizes = original_registry.parent / "sizes.tsv"
    write_json(run / "00_metadata/inputs.json", {str(path): {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in (original_registry, original_sizes)})
    relocated = tmp_path / "relocated"
    relocated.mkdir()
    shutil.copyfile(original_sizes, relocated / "sizes.tsv")
    registry = json.loads(original_registry.read_text())
    if changed:
        (relocated / "sizes.tsv").write_text("chr1\t248956423\nchr2\t242193529\n")
        registry["chromsizes_sha256"] = hashlib.sha256((relocated / "sizes.tsv").read_bytes()).hexdigest()
    write_json(relocated / "registry.json", registry)
    if not original_available:
        original_registry.unlink()
    active = dict(config, REFERENCE_REGISTRY=str(relocated / "registry.json"), BROWSER_OUTPUT_DIR=str(tmp_path / "relocated_bundle"))
    if changed:
        with pytest.raises(ValueError, match="frozen"):
            export_browser(active)
    else:
        assert any(path.name == "export_audit.json" for path in export_browser(active))
