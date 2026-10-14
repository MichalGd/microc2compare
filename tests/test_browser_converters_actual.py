"""Real optional Linux converters: named vectors and chromosome-offset DB geometry."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sqlite3

import pytest

from microc2compare.browser_common import ExportContext
from microc2compare.browser_export import _profiles, _higlass
from microc2compare.browser_runtime import probe


def require_runtime():
    available = all(shutil.which(name) for name in ("bedGraphToBigWig", "clodius")) and importlib.util.find_spec("pyBigWig")
    if not available:
        if os.environ.get("M2C_REQUIRE_BROWSER_CONVERTERS") == "1":
            pytest.fail("The actual converter gate requires bedGraphToBigWig, clodius and pyBigWig")
        pytest.skip("Optional installed browser converters unavailable")


def assert_interval_database(path, raw, profile, paired):
    ordered = profile["chromosomes"]
    offsets, position = {}, 0
    for chrom, length in ordered:
        offsets[chrom] = position
        position += length
    expected = {line.strip() for line in raw.read_text().splitlines() if line.strip()}
    with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as database:
        database.row_factory = sqlite3.Row
        assert database.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        header = database.execute("SELECT chrom_names, chrom_sizes FROM tileset_info").fetchone()
        assert list(zip(header[0].split("\t"), map(int, header[1].split("\t")))) == [tuple(row) for row in ordered]
        records = database.execute("SELECT * FROM intervals").fetchall()
        assert {r["fields"].strip() for r in records} == expected
        assert len(records) == len(expected)
        for row in records:
            fields = row["fields"].strip().split("\t")
            index = database.execute("SELECT * FROM position_index WHERE id=?", (row["id"],)).fetchone()
            assert index is not None
            assert row["chrOffset"] == offsets[fields[0]]
            if paired:
                coordinates = (offsets[fields[0]] + int(fields[1]), offsets[fields[0]] + int(fields[2]),
                               offsets[fields[3]] + int(fields[4]), offsets[fields[3]] + int(fields[5]))
                assert tuple(row[key] for key in ("fromX", "toX", "fromY", "toY")) == coordinates
                # Rtree bounds may round outward at large coordinates, but must
                # contain the original exact interval and refer to the same row.
                assert index["rFromX"] <= coordinates[0] <= coordinates[1] <= index["rToX"]
                assert index["rFromY"] <= coordinates[2] <= coordinates[3] <= index["rToY"]
            else:
                start, end = offsets[fields[0]] + int(fields[1]), offsets[fields[0]] + int(fields[2])
                assert (row["startPos"], row["endPos"]) == (start, end)
                assert index["rStartPos"] <= start <= end <= index["rEndPos"]


def test_real_export_conversion_preserves_vectors_rectangles_and_axis_order(tmp_path):
    require_runtime()
    import pyBigWig
    evidence = probe()
    run = tmp_path / "fresh-native-metadata"
    scope = run / "00_metadata/masks/10.scope.json"
    scope.parent.mkdir(parents=True)
    samples = [{"sample_id": f"Fresh{i+1}", "condition": "first" if i < 2 else "second", "reference_id": "synthetic_coordinates_v2"} for i in range(5)]
    scope.write_text(json.dumps({"resolution": 10, "samples": [
        {"sample_id": s["sample_id"], "native_chromosomes": ["chr1", "chr2"] if i < 2 else ["chr2", "chr1"]} for i, s in enumerate(samples)]}))
    destination = os.environ.get("M2C_BROWSER_ACCEPTANCE_DIR")
    out = Path(destination).resolve() / "bundle" if destination else tmp_path / "export"
    out.mkdir(parents=True, exist_ok=False)
    lengths = {"chr1": 200000003, "chr2": 120000007}
    ctx = ExportContext(run, out, {"REFERENCE_ID": "synthetic_coordinates_v2", "EXPORT_UCSC": True, "EXPORT_HIGLASS": True}, samples, lengths)
    ctx.registered_orders = {10: ["chr1", "chr2"]}
    profiles, default = _profiles(ctx, lengths)
    assert len(profiles) == 2
    assert len({p["coord_system"] for p in profiles.values()}) == 2
    ctx.bedgraph("signal.bedgraph", [("chr1", 0, 10, 0), ("chr1", 10, 20, -2.5), ("chr2", 10, 20, 1.25)],
                 "Saved finite signal", "insulation", "individual", "Fresh1", "fresh saved records", 10, 20)
    ctx.bed("domains.bed", [("chr1", 10, 20, "exact_domain_1"), ("chr2", 0, 30, "exact_domain_2")],
            "Exact saved domains", "tads", "individual", "Fresh1", "fresh saved domains", 10, 20)
    ctx.interactions("loop", [{"feature_id": "native_trans_endpoint", "chrom1": "chr1", "start1": 10, "end1": 20,
                               "chrom2": "chr2", "start2": 15, "end2": 25}],
                     "Exact native endpoint geometry", "loops", "individual", "Fresh1", "fresh native calls", 10)
    _higlass(ctx, profiles, default)
    ctx.verify_inputs()
    assert all(t["higlass_status"] == "converted; not ingested or render-tested" for t in ctx.tracks)
    receipts = json.loads((out / "higlass/conversion_receipts.json").read_text())
    assert {r["filetype"] for r in receipts} == {"bigwig", "beddb", "bed2ddb"}
    assert all(not r["ingested"] and not r["render_tested"] for r in receipts)
    for receipt in receipts:
        path = out / receipt["path"]
        profile = profiles[receipt["profile_id"]]
        command = receipt["command"]
        raw = Path(command[1] if receipt["filetype"] == "bigwig" else command[-1])
        if receipt["filetype"] == "bigwig":
            assert "roundtrip_checked" in receipt["validation"]
            with pyBigWig.open(str(path)) as reader:
                assert reader.values("chr1", 0, 1) == [0.0]
                assert reader.values("chr1", 10, 11) == [-2.5]
                assert reader.values("chr2", 10, 11) == [1.25]
        else:
            assert_interval_database(path, raw, profile, receipt["filetype"] == "bed2ddb")
    commands = json.loads((out / "higlass/ingestion_commands.json").read_text())["commands"]
    profile_ids = set(profiles)
    assert {cmd[cmd.index("--assembly") + 1] for cmd in commands} == profile_ids
    evidence["format_acceptance"] = {"filetypes": ["bigwig", "beddb", "bed2ddb"], "converted_files": len(receipts),
                                     "exact_source_fields_verified": True, "chromosome_offsets_and_rtree_verified": True,
                                     "profile_ids": sorted(profile_ids), "higlass_ingested": False, "higlass_render_tested": False}
    if destination:
        (out.parent / "converter_acceptance.json").write_text(json.dumps(evidence, indent=2) + "\n")
