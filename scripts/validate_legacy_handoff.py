#!/usr/bin/env python3
"""Check legacy input discovery using new tiny fixtures; never read real inputs."""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE / "src"))
from microc2compare.contracts import ContractError, file_hash, write_table
from microc2compare.handoff import discover, write_drafts


def validate(destination):
    out = Path(destination).resolve()
    if out.exists():
        raise ValueError("Choose a new private validation directory; existing paths are never overwritten")
    out.mkdir(parents=True, exist_ok=False)
    root = out / "upstream"
    root.mkdir()
    rows = [{"sample": f"S{i}", "assay": "microc", "condition": "basal",
             "biological_replicate": str(i), "technical_replicate": "1"} for i in range(1, 5)]
    sheet = out / "samplesheet.csv"
    with sheet.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    recorded = [dict(row, biological_replicate="1", technical_replicate=row["biological_replicate"],
                     merge_group="historical_pool") for row in rows if row["sample"] != "S3"]
    write_table(root / "run_metadata/sample_manifest.tsv", recorded, tuple(recorded[0]))
    with (root / "run_metadata/technical_replicates.tsv").open("x", encoding="utf-8", newline="") as handle:
        for row in recorded:
            handle.write(f"historical_pool\t{row['sample']}\t/synthetic/{row['sample']}.valid.pairs.gz\n")
    for row in rows:
        path = root / row["sample"] / "04_matrices" / f"{row['sample']}.norm.mcool"
        path.parent.mkdir(parents=True)
        for resolution in (10000, 100000):
            bins = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, resolution],
                                 "end": [resolution, resolution * 2]})
            pixels = pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [5]})
            cooler.create_cooler(f"{path}::/resolutions/{resolution}", bins, pixels,
                                 ordered=True, assembly="mm39", mode="a")
            with h5py.File(path, "r+") as handle:
                weight = handle[f"resolutions/{resolution}/bins"].create_dataset("weight", data=np.ones(2))
                weight.attrs.update(converged=True, ignore_diags=2, tol=1e-5)
    before = {str(path): file_hash(path) for path in root.rglob("*") if path.is_file()}
    report = discover(sheet, root, legacy_reference_id="mm39", individual_matrices=True)
    assert len(report["units"]) == 4
    assert all(unit["status"] == "needs_review" and unit["selection"] == "individual" for unit in report["units"])
    by_id = {unit["sample_id"]: unit for unit in report["units"]}
    assert any("absent from" in issue for issue in by_id["S3"]["review_required"])
    assert by_id["S4"]["recorded_upstream_metadata"][0]["biological_replicate"] == "1"
    assert by_id["S4"]["biological_replicate"] == "4"
    assert len(by_id["S4"]["recorded_technical_membership"]) == 1
    output = write_drafts(report, out / "drafts")
    sidecar = json.loads((output / "provenance/S4.json").read_text())
    assert sidecar["draft"] and sidecar["count_units"] is None
    assert sidecar["handoff_evidence"]["recorded_upstream_metadata"][0]["merge_group"] == "historical_pool"
    missing = discover(sheet, root, individual_matrices=True)
    assert all(unit["status"] == "blocked" and not unit["reference_id"] for unit in missing["units"])
    ordinary = discover(sheet, root, legacy_reference_id="mm39")
    assert any(unit["status"] == "blocked" for unit in ordinary["units"])
    try:
        write_drafts(report, output)
    except ContractError:
        pass
    else:
        raise AssertionError("Existing drafts were overwritten")
    completed = subprocess.run([sys.executable, "-I", "-B", str(SOURCE / "scripts/prepare_microc2compare_inputs.py"),
                                "--samplesheet", str(sheet), "--upstream-results", str(root),
                                "--legacy-reference-id", "mm39", "--individual-matrices", "--json"],
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert len(json.loads(completed.stdout)["units"]) == 4
    assert before == {str(path): file_hash(path) for path in root.rglob("*") if path.is_file()}
    print("Legacy handoff acceptance passed: four individual inputs; old grouping and missing provenance preserved; input hashes unchanged.")
    print(f"Synthetic audit: {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="New private synthetic validation directory")
    args = parser.parse_args()
    try:
        validate(args.out)
    except (OSError, ValueError, AssertionError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
