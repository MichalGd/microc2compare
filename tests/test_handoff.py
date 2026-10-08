"""Real Cooler fixtures test source selection, evidence gaps and immutability."""
import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

h5py = pytest.importorskip("h5py")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
cooler = pytest.importorskip("cooler")

from microc2compare.contracts import ContractError, file_hash, read_table, write_table
from microc2compare.handoff import discover, write_drafts


def write_sheet(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def matrix(path, assembly="mm39", resolutions=(10000, 100000), converged=True, invalid_weight=False, floating=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    for resolution in resolutions:
        starts = np.arange(0, 300000, resolution)
        bins = pd.DataFrame({"chrom": "chr1", "start": starts, "end": starts + resolution})
        pixels = pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [8.5 if floating else 8]})
        uri = str(path) + f"::/resolutions/{resolution}"
        cooler.create_cooler(uri, bins, pixels, ordered=True, assembly=assembly, mode="a",
                             dtypes={"count": "float64" if floating else "int32"})
        with h5py.File(path, "r+") as handle:
            values = np.ones(len(bins))
            values[-1] = np.nan
            if invalid_weight:
                values[0] = np.inf
            weight = handle[f"resolutions/{resolution}/bins"].create_dataset("weight", data=values)
            if converged is not None:
                weight.attrs["converged"] = converged
            weight.attrs["ignore_diags"] = 2
            weight.attrs["tol"] = 1e-5
    return path


@pytest.fixture
def upstream(tmp_path):
    root = tmp_path / "upstream results"
    root.mkdir()
    rows = [
        {"sample": "C1_A", "assay": "microc", "reference_genome": "mouse", "condition": "control",
         "biological_replicate": "1", "technical_replicate": "1", "fastq_r1": "/reads/a, r1.fastq.gz", "fastq_r2": "/reads/a r2.fastq.gz"},
        {"sample": "C1_B", "assay": "microc", "reference_genome": "mm39", "condition": "control",
         "biological_replicate": "1", "technical_replicate": "2", "fastq_r1": "/reads/b r1.fastq.gz", "fastq_r2": "/reads/b r2.fastq.gz"},
        {"sample": "C2", "assay": "microc", "reference_genome": "mm39", "condition": "control",
         "biological_replicate": "2", "technical_replicate": "1", "fastq_r1": "/reads/c2 r1.fastq.gz", "fastq_r2": "/reads/c2 r2.fastq.gz"},
    ]
    sheet = tmp_path / "upstream samples.csv"
    write_sheet(sheet, rows)
    metadata = [{**row, "reference_id": "mm39", "merge_group": f"control_microc_B{row['biological_replicate']}_tech_merged"} for row in rows]
    write_table(root / "run_metadata/sample_manifest.tsv", metadata,
                ("sample", "assay", "reference_id", "condition", "biological_replicate", "technical_replicate", "merge_group"))
    write_table(root / "run_metadata/technical_replicates.tsv", metadata,
                ("merge_group", "sample", "reference_id", "assay", "condition", "biological_replicate"))
    write_table(root / "run_metadata/merge_manifest.tsv", metadata[:1], ("merge_group", "reference_id"))
    merged = matrix(root / "merged/control_microc_B1_tech_merged/04_matrices/control_microc_B1_tech_merged.norm.mcool")
    single = matrix(root / "C2/04_matrices/C2.norm.mcool")
    per_row = matrix(root / "C1_A/04_matrices/C1_A.norm.mcool")
    inventory = [{"category": "merged_norm_mcool", "path": str(merged)}, {"category": "sample_norm_mcool", "path": str(single)}]
    write_table(root / "final_report/microc2tracks_file_manifest.tsv", inventory, ("category", "path"))
    pairs = root / "C2/03_pairs/C2.valid.mapq30.pairs.gz"
    pairs.parent.mkdir()
    pairs.write_bytes(b"fixture-contact-file-not-read")
    hic = root / "C2/04_matrices/C2.norm.hic"
    hic.write_bytes(b"fixture-browser-file-not-read")
    return root, sheet, rows, merged, single, per_row


def test_selects_consolidated_and_single_maps_without_writes(upstream):
    root, sheet, _, merged, single, _ = upstream
    before = {str(path): file_hash(path) for path in root.rglob("*") if path.is_file()}
    report = discover(sheet, root)
    units = {row["sample_id"]: row for row in report["units"]}
    assert len(units) == 2
    assert units["control_microc_B1_tech_merged"]["upstream_rows"] == ["C1_A", "C1_B"]
    assert units["control_microc_B1_tech_merged"]["matrix_path"] == str(merged)
    assert units["C2"]["matrix_path"] == str(single)
    assert all(unit["matrix_origin"] == "file_manifest" and unit["status"] == "needs_review" for unit in units.values())
    assert units["C2"]["inspection"]["resolutions"]["10000"]["excluded_weight_bins"] == 1
    assert units["C2"]["inspection"]["resolutions"]["10000"]["balanced_inputs_available"]
    assert any(asset["path"].endswith(".pairs.gz") for asset in units["C2"]["related_outputs"])
    assert any(asset["path"].endswith(".hic") for asset in units["C2"]["related_outputs"])
    assert before == {str(path): file_hash(path) for path in root.rglob("*") if path.is_file()}


def test_missing_merged_matrix_never_falls_back_to_technical_row(upstream):
    root, sheet, _, merged, _, per_row = upstream
    merged.unlink()
    unit = next(unit for unit in discover(sheet, root)["units"] if len(unit["upstream_rows"]) == 2)
    assert unit["status"] == "blocked" and unit["matrix_path"] != str(per_row)
    assert unit["inspection"]["status"] == "missing"


@pytest.mark.parametrize("converged,invalid,floating,expected", [
    (False, False, False, "convergence"), (None, False, False, "convergence"),
    (True, True, False, "infinite"), (True, False, True, "integer storage")])
def test_failed_unknown_or_invalid_matrix_evidence(upstream, converged, invalid, floating, expected):
    root, sheet, _, _, single, _ = upstream
    single.unlink()
    matrix(single, converged=converged, invalid_weight=invalid, floating=floating)
    unit = next(unit for unit in discover(sheet, root)["units"] if unit["sample_id"] == "C2")
    assert unit["status"] == "blocked"
    assert any(expected in issue for issue in unit["errors"])
    assert unit["inspection"]["resolutions"]["10000"]["converged"] is converged


def test_missing_resolution_is_listed_without_creating_it(upstream):
    root, sheet, _, _, single, _ = upstream
    before = file_hash(single)
    unit = discover(sheet, root, resolutions=(5000, 10000))["units"][0]
    assert unit["status"] == "blocked"
    assert unit["inspection"]["available_resolutions"] == [10000, 100000]
    assert not unit["inspection"]["resolutions"]["5000"]["available"]
    assert file_hash(single) == before


def test_partial_sheet_cannot_masquerade_as_complete_merge(upstream):
    root, sheet, rows, *_ = upstream
    write_sheet(sheet, rows[1:])
    unit = next(unit for unit in discover(sheet, root)["units"] if unit["upstream_rows"] == ["C1_B"])
    assert unit["status"] == "blocked"
    assert any("membership differs" in issue for issue in unit["errors"])


def test_stale_sheet_reference_and_matrix_assembly_are_flagged(upstream):
    root, sheet, rows, _, single, _ = upstream
    rows[-1]["reference_genome"] = "human"
    write_sheet(sheet, rows)
    unit = next(unit for unit in discover(sheet, root)["units"] if unit["sample_id"] == "C2")
    assert any("reference differs" in issue for issue in unit["errors"])
    assert any("assembly metadata" in issue for issue in unit["errors"])


def test_drafts_preserve_gaps_and_refuse_overwriting(upstream, tmp_path):
    root, sheet, _, merged, single, _ = upstream
    hashes = [file_hash(merged), file_hash(single)]
    output = write_drafts(discover(sheet, root, project="study"), tmp_path / "drafts")
    assert (output / "samples.csv").read_text().splitlines()[0].startswith("sample_id,project,")
    samples = read_table(output / "samples.csv")
    assert all(not row["subject_id"] and row["project"] == "study" for row in samples)
    membership = read_table(output / "technical_membership.tsv")
    assert len(membership) == 3 and all(not row["library_id"] and not row["source_checksum"] for row in membership)
    assert membership[0]["fastq_r1"] == "/reads/a, r1.fastq.gz"
    assert read_table(output / "comparisons.tsv") == []
    provenance = json.loads((output / "provenance/C2.json").read_text())
    assert provenance["count_units"] is None and provenance["filters"]["mapq"] is None
    assert provenance["library_ids"] == [] and provenance["dedup_scope"] == "unknown"
    assert provenance["balancing"]["10000"]["converged"] is True
    assert provenance["balancing"]["10000"]["recipe"]["ignore_diags"] == 2
    with pytest.raises(ContractError, match="already exists"):
        write_drafts(discover(sheet, root), output)
    assert hashes == [file_hash(merged), file_hash(single)]


def test_drafts_cannot_be_written_inside_upstream(upstream):
    root, sheet, *_ = upstream
    destination = root / "new drafts"
    with pytest.raises(ContractError, match="outside upstream"):
        write_drafts(discover(sheet, root), destination)
    assert not destination.exists()


def test_explicit_subject_and_invalid_lane_relationships(upstream):
    root, sheet, rows, *_ = upstream
    for row in rows[:2]:
        row.update(subject_id="culture1", library_id="same_library", source_kind="lane", dedup_scope="per_row", dedup_group_id="")
    rows[-1]["subject_id"] = "culture2"
    write_sheet(sheet, rows)
    unit = next(unit for unit in discover(sheet, root)["units"] if len(unit["upstream_rows"]) == 2)
    assert unit["subject_id"] == "culture1"
    assert any("split-library" in issue for issue in unit["errors"])


def test_missing_inspection_dependency_is_visible(upstream, monkeypatch):
    root, sheet, *_ = upstream
    monkeypatch.setitem(sys.modules, "h5py", None)
    report = discover(sheet, root)
    assert all(unit["inspection"]["status"] == "inspection_unavailable" and unit["status"] == "blocked" for unit in report["units"])


def test_legacy_reference_requires_explicit_evidence(upstream):
    root, sheet, rows, *_ = upstream
    for path in (root / "run_metadata").iterdir():
        path.unlink()
    for row in rows:
        row.pop("reference_genome")
    write_sheet(sheet, rows)
    assert all(unit["status"] == "blocked" and not unit["reference_id"] for unit in discover(sheet, root)["units"])
    assert all(unit["reference_id"] == "mm39" and unit["status"] == "needs_review"
               for unit in discover(sheet, root, legacy_reference_id="mouse")["units"])


def test_other_run_manifest_path_is_not_selected(upstream, tmp_path):
    root, sheet, _, _, single, _ = upstream
    foreign = matrix(tmp_path / "another run/C2.norm.mcool", assembly="hg38")
    write_table(root / "final_report/microc2tracks_file_manifest.tsv",
                [{"category": "sample_norm_mcool", "path": str(foreign)}], ("category", "path"))
    unit = next(unit for unit in discover(sheet, root)["units"] if unit["sample_id"] == "C2")
    assert unit["matrix_path"] == str(single) and unit["matrix_origin"] == "standard_layout"


def test_cli_defaults_to_read_only_and_handles_malformed_input(upstream):
    root, sheet, rows, *_ = upstream
    command = [sys.executable, str(Path(__file__).resolve().parents[1] / "scripts/prepare_microc2compare_inputs.py"),
               "--samplesheet", str(sheet), "--upstream-results", str(root), "--json"]
    before = {str(path) for path in root.rglob("*")}
    completed = subprocess.run(command, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert len(json.loads(completed.stdout)["units"]) == 2
    assert before == {str(path) for path in root.rglob("*")}
    write_sheet(sheet, rows + rows[:1])
    failed = subprocess.run(command, capture_output=True, text=True)
    assert failed.returncode == 2 and "Duplicate upstream sample" in failed.stderr


def test_mixed_reference_export_requires_selection(upstream, tmp_path):
    root, sheet, rows, _, single, _ = upstream
    rows[-1]["reference_genome"] = "human"
    write_sheet(sheet, rows)
    for name in ("sample_manifest.tsv", "technical_replicates.tsv"):
        path = root / "run_metadata" / name
        metadata = read_table(path)
        for row in metadata:
            if row["sample"] == "C2":
                row["reference_id"] = "hg38"
        write_table(path, metadata, tuple(metadata[0]))
    with h5py.File(single, "r+") as handle:
        for resolution in (10000, 100000):
            handle[f"resolutions/{resolution}"].attrs["genome-assembly"] = "hg38"
    report = discover(sheet, root)
    assert {unit["reference_id"] for unit in report["units"]} == {"mm39", "hg38"}
    with pytest.raises(ContractError, match="one reference"):
        write_drafts(report, tmp_path / "mixed drafts")
    assert not (tmp_path / "mixed drafts").exists()
    selected = discover(sheet, root, reference_id="GRCh38")
    assert len(selected["units"]) == 1 and selected["units"][0]["status"] == "needs_review"
    assert selected["excluded_upstream_samples"] == ["C1_A", "C1_B"]
    write_drafts(selected, tmp_path / "human drafts")


def test_relative_file_manifest_and_chromosome_mismatch(upstream):
    root, sheet, _, merged, single, _ = upstream
    inventory = [{"category": "merged_norm_mcool", "path": str(merged.relative_to(root.parent))},
                 {"category": "sample_norm_mcool", "path": str(single.relative_to(root))}]
    write_table(root / "final_report/microc2tracks_file_manifest.tsv", inventory, ("category", "path"))
    assert all(unit["matrix_origin"] == "file_manifest" for unit in discover(sheet, root)["units"])
    with h5py.File(single, "r+") as handle:
        handle["resolutions/10000/chroms/name"][0] = b"chrX"
    assert all(any("Chromosome names/order/length differ" in issue for issue in unit["errors"])
               for unit in discover(sheet, root)["units"])


def test_shared_subject_or_library_is_not_extra_replication(upstream):
    root, sheet, rows, *_ = upstream
    for row in rows:
        row["subject_id"] = "one_culture"
        row["library_id"] = "one_library"
    write_sheet(sheet, rows)
    report = discover(sheet, root)
    assert all(unit["status"] == "blocked" for unit in report["units"])
    assert all(any("same declared subject" in issue for issue in unit["errors"]) for unit in report["units"])
    assert all(any("declared library" in issue for issue in unit["errors"]) for unit in report["units"])


def test_different_merge_groups_cannot_silently_pick_one_map(upstream):
    root, sheet, rows, *_ = upstream
    for path in (root / "run_metadata").iterdir():
        path.unlink()
    rows[0]["merge_group"] = "group_A"
    rows[1]["merge_group"] = "group_B"
    write_sheet(sheet, rows)
    unit = next(unit for unit in discover(sheet, root)["units"] if len(unit["upstream_rows"]) == 2)
    assert unit["status"] == "blocked" and not unit["matrix_path"]
    assert any("different upstream merge groups" in issue for issue in unit["errors"])


def test_preview_does_not_read_contact_pixel_values(upstream, monkeypatch):
    root, sheet, *_ = upstream
    original = h5py.Dataset.__getitem__

    def guarded(dataset, key):
        if "/pixels/" in dataset.name:
            raise AssertionError("Preview must not scan contact pixels")
        return original(dataset, key)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded)
    report = discover(sheet, root)
    assert all(unit["status"] == "needs_review" for unit in report["units"])
    assert all(unit["inspection"]["resolutions"]["10000"]["potential_uses"] for unit in report["units"])


@pytest.mark.parametrize("conflict", ["subjects", "independent_libraries"])
def test_known_biological_or_library_conflict_blocks_consolidation(upstream, conflict):
    root, sheet, rows, *_ = upstream
    if conflict == "subjects":
        rows[0]["subject_id"] = "donor_A"
        rows[1]["subject_id"] = "donor_B"
        expected = "different declared subject"
    else:
        for row in rows[:2]:
            row.update(library_id="one_library", source_kind="independent_library")
        expected = "reused for sources declared independent"
    write_sheet(sheet, rows)
    unit = next(unit for unit in discover(sheet, root)["units"] if len(unit["upstream_rows"]) == 2)
    assert unit["status"] == "blocked" and any(expected in issue for issue in unit["errors"])
