"""Suggest biological-unit inputs from microc2tracks without changing its files.

Discovery uses the standard library. Matrix inspection imports h5py/numpy only
when available. This is an input inventory, not scientific preflight acceptance.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from .contracts import IDENTIFIER, SAMPLE_FIELDS, atomic_json, read_table, require, write_table

ALIASES = {"mouse": "mm39", "mm39": "mm39", "GRCm39": "mm39",
           "human": "hg38", "hg38": "hg38", "GRCh38": "hg38"}
MEMBERSHIP_FIELDS = ("sample_id", "library_id", "source_id", "source_kind",
                     "dedup_group_id", "dedup_scope", "source_checksum",
                     "upstream_sample", "fastq_r1", "fastq_r2")
COMPARISON_FIELDS = ("comparison_id", "case_condition", "reference_condition",
                     "design", "inference", "modules")


def canonical_reference(value):
    value = str(value or "").strip()
    if value.lower() in ("unknown", "undefined", "none", "null", "n/a"):
        return ""
    return ALIASES.get(value, value)


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.-")


def optional_table(path, fields):
    # Upstream can emit an empty file when a report has no rows.
    return read_table(path, fields) if path.is_file() and path.stat().st_size else []


def upstream_table(path, fields, technical=False):
    """Read current or known legacy metadata without inventing reference evidence."""
    if not path.is_file() or not path.stat().st_size:
        return [], []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first_line = handle.readline()
        first = next(csv.reader([first_line], delimiter="\t" if "\t" in first_line else ","))
        handle.seek(0)
        if technical and "sample" not in first and "merge_group" not in first:
            # The original shell helper wrote group, sample, valid-pairs path.
            records = []
            for line, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
                if not values:
                    continue
                require(len(values) == 3 and IDENTIFIER.fullmatch(values[0])
                        and values[1] and values[2].endswith((".pairs", ".pairs.gz")),
                        f"{path}:{line}: expected legacy merge_group, sample, valid-pairs path")
                records.append(dict(zip(("merge_group", "sample", "pairs_path"), values)))
            return [dict(row, reference_id="") for row in records], [
                f"Legacy headerless technical membership read from {path}; reference identity is not recorded there."]
    records = read_table(path, fields)
    missing_reference = any(not canonical_reference(row.get("reference_id")) for row in records)
    notes = [f"Upstream reference identity is missing in {path}; the sheet/explicit legacy declaration needs provenance review."] if missing_reference else []
    return [dict(row, reference_id=row.get("reference_id", "")) for row in records], notes


def scalar(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "tolist"):
        return scalar(value.tolist())
    if isinstance(value, (list, tuple)):
        return [scalar(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def inspect_matrix(path, resolutions):
    """Read Cooler headers and requested bin weights, never contact pixels."""
    if not path or not Path(path).is_file():
        return {"status": "missing", "available_resolutions": [], "resolutions": {},
                "errors": ["Required biological-unit matrix is missing."]}
    try:
        import h5py
        import numpy as np
    except ImportError:
        return {"status": "inspection_unavailable", "available_resolutions": [],
                "resolutions": {}, "errors": ["Use the checked core Python with h5py/numpy to inspect matrices."]}
    result = {"status": "inspected", "available_resolutions": [], "resolutions": {}, "errors": []}
    try:
        with h5py.File(path, "r") as handle:
            groups = {}
            if "resolutions" in handle:
                for key in handle["resolutions"]:
                    if key.isdigit():
                        groups[int(key)] = f"/resolutions/{key}"
            elif handle.attrs.get("bin-size") is not None:
                groups[int(handle.attrs["bin-size"])] = "/"
            result["available_resolutions"] = sorted(groups)
            for resolution in resolutions:
                info = {"available": resolution in groups, "raw_integer_storage": False,
                        "balanced_inputs_available": False, "converged": None,
                        "weight_column": "weight", "convention": None, "recipe": {}, "errors": []}
                result["resolutions"][str(resolution)] = info
                if resolution not in groups:
                    info["errors"].append(f"Missing {resolution}-bp resolution; request it upstream.")
                    continue
                group = handle[groups[resolution]]
                info["uri"] = str(Path(path).resolve()) + "::" + groups[resolution]
                info["assembly"] = scalar(group.attrs.get("genome-assembly", group.attrs.get("assembly", "")))
                info["storage_mode"] = scalar(group.attrs.get("storage-mode", "symmetric-upper"))
                required = ("chroms/name", "chroms/length", "bins/chrom", "bins/start", "bins/end",
                            "pixels/count", "pixels/bin1_id", "pixels/bin2_id",
                            "indexes/chrom_offset", "indexes/bin1_offset")
                if (scalar(group.attrs.get("format", "")) != "HDF5::Cooler"
                        or group.attrs.get("bin-size") != resolution
                        or any(key not in group for key in required)):
                    info["errors"].append("Invalid or inconsistent Cooler layout; scientific preflight is required.")
                    continue
                info["chromosomes"] = [[scalar(name), int(length)] for name, length in
                                       zip(group["chroms/name"][:], group["chroms/length"][:])]
                info["chromosome_signature"] = hashlib.sha256(
                    json.dumps(info["chromosomes"], separators=(",", ":")).encode()).hexdigest()
                info["raw_integer_storage"] = group["pixels/count"].dtype.kind in "iu"
                info["n_bins"] = len(group["bins/start"])
                info["n_stored_pixels"] = len(group["pixels/count"])
                if not info["raw_integer_storage"]:
                    info["errors"].append("pixels/count is not integer storage; do not substitute balanced values.")
                if info["storage_mode"] not in ("symmetric-upper", "square"):
                    info["errors"].append("Unsupported matrix storage mode.")
                if "bins/weight" not in group:
                    info["errors"].append("Stored balancing weights are missing; no downstream repair is attempted.")
                    continue
                weights = group["bins/weight"]
                if weights.ndim != 1 or len(weights) != info["n_bins"] or weights.dtype.kind not in "fiu":
                    info["errors"].append("Invalid balancing-weight column shape or type.")
                    continue
                attributes = {key: scalar(value) for key, value in weights.attrs.items()}
                convergence = attributes.get("converged")
                if isinstance(convergence, (bool, int)) and convergence in (0, 1):
                    info["converged"] = bool(convergence)
                divisive = attributes.get("divisive_weights", False)
                if not isinstance(divisive, (bool, int)) or divisive not in (0, 1):
                    info["errors"].append("Unrecognized weight convention metadata.")
                else:
                    info["convention"] = "divisive" if divisive else "multiplicative"
                # Keep recipe parameters separate from sample-specific fit results.
                recipe_keys = ("cis_only", "trans_only", "ignore_diags", "mad_max", "min_count", "min_nnz", "tol")
                info["recipe"] = {key: attributes[key] for key in recipe_keys if key in attributes}
                info["weight_attributes"] = attributes
                valid = excluded = invalid = 0
                for start in range(0, len(weights), 100000):
                    values = weights[start:start + 100000]
                    valid += int((np.isfinite(values) & (values > 0)).sum())
                    excluded += int(np.isnan(values).sum())
                    invalid += int((np.isinf(values) | (np.isfinite(values) & (values <= 0))).sum())
                info.update(valid_weight_bins=valid, excluded_weight_bins=excluded, invalid_weight_bins=invalid)
                if invalid or not valid:
                    info["errors"].append("Weights contain infinite/nonpositive values or no usable bins.")
                if info["converged"] is not True:
                    info["errors"].append("Balancing convergence failed or is unknown; review actual upstream evidence.")
                info["balanced_inputs_available"] = not info["errors"]
    except (OSError, ValueError, TypeError, KeyError) as error:
        result["status"] = "invalid"
        result["errors"].append(f"Cannot inspect matrix: {error}")
    for resolution, info in result["resolutions"].items():
        info["potential_uses"] = []
        if info["raw_integer_storage"]:
            info["potential_uses"].append("Raw contact summaries/HiCRep (count/reference evidence and preflight still required)")
        if info["balanced_inputs_available"]:
            info["potential_uses"].append("Balanced P(s) and feature analysis at a compatible configured grid")
            if resolution == "10000":
                info["potential_uses"].append("Default insulation/loop/aggregate input grid; feature/design gates still apply")
            if resolution == "100000":
                info["potential_uses"].append("Default compartment/saddle input grid; registered phasing track still required")
    return result


def manifest_paths(value, root):
    """Resolve relative upstream paths and reject paths from another run root."""
    path = Path(value)
    candidates = [path] if path.is_absolute() else [root / path, root.parent / path, Path.cwd() / path]
    return {p.resolve() for p in candidates if p.resolve().is_relative_to(root) and p.is_file()}


def choose_matrix(root, stem, merged, inventory):
    folder = root / "merged" / stem if merged else root / stem
    expected = folder / "04_matrices" / f"{stem}.norm.mcool"
    category = "merged_norm_mcool" if merged else "sample_norm_mcool"
    listed = [row for row in inventory if row["category"] == category
              and Path(row["path"]).name == expected.name]
    candidates = set()
    for row in listed:
        candidates.update(manifest_paths(row["path"], root))
    if len(candidates) > 1:
        return "", "ambiguous", ["Multiple manifest matrices match this unit; select a reviewed upstream run."]
    if candidates:
        return str(next(iter(candidates))), "file_manifest", []
    return str(expected.resolve()), "standard_layout", []


def related_outputs(root, stem, merged, sources, inventory):
    records = {}
    for owner, is_merged in [(stem, merged)] + [(row["sample"], False) for row in sources]:
        prefix = "merged_" if is_merged else "sample_"
        for row in inventory:
            if row["category"].startswith(prefix) and Path(row["path"]).name.startswith(owner + "."):
                for path in manifest_paths(row["path"], root):
                    records[str(path)] = {"category": row["category"], "owner": owner, "path": str(path)}
        folder = root / "merged" / owner if is_merged else root / owner
        for stage in ("01_qc", "03_pairs", "04_matrices", "05_downstream", "logs"):
            directory = folder / stage
            if not directory.is_dir():
                continue
            for path in sorted(directory.glob(owner + ".*")):
                if not path.is_file() or path.suffix not in (".gz", ".hic", ".json", ".html", ".tsv", ".txt", ".bed", ".bedgraph"):
                    continue
                records.setdefault(str(path.resolve()), {"category": stage, "owner": owner, "path": str(path.resolve())})
    return sorted(records.values(), key=lambda row: row["path"])


def discover(samplesheet, upstream_results, project="project", resolutions=(10000, 100000),
             reference_id=None, legacy_reference_id=None, individual_matrices=False):
    sheet, root = Path(samplesheet).resolve(), Path(upstream_results).resolve()
    require(root.is_dir(), f"Missing microc2tracks results directory: {root}")
    require(IDENTIFIER.fullmatch(project), "Use a safe project identifier: letters, digits, _, - or .")
    require(resolutions and all(isinstance(value, int) and value > 0 for value in resolutions), "Resolutions must be positive integers")
    rows = read_table(sheet, ("sample", "assay", "condition", "biological_replicate"))
    require(rows, "The upstream sample sheet contains no samples")
    rows = [{key: value.strip() for key, value in row.items()} for row in rows]
    identifiers = [row["sample"] for row in rows]
    require(len(set(identifiers)) == len(identifiers), "Duplicate upstream sample IDs; each FASTQ source needs its upstream unique sample name")
    require(all(IDENTIFIER.fullmatch(value) for value in identifiers), "Unsafe upstream sample identifier")
    metadata, metadata_notes = upstream_table(root / "run_metadata/sample_manifest.tsv", ("sample", "merge_group"))
    require(len({row["sample"] for row in metadata}) == len(metadata), "Duplicate samples in upstream sample_manifest.tsv")
    by_sample = {row["sample"]: row for row in metadata}
    technical, technical_notes = upstream_table(root / "run_metadata/technical_replicates.tsv", ("sample", "merge_group"), technical=True)
    merge_rows, merge_notes = upstream_table(root / "run_metadata/merge_manifest.tsv", ("merge_group",))
    inventory = optional_table(root / "final_report/microc2tracks_file_manifest.tsv", ("category", "path"))
    groups, excluded = defaultdict(list), []
    for original in rows:
        row = dict(original)
        require(row["condition"] and row["biological_replicate"], "Condition and biological_replicate cannot be empty")
        row["assay"] = row["assay"].lower()
        require(row["assay"] in ("microc", "hic"), f"Unsupported assay for {row['sample']}")
        upstream = by_sample.get(row["sample"], {})
        declared = canonical_reference(row.get("reference_genome") or row.get("reference_id"))
        recorded = canonical_reference(upstream.get("reference_id"))
        row["reference_id"] = declared or recorded or canonical_reference(legacy_reference_id)
        row["issues"] = []
        row["review_issues"] = []
        row["recorded_upstream_metadata"] = upstream
        if metadata and not upstream:
            row["review_issues"].append("Sample is absent from the existing sample manifest; locate its actual processing provenance.")
        if declared and recorded and declared != recorded:
            row["issues"].append("Sample-sheet reference differs from the upstream run manifest.")
        for item in technical:
            recorded_member_reference = canonical_reference(item["reference_id"])
            if item["sample"] == row["sample"] and recorded_member_reference and recorded_member_reference != row["reference_id"]:
                row["issues"].append("Technical membership reference differs from the selected sample reference.")
        for field in ("assay", "condition", "biological_replicate", "technical_replicate"):
            if upstream.get(field) and row.get(field) and upstream[field] != row[field]:
                target = row["review_issues"] if individual_matrices and field in ("biological_replicate", "technical_replicate") else row["issues"]
                target.append(f"Sample sheet and run manifest disagree on {field}." +
                              (" Individual-matrix selection uses the sheet's declared unit; verify the experimental identity." if target is row["review_issues"] else ""))
        if row.get("merge_group") and upstream.get("merge_group") and row["merge_group"] != upstream["merge_group"]:
            target = row["review_issues"] if individual_matrices else row["issues"]
            target.append("Sample-sheet merge_group differs from upstream metadata.")
        row["merge_group"] = upstream.get("merge_group") or row.get("merge_group") or safe_name(
            f"{row['condition']}_{row['assay']}_B{row['biological_replicate']}_tech_merged")
        if individual_matrices:
            if upstream.get("merge_group"):
                row["review_issues"].append(f"Historical merge group {upstream['merge_group']} is retained as evidence; its pooled matrix is not selected.")
            row["merge_group"] = safe_name(f"individual_{row['sample']}")
        require(IDENTIFIER.fullmatch(row["merge_group"]), "Unsafe or empty upstream merge_group")
        if reference_id and row["reference_id"] != canonical_reference(reference_id):
            excluded.append(row["sample"])
            continue
        groups[(row["reference_id"], row["assay"], row["condition"], row["biological_replicate"])].append(row)
    require(groups, "No samples match the requested reference")
    require(not individual_matrices or all(len(sources) == 1 for sources in groups.values()),
            "--individual-matrices requires exactly one sheet source per biological unit; technical inputs still need upstream merging")
    units = []
    for (reference, assay, condition, replicate), sources in sorted(groups.items()):
        sources.sort(key=lambda row: row["sample"])
        merge_names = {row["merge_group"] for row in sources}
        merged = len(sources) > 1
        stem = sources[0]["merge_group"] if merged else sources[0]["sample"]
        errors = [issue for row in sources for issue in row["issues"]]
        if len(merge_names) != 1:
            errors.append("Sources for this biological unit belong to different upstream merge groups; a consolidated map is required.")
        members = {row["sample"] for row in sources}
        for merge_name in merge_names:
            if individual_matrices:
                continue
            recorded_members = {row["sample"] for row in technical + metadata if row["merge_group"] == merge_name}
            if recorded_members and recorded_members != members:
                errors.append("Upstream merge membership differs from the selected sheet; do not import a partial or mixed group.")
            for item in merge_rows + technical + metadata:
                if item["merge_group"] == merge_name:
                    if canonical_reference(item["reference_id"]) and canonical_reference(item["reference_id"]) != reference:
                        errors.append("Upstream merge group crosses reference identities.")
                    for field, expected in (("assay", assay), ("condition", condition), ("biological_replicate", replicate)):
                        if item.get(field) and item[field] != expected:
                            errors.append(f"Upstream merge group crosses {field} identities.")
        if not reference:
            errors.append("Reference identity is unknown; supply actual run metadata or --legacy-reference-id explicitly.")
        path, origin, selection_errors = choose_matrix(root, stem, merged, inventory)
        errors += selection_errors
        if len(merge_names) != 1:
            path = ""
        inspection = inspect_matrix(path, sorted(set(resolutions)))
        errors += inspection["errors"]
        for info in inspection["resolutions"].values():
            errors += info["errors"]
            assembly = canonical_reference(info.get("assembly"))
            if assembly and reference and assembly != reference:
                errors.append("Matrix assembly metadata differs from the selected reference.")
        subject_values = {row.get("subject_id", "") for row in sources}
        batch_values = {row.get("batch", "") for row in sources}
        if len(subject_values - {""}) > 1:
            errors.append("Upstream group contains different declared subject identities; do not merge biological replicates.")
        review = ["Confirm actual independent biological units and the intended contrasts.",
                  "Supply source/library identities, source checksums, filters and deduplication evidence."]
        review += metadata_notes + technical_notes + merge_notes
        review += [issue for row in sources for issue in row["review_issues"]]
        if any(not canonical_reference(info.get("assembly")) for info in inspection["resolutions"].values()):
            review.append("Matrix assembly metadata is missing/unknown; verify reference provenance and selected chromosome/bin coordinates during preflight.")
        if len(subject_values) != 1 or "" in subject_values:
            review.append("Supply a consistent subject_id (actual donor/culture identity) for this unit.")
        if len(batch_values) > 1:
            review.append("Resolve inconsistent batch labels within this biological unit.")
        if merged:
            review.append("Merged inputs were deduplicated per upstream row; same-library lanes/resequencing need joint-library evidence or upstream reprocessing.")
        libraries = defaultdict(list)
        for row in sources:
            if row.get("library_id"):
                libraries[row["library_id"]].append(row)
        for library, members_in_library in libraries.items():
            if len(members_in_library) > 1 and any(row.get("source_kind") == "independent_library" for row in members_in_library):
                errors.append(f"Library {library}: one library ID was reused for sources declared independent libraries.")
            if len(members_in_library) > 1 and any(row.get("source_kind") in ("lane", "resequencing") for row in members_in_library):
                if (any(row.get("dedup_scope") != "joint_library" for row in members_in_library)
                        or not all(row.get("dedup_group_id") for row in members_in_library)
                        or len({row.get("dedup_group_id") for row in members_in_library}) != 1):
                    errors.append(f"Library {library}: split-library sources lack consistent joint-library deduplication evidence.")
        units.append({"sample_id": stem, "reference_id": reference, "assay": assay, "condition": condition,
                      "biological_replicate": replicate, "subject_id": next(iter(subject_values)) if len(subject_values) == 1 else "",
                      "batch": next(iter(batch_values)) if len(batch_values) == 1 else "",
                      "upstream_rows": sorted(members), "sources": sources, "selection": "merged" if merged else "individual",
                      "matrix_path": path, "matrix_origin": origin, "inspection": inspection,
                      "related_outputs": related_outputs(root, stem, merged, sources, inventory),
                      "recorded_upstream_metadata": [row["recorded_upstream_metadata"] for row in sources],
                      "recorded_technical_membership": [row for row in technical if row["sample"] in members],
                      "errors": sorted(set(errors)), "review_required": review})
    id_counts = defaultdict(int)
    path_counts = defaultdict(int)
    group_keys = defaultdict(set)
    library_owners = defaultdict(set)
    subject_observations = defaultdict(set)
    for unit in units:
        id_counts[unit["sample_id"]] += 1
        if unit["matrix_path"]:
            path_counts[unit["matrix_path"]] += 1
        for row in unit["sources"]:
            group_keys[row["merge_group"]].add((unit["reference_id"], unit["assay"], unit["condition"], unit["biological_replicate"]))
            if row.get("library_id"):
                library_owners[row["library_id"]].add(unit["sample_id"])
        if unit["subject_id"]:
            subject_observations[(unit["reference_id"], unit["subject_id"], unit["condition"])].add(unit["sample_id"])
    for unit in units:
        if unit["sample_id"] in ("mask", "labels") or unit["sample_id"].startswith("pool__"):
            unit["errors"].append("Selected sample ID is reserved by microc2compare; choose an explicit safe identity before import.")
        if id_counts[unit["sample_id"]] > 1 or path_counts[unit["matrix_path"]] > 1:
            unit["errors"].append("One sample ID or matrix was selected for multiple biological units.")
        if any(len(group_keys[row["merge_group"]]) > 1 for row in unit["sources"]):
            unit["errors"].append("A merge-group name represents multiple biological identities; correct upstream grouping.")
        if any(len(library_owners[row["library_id"]]) > 1 for row in unit["sources"] if row.get("library_id")):
            unit["errors"].append("One declared library belongs to multiple biological observations.")
        if unit["subject_id"] and len(subject_observations[(unit["reference_id"], unit["subject_id"], unit["condition"])]) > 1:
            unit["errors"].append("Multiple rows represent the same declared subject/condition observation; consolidate upstream.")
        for resolution, info in unit["inspection"]["resolutions"].items():
            signatures = {other["inspection"]["resolutions"].get(resolution, {}).get("chromosome_signature")
                          for other in units if other["reference_id"] == unit["reference_id"]}
            signatures.discard(None)
            if len(signatures) > 1:
                axis = dict(info.get("chromosomes", []))
                for other in units:
                    if other["reference_id"] != unit["reference_id"]: continue
                    other_axis = dict(other["inspection"]["resolutions"].get(resolution, {}).get("chromosomes", []))
                    if axis and other_axis and not axis.keys() & other_axis.keys():
                        unit["errors"].append(f"Chromosome names/order/length differ with no shared chromosomes between {resolution}-bp inputs.")
                    if any(axis[chrom] != other_axis[chrom] for chrom in axis.keys() & other_axis.keys()):
                        unit["errors"].append(f"Shared chromosome lengths differ between {resolution}-bp inputs.")
                unit["review_required"].append(f"Chromosome sets/order differ between {resolution}-bp inputs; declare common analysis views. Preflight must verify their exact bin coordinates. Stored weights retain upstream balancing scope.")
        unit["errors"] = sorted(set(unit["errors"]))
        unit["status"] = "blocked" if unit["errors"] else "needs_review"
    return {"schema_version": 1, "status": "draft_inventory_not_preflight_acceptance",
            "project": project, "samplesheet": str(sheet), "upstream_results": str(root),
            "requested_resolutions": sorted(set(resolutions)), "excluded_upstream_samples": excluded,
            "individual_matrices": bool(individual_matrices),
            "units": units,
            "limitations": ["No contact-pixel scan, full bin-axis comparison or scientific preflight was performed.",
                            "Matrix existence and metadata cannot establish independent biological replication.",
                            "No input files, counts, weights or upstream configurations are changed."]}


def summary(report):
    lines = ["microc2tracks -> microc2compare input inventory (DRAFT)",
             "Status is discovery/review status, not statistical eligibility.",
             "sample_id\tcondition\tbiological_replicate\treference\tsources\tstatus\tmatrix"]
    for unit in report["units"]:
        lines.append("\t".join([unit["sample_id"], unit["condition"], unit["biological_replicate"],
                                unit["reference_id"] or "UNKNOWN", ";".join(unit["upstream_rows"]),
                                unit["status"], unit["matrix_path"] or "AMBIGUOUS"]))
        for resolution, info in unit["inspection"]["resolutions"].items():
            lines.append(f"  {resolution} bp: raw_integer_storage={info['raw_integer_storage']}; "
                         f"balanced_inputs_available={info['balanced_inputs_available']}; converged={info['converged']}")
            for use in info.get("potential_uses", []):
                lines.append("    Potential use: " + use)
        lines.append("  Available resolutions: " + ", ".join(map(str, unit["inspection"]["available_resolutions"])))
        for issue in unit["errors"] + unit["review_required"]:
            lines.append("  Review: " + issue)
        for asset in unit["related_outputs"]:
            lines.append(f"  Related {asset['category']}: {asset['path']}")
    if report["excluded_upstream_samples"]:
        lines.append("Excluded by reference selection: " + ", ".join(report["excluded_upstream_samples"]))
    return "\n".join(lines)


def write_drafts(report, destination):
    """Write only into a new directory outside upstream results."""
    output = Path(destination).resolve()
    require(not output.is_relative_to(Path(report["upstream_results"])), "Draft destination must be outside upstream results")
    require(not output.exists(), "Draft destination already exists; choose a new directory (no overwriting)")
    require(len({unit["reference_id"] for unit in report["units"]}) == 1,
            "Draft one reference at a time: use --reference-id; resolve unknown legacy references explicitly")
    require(len({unit["sample_id"] for unit in report["units"]}) == len(report["units"]),
            "Ambiguous duplicate biological sample IDs; correct upstream grouping before exporting")
    output.mkdir(parents=True, exist_ok=False)
    atomic_json(output / "input_inventory.json", report)
    (output / "READ_BEFORE_RUNNING.md").write_text(
        "# Draft upstream handoff\n\nThese files are not ready to run. Blank/null fields are unresolved evidence, not defaults.\n\n"
        "Review input_inventory.json and input_inventory.tsv. Complete subject/library/source identity, deduplication, "
        "source checksums, count units, filters and reference assets. Inspect stored-weight recipes and convergence. "
        "Choose comparisons explicitly; comparisons.tsv is intentionally header-only. Run normal scientific preflight "
        "after completing the drafts. Matrix inspection does not establish biological independence or formal inference.\n\n"
        + summary(report) + "\n", encoding="utf-8", newline="\n")
    samples, membership, assets, inventories = [], [], [], []
    for unit in report["units"]:
        sid = unit["sample_id"]
        samples.append({key: unit.get(key, "") for key in SAMPLE_FIELDS} | {
            "project": report["project"], "provenance_path": f"provenance/{sid}.json", "batch": unit["batch"]})
        balancing = {resolution: {key: info[key] for key in ("weight_column", "convention", "converged", "recipe")}
                     for resolution, info in unit["inspection"]["resolutions"].items() if info["available"]}
        # Only explicit upstream sheet fields are carried forward; no library/subject inference.
        for row in unit["sources"]:
            membership.append({"sample_id": sid, "upstream_sample": row["sample"], "source_id": row["sample"],
                               **{key: row.get(key, "") for key in ("library_id", "source_kind", "dedup_group_id", "source_checksum", "fastq_r1", "fastq_r2")},
                               "dedup_scope": row.get("dedup_scope") or "unknown"})
        scopes = {row.get("dedup_scope", "unknown") or "unknown" for row in unit["sources"]}
        atomic_json(output / "provenance" / f"{sid}.json", {
            "draft": True, "reference_id": unit["reference_id"], "count_units": None,
            "filters": {"mapq": None, "min_distance_bp": None, "chromosome_policy": None},
            "dedup_scope": next(iter(scopes)) if len(scopes) == 1 else "unknown",
            "library_ids": sorted({row.get("library_id") for row in unit["sources"] if row.get("library_id")}),
            "source_ids": unit["upstream_rows"], "upstream_commit": "", "balancing": balancing,
            "handoff_evidence": {"samplesheet": report["samplesheet"], "upstream_rows": unit["upstream_rows"],
                                 "matrix_path": unit["matrix_path"], "matrix_origin": unit["matrix_origin"],
                                 "individual_matrices": report.get("individual_matrices", False),
                                 "recorded_upstream_metadata": unit.get("recorded_upstream_metadata", []),
                                 "recorded_technical_membership": unit.get("recorded_technical_membership", []),
                                 "review_required": unit["review_required"]}})
        assets += [{"sample_id": sid, **asset} for asset in unit["related_outputs"]]
        inventories.append({"sample_id": sid, "condition": unit["condition"], "biological_replicate": unit["biological_replicate"],
                            "reference_id": unit["reference_id"], "upstream_rows": ";".join(unit["upstream_rows"]),
                            "matrix_path": unit["matrix_path"], "status": unit["status"], "issues": "; ".join(unit["errors"] + unit["review_required"])})
    with (output / "samples.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=(*SAMPLE_FIELDS, "batch"), lineterminator="\n")
        writer.writeheader()
        writer.writerows(samples)
    write_table(output / "technical_membership.tsv", membership, MEMBERSHIP_FIELDS)
    write_table(output / "comparisons.tsv", [], COMPARISON_FIELDS)
    write_table(output / "related_outputs.tsv", assets, ("sample_id", "owner", "category", "path"))
    write_table(output / "input_inventory.tsv", inventories, tuple(inventories[0]))
    atomic_json(output / "references.json", {"draft": True, "reference_id": report["units"][0]["reference_id"],
                                             "chromsizes": "", "chromsizes_sha256": "", "views": "", "views_sha256": "",
                                             "phasing_track": "", "phasing_track_sha256": ""})
    return output
