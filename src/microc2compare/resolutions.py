"""Resolution policy and optional per-sample matrix URI overrides."""
from pathlib import Path
from .contracts import load_structured, require, resolve


def loop_resolutions(c):
    requested = c.get("LOOP_DISCOVERY_RESOLUTIONS", "")
    result = [int(x) for x in requested.split(",") if x] if requested else [c["LOOP_RESOLUTION"]]
    sensitivity = c.get("LOOP_SENSITIVITY_RESOLUTIONS", "")
    result += [int(x) for x in sensitivity.split(",") if x]
    return sorted(set(result))


def quantification_resolutions(c):
    value = c.get("LOOP_QUANTIFICATION_RESOLUTIONS", "")
    return sorted({int(x) for x in value.split(",") if x} or {c["LOOP_RESOLUTION"]})


def chromosight_discovery_resolutions(c):
    value = c.get("CHROMOSIGHT_DISCOVERY_RESOLUTIONS", "")
    return sorted({int(x) for x in value.split(",") if x}) if value else loop_resolutions(c)


def sample_matrix(c, sample, resolution):
    registry_path = c.get("MATRIX_REGISTRY", "")
    if not registry_path:
        return sample["matrix_path"]
    registry = load_structured(registry_path)
    require(registry.get("reference_id") == c["REFERENCE_ID"], "Matrix registry reference mismatch")
    entry = registry.get("samples", {}).get(sample["sample_id"], {}).get(str(resolution))
    if not entry:
        return sample["matrix_path"]
    require(isinstance(entry, str), "Matrix registry values must be matrix paths/URIs")
    return resolve(entry, Path(registry_path).parent)


def input_capabilities(c, samples, requested=None):
    """Inspect native grid metadata without reading pixels or modifying inputs.

    This discovery report is deliberately separate from scientific admission:
    inherited weights, depth and reproducibility still require validation.
    Explicit URI registrations are honored; other native grids can be proposed
    as preparation sources but are not silently substituted for active inputs.
    """
    import cooler
    import numpy as np
    from .matrix import matrix_uri

    grids = sorted(set(requested if requested is not None else
                       loop_resolutions(c)+quantification_resolutions(c)+
                       (chromosight_discovery_resolutions(c) if c.get("RUN_CHROMOSIGHT_DISCOVERY") else [])))
    require(grids and all(isinstance(r, int) and r > 0 for r in grids), "Requested grids must be positive integers")
    registry_path = c.get("MATRIX_REGISTRY", "")
    registry = load_structured(registry_path) if registry_path else {"samples": {}}
    if registry_path:
        require(registry.get("reference_id") == c["REFERENCE_ID"], "Matrix registry reference mismatch")
    rows, inventories = [], {}
    for sample in samples:
        sid = sample["sample_id"]
        require(sample.get("reference_id", c["REFERENCE_ID"]) == c["REFERENCE_ID"], f"{sid}: reference mismatch")
        paths = [sample["matrix_path"]]
        for entry in registry.get("samples", {}).get(sid, {}).values():
            require(isinstance(entry, str), "Matrix registry values must be matrix paths/URIs")
            paths.append(resolve(entry, Path(registry_path).parent))
        inventory = []
        seen = set()
        for path in dict.fromkeys(paths):
            filename, sep, group = str(path).partition("::")
            if not Path(filename).is_file():
                inventory.append({"uri": str(path), "status": "missing_file"})
                continue
            groups = [group] if sep else cooler.fileops.list_coolers(filename)
            for native_group in groups:
                uri = filename if native_group == "/" else filename+"::"+native_group
                if uri in seen:
                    continue
                seen.add(uri)
                try:
                    clr = cooler.Cooler(uri)
                    with clr.open("r") as handle:
                        count_dtype = handle["pixels/count"].dtype
                        weight = c.get("WEIGHT_COLUMN", "weight")
                        weight_present = "bins/"+weight in handle
                        attrs = dict(handle["bins/"+weight].attrs) if weight_present else {}
                    assembly = clr.info.get("genome-assembly")
                    inventory.append({"uri": uri, "resolution": clr.binsize, "assembly": assembly,
                                      "status": "native_metadata", "integer_counts": bool(np.issubdtype(count_dtype, np.integer)),
                                      "weight_column_present": weight_present,
                                      "weight_converged": bool(np.all(attrs.get("converged", False))),
                                      "assembly_matches": assembly in (None, "", c["REFERENCE_ID"])})
                except (OSError, KeyError, ValueError) as exc:
                    inventory.append({"uri": uri, "status": "unreadable", "reason": str(exc)})
        inventories[sid] = inventory
        native = [r for r in inventory if r.get("resolution") and r.get("assembly_matches") and r.get("integer_counts")]
        pairs = sample.get("pairs_path", "")
        pairs_record = {}
        if sample.get("provenance_path") and Path(sample["provenance_path"]).is_file():
            pairs_record = load_structured(sample["provenance_path"]).get("pairs_provenance", {})
        retained_minimum = pairs_record.get("upstream_filters", {}).get("minimum_cis_distance_bp")
        for grid in grids:
            active_uri, active_error = None, None
            try:
                active_uri = matrix_uri(sample_matrix(c, sample, grid), grid)
            except (OSError, KeyError, ValueError) as exc:
                active_error = str(exc)
            source = max((r for r in native if grid % r["resolution"] == 0),
                         key=lambda r: r["resolution"], default=None)
            if active_uri:
                detail = next((r for r in native if r["uri"] == active_uri), None)
                status = "available" if detail else "invalid_native_metadata"
                action = "Validate inherited weight conventions, assembly, common masks and scientific resolution eligibility."
            elif source:
                status = "coarsenable"
                action = "Explicitly prepare a NEW derived grid with conserved integer counts, a recorded balancing recipe and reviewed provenance; register its URI."
            else:
                status = "retained_pairs_required" if pairs and Path(pairs).is_file() else "missing_native_grid"
                action = ("Supply retained deduplicated exact-position pairs (plus assembly, chromosome sizes, filtering/dedup statistics and library membership), or an authentic finer grid that divides the target; do not split coarse counts or copy coarse weights.")
            rows.append({"sample_id": sid, "reference_id": c["REFERENCE_ID"], "resolution": grid,
                         "status": status, "active_uri": active_uri, "active_error": active_error,
                         "coarsening_source_uri": source["uri"] if source else None,
                         "pairs_path": pairs or None, "pairs_present": bool(pairs and Path(pairs).is_file()),
                         "pairs_declared_minimum_cis_distance_bp": retained_minimum,
                         "pairs_declared_deduplicated": pairs_record.get("deduplicated"),
                         "pairs_completeness_declared": pairs_record.get("complete"),
                         "pairs_coordinate_system": pairs_record.get("coordinate_system"),
                         "short_range_retention_note": f"Contacts below {retained_minimum} bp were removed upstream and cannot be recovered by finer binning." if retained_minimum else None,
                         "fine_grid_review_required": grid <= 1000, "action": action})
    return {"reference_id": c["REFERENCE_ID"], "requested_resolutions": grids,
            "sample_count": len(samples), "native_inputs": inventories, "rows": rows,
            "scientific_admission": "not_performed", "input_mutation": False}
