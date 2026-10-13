"""Explicit preparation of authentic NEW grids from retained exact-position pairs.

This utility is never called by a scientific stage. It trusts and records
explicit upstream deduplication/filter declarations; it does not deduplicate,
guess coordinate conventions, interpolate contacts, or overwrite source maps.
"""
from __future__ import annotations

import copy
import gzip
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

from .contracts import atomic_json, file_hash, load_structured, require, resolve, validate_metadata, write_table
from .resolutions import loop_resolutions, quantification_resolutions, chromosight_discovery_resolutions


def _header(path):
    opener = gzip.open if str(path).lower().endswith((".gz", ".bgz")) else open
    fields, assembly, chromosome_sizes = None, None, {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.startswith("#columns:"):
                require(fields is None, "Pairs header has duplicate #columns declarations")
                fields = line.split(":", 1)[1].split()
            elif line.startswith("#genome_assembly:"):
                assembly = line.split(":", 1)[1].strip()
            elif line.startswith("#chromsize:"):
                values = line.split(":", 1)[1].split()
                require(len(values) == 2 and values[1].isdigit() and values[0] not in chromosome_sizes,
                        "Pairs header has invalid/duplicate #chromsize declarations")
                chromosome_sizes[values[0]] = int(values[1])
    require(fields and len(fields) == len(set(fields)) and
            {"chrom1", "pos1", "chrom2", "pos2", "pair_type"} <= set(fields),
            f"{path}: require a unique #columns header containing chrom1,pos1,chrom2,pos2,pair_type")
    return fields, assembly, chromosome_sizes


def _pairs_pixels(path, columns, record, sizes, resolution, chunksize, stats):
    """Stream, aggregate within chunks, and leave global merging to Cooler."""
    counts = (sizes+resolution-1)//resolution
    offsets = counts.cumsum().shift(1, fill_value=0).to_dict()
    coordinate_shift = 1 if record["coordinate_system"] == "1-based" else 0
    allowed = set(record["upstream_filters"]["retained_pair_types"])
    minimum = record["upstream_filters"]["minimum_cis_distance_bp"]
    compression = "gzip" if str(path).lower().endswith((".gz", ".bgz")) else None
    iterator = pd.read_csv(path, sep="\t", comment="#", names=columns, dtype=str,
                           keep_default_na=False, na_filter=False, compression=compression, chunksize=chunksize)
    for frame in iterator:
        stats["records"] += len(frame)
        selected = frame.pair_type.isin(allowed)
        stats["excluded_pair_type"] += int((~selected).sum())
        frame = frame.loc[selected].copy()
        mapped = ~frame.chrom1.isin(("!", "*")) & ~frame.chrom2.isin(("!", "*"))
        stats["excluded_unmapped"] += int((~mapped).sum())
        frame = frame.loc[mapped].copy()
        if frame.empty:
            continue
        require(frame.chrom1.isin(sizes.index).all() and frame.chrom2.isin(sizes.index).all(),
                "Pairs contain chromosomes absent from registered chromosome sizes; correct assembly metadata")
        for anchor in (1, 2):
            column = f"pos{anchor}"
            require(frame[column].str.fullmatch(r"[0-9]+").all(), "Pairs positions must be exact nonnegative integers")
            try:
                frame[column] = frame[column].astype(np.int64)-coordinate_shift
            except (ValueError, OverflowError) as exc:
                raise ValueError("Pairs positions exceed supported exact integer coordinates") from exc
            lengths = frame[f"chrom{anchor}"].map(sizes)
            require(((frame[column] >= 0) & (frame[column] < lengths)).all(),
                    "Pairs position outside chromosome bounds; check declared coordinate_system")
        same = frame.chrom1 == frame.chrom2
        distance = abs(frame.pos2-frame.pos1)
        require((~same | (distance >= minimum)).all(),
                "Pairs contradict declared upstream minimum_cis_distance_bp; supply truthful retained-pairs provenance")
        i = (frame.chrom1.map(offsets)+frame.pos1//resolution).to_numpy(dtype=np.int64)
        j = (frame.chrom2.map(offsets)+frame.pos2//resolution).to_numpy(dtype=np.int64)
        pixels = pd.DataFrame({"bin1_id": np.minimum(i,j), "bin2_id": np.maximum(i,j), "count": np.ones(len(i), dtype=np.int64)})
        stats["retained"] += len(pixels)
        yield pixels.groupby(["bin1_id", "bin2_id"], sort=True, as_index=False)["count"].sum()


def _recipe(record):
    recipe = record.get("balancing_recipe", {})
    required = {"cis_only", "ignore_diags", "mad_max", "min_nnz", "min_count", "tol", "max_iters"}
    allowed = required | {"rescale_marginals"}
    require(isinstance(recipe, dict) and required <= set(recipe) <= allowed,
            "pairs_provenance.balancing_recipe requires cis_only,ignore_diags,mad_max,min_nnz,min_count,tol,max_iters; optional rescale_marginals")
    require(isinstance(recipe["cis_only"], bool) and all(isinstance(recipe[k], int) and not isinstance(recipe[k], bool) and recipe[k] >= 0
            for k in ("ignore_diags", "mad_max", "min_nnz", "min_count", "max_iters")) and recipe["max_iters"] > 0
            and isinstance(recipe["tol"], (int, float)) and np.isfinite(recipe["tol"]) and recipe["tol"] > 0,
            "Invalid explicit pairs balancing recipe")
    if "rescale_marginals" in recipe:
        require(isinstance(recipe["rescale_marginals"], bool), "rescale_marginals must be boolean")
    return dict(recipe, rescale_marginals=recipe.get("rescale_marginals", True))


def prepare_from_pairs(c, out, resolutions=None):
    """Return a registry/result dict after every requested NEW map converges.

    ``pairs_provenance`` in each sample's provenance JSON must declare a matching
    source SHA256, reference, complete=True, deduplicated=True, dedup_scope,
    coordinate_system (1-based or 0-based), explicit upstream_filters and an
    explicit balancing_recipe. Missing declarations fail before preparation.
    """
    root = Path(out).resolve()
    require(not root.exists(), "Pairs preparation destination exists; choose a NEW directory")
    grids = sorted(set(resolutions if resolutions is not None else
                       loop_resolutions(c)+quantification_resolutions(c)+
                       (chromosight_discovery_resolutions(c) if c.get("RUN_CHROMOSIGHT_DISCOVERY") else [])))
    require(grids and all(isinstance(r, int) and r > 0 for r in grids), "Pairs grids must be positive integers")
    metadata = validate_metadata(c)
    samples = metadata["samples"]
    require(len({s.get("pairs_path", "") for s in samples}) == len(samples),
            "Same retained-pairs file cannot supply two biological observations")
    reference = load_structured(c["REFERENCE_REGISTRY"])
    require(reference.get("reference_id") == c["REFERENCE_ID"], "Pairs preparation reference mismatch")
    sizes_path = resolve(reference.get("chromsizes", ""), Path(c["REFERENCE_REGISTRY"]).parent)
    require(sizes_path and Path(sizes_path).is_file(), "Registered chromosome sizes are required")
    sizes_hash = file_hash(sizes_path)
    require(reference.get("chromsizes_sha256") == sizes_hash, "Registered chromosome size SHA256 mismatch")
    sizes_frame = pd.read_csv(sizes_path, sep="\t", header=None, names=["chrom", "size"], dtype={"chrom": str, "size": np.int64})
    require(sizes_frame.chrom.is_unique and (sizes_frame["size"] > 0).all(), "Invalid chromosome sizes")
    sizes = sizes_frame.set_index("chrom")["size"]
    inputs, prepared, pair_hashes = {}, [], set()
    for sample in samples:
        path = sample.get("pairs_path", "")
        require(path and Path(path).is_file(), f"{sample['sample_id']}: retained pairs_path is required")
        provenance = load_structured(sample["provenance_path"])
        record = provenance.get("pairs_provenance", {})
        require(isinstance(record, dict) and record.get("reference_id") == c["REFERENCE_ID"]
                and provenance.get("reference_id") == c["REFERENCE_ID"], "Pairs provenance reference mismatch")
        require(record.get("complete") is True and record.get("deduplicated") is True and
                record.get("dedup_scope") in ("joint_library", "consolidated_biological_unit"),
                "Retained pairs need explicit completeness, deduplication and dedup_scope evidence")
        require(record.get("coordinate_system") in ("0-based", "1-based"), "Declare exact pairs coordinate_system (0-based or 1-based)")
        filters = record.get("upstream_filters", {})
        require(isinstance(filters, dict) and isinstance(filters.get("minimum_cis_distance_bp"), int)
                and not isinstance(filters["minimum_cis_distance_bp"], bool) and filters["minimum_cis_distance_bp"] >= 0
                and isinstance(filters.get("retained_pair_types"), list) and filters["retained_pair_types"]
                and all(isinstance(t, str) and re.fullmatch(r"[A-Z]{2}", t) for t in filters["retained_pair_types"]),
                "Declare upstream_filters minimum_cis_distance_bp and retained_pair_types; no filters are guessed")
        require(all(r["dedup_scope"] != "unknown" for r in metadata["membership"] if r["sample_id"] == sample["sample_id"]),
                "Unknown technical-library deduplication scope blocks pairs preparation")
        recipe = _recipe(record)
        columns, header_assembly, header_sizes = _header(path)
        require(header_assembly in (None, "", c["REFERENCE_ID"]), "Pairs header assembly contradicts registered reference")
        require(all(chrom in sizes and length == sizes[chrom] for chrom, length in header_sizes.items()),
                "Pairs header chromosome sizes contradict registered reference")
        source_hash = file_hash(path)
        require(record.get("source_sha256") == source_hash, "Retained pairs source_sha256 mismatch")
        require(source_hash not in pair_hashes, "Identical retained-pairs content cannot supply two biological observations")
        pair_hashes.add(source_hash)
        for input_path in (path, sample["provenance_path"], sample["matrix_path"].split("::")[0]):
            if Path(input_path).is_file():
                inputs[input_path] = file_hash(input_path)
        prepared.append((sample, provenance, record, recipe, columns))
    inputs[sizes_path] = sizes_hash
    inputs[c["REFERENCE_REGISTRY"]] = file_hash(c["REFERENCE_REGISTRY"])
    for key in ("SAMPLESHEET", "TECHNICAL_MEMBERSHIP"):
        inputs[c[key]] = file_hash(c[key])
    root.mkdir(parents=True)
    registry = {"reference_id": c["REFERENCE_ID"], "samples": {}}
    outputs, summaries, updated = [], [], []
    threads = min(int(c.get("THREADS_BALANCE", 4)), int(c.get("TOTAL_CPUS", 4)))
    budget = int(c.get("MEMORY_MB_BALANCE", 16000))*1000000
    require(threads > 0 and budget > 0, "Pairs preparation resources must be positive")
    # Bound string parsing and merge buffers conservatively within reservation.
    chunksize = min(int(c.get("PIXEL_CHUNKSIZE", 1000000)), max(1000, int(budget*.1)//512))
    try:
        for sample, original, record, recipe, columns in prepared:
            sid = sample["sample_id"]
            registry["samples"][sid] = {}
            provenance = copy.deepcopy(original)
            previous_stats = None
            for grid in grids:
                bins = cooler.binnify(sizes, grid)
                require(len(bins)*8*16+chunksize*256*threads <= budget*.8,
                        "Pairs bin vectors/balance chunks exceed MEMORY_MB_BALANCE; raise its reservation or lower PIXEL_CHUNKSIZE/THREADS_BALANCE")
                stats = {"records": 0, "retained": 0, "excluded_pair_type": 0, "excluded_unmapped": 0}
                target = root/f"{sid}.{grid}.cool"
                cooler.create_cooler(str(target), bins, _pairs_pixels(sample["pairs_path"], columns, record, sizes, grid, chunksize, stats),
                                     ordered=False, symmetric_upper=True, assembly=c["REFERENCE_ID"],
                                     dtypes={"count": np.int64}, mergebuf=chunksize, temp_dir=str(root), max_merge=20)
                clr = cooler.Cooler(str(target))
                require(stats["records"] == stats["retained"]+stats["excluded_pair_type"]+stats["excluded_unmapped"], "Pairs record accounting failed")
                require(int(clr.info["sum"]) == stats["retained"] and stats["retained"] > 0, "Pairs count conservation failed or no retained contacts")
                if record.get("records_total") is not None:
                    require(record["records_total"] == stats["records"], "Pairs completeness records_total mismatch")
                require(previous_stats in (None, stats), "Retained pairs changed between grids")
                previous_stats = dict(stats)
                with ThreadPoolExecutor(max_workers=threads) as executor:
                    weights, balancing = cooler.balance_cooler(clr, store=False, map=executor.map, chunksize=chunksize, **recipe)
                require(bool(np.asarray(balancing["converged"]).all()) and np.isfinite(weights).any(),
                        f"{sid} {grid} bp NEW balancing did not converge; no usable registry published")
                with h5py.File(target, "r+") as handle:
                    ds = handle["bins"].create_dataset("weight", data=weights)
                    ds.attrs.update(balancing)
                    ds.attrs.update(recipe)
                    ds.attrs["divisive_weights"] = False
                provenance.setdefault("balancing", {})[str(grid)] = {"weight_column": "weight", "convention": "multiplicative", "converged": True, "recipe": recipe}
                provenance.setdefault("derived_resolutions", {})[str(grid)] = {
                    "source_pairs": sample["pairs_path"], "source_pairs_sha256": record["source_sha256"],
                    "source_coordinate_system": record["coordinate_system"], "counts": stats,
                    "generated_by": f"cooler {cooler.__version__}", "method": "direct_exact_position_binning", "authentic_upsampling": False}
                registry["samples"][sid][str(grid)] = str(target)
                summaries.append({"sample_id": sid, "resolution": grid, "uri": str(target), **stats,
                                  "balancing_converged": True, "finite_weights": int(np.isfinite(weights).sum())})
                outputs.append(str(target))
            provpath = root/f"{sid}.provenance.json"
            provenance["preparation"] = {"source_provenance_sha256": inputs[sample["provenance_path"]], "inputs_modified": False, "new_maps_only": True}
            atomic_json(provpath, provenance)
            updated.append(dict(sample, matrix_path=registry["samples"][sid][str(grids[0])], provenance_path=str(provpath)))
            outputs.append(str(provpath))
        require(all(file_hash(path) == before for path, before in inputs.items()), "An upstream input changed during pairs preparation")
        registry_path = root/"matrix_registry.json"
        atomic_json(registry_path, registry)
        fields = list(dict.fromkeys(k for sample in updated for k in sample))
        write_table(root/"samples.tsv", updated, fields)
        atomic_json(root/"preparation.json", {"status": "prepared_not_scientifically_admitted", "reference_id": c["REFERENCE_ID"],
                    "resolutions": grids, "sample_count": len(samples), "rows": summaries, "inputs": inputs,
                    "inputs_modified": False, "deduplication": "upstream_declared_not_recomputed",
                    "next_action": "Run preflight and project-specific coverage/contact-retention/reproducibility review before fine-grid analysis."})
        outputs.extend(map(str, (registry_path, root/"samples.tsv", root/"preparation.json")))
        return {"matrix_registry": str(registry_path), "samplesheet": str(root/"samples.tsv"),
                "resolutions": grids, "sample_count": len(samples), "outputs": outputs,
                "scientific_admission": "not_performed"}
    except Exception as exc:
        atomic_json(root/"preparation.failed.json", {"status": "failed_no_usable_registry", "reason": str(exc), "inputs_modified": False})
        raise
