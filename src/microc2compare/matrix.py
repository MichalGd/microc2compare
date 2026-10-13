"""Read-only, chunked Cooler access and inherited-weight analysis masks."""
from __future__ import annotations

from pathlib import Path

import cooler
from cooler.core import RangeSelector1D
import numpy as np
import pandas as pd

from .contracts import (digest, file_hash, load_structured, read_json, require, resolve)
from .analyses import resolution_set
from .resolutions import sample_matrix
from .chromosomes import chromosome_selection


def matrix_uri(path, resolution):
    file, sep, group = str(path).partition("::")
    if sep:
        clr = cooler.Cooler(path)
        require(clr.binsize == resolution, f"{path}: requested resolution {resolution} unavailable")
        return str(path)
    groups = cooler.fileops.list_coolers(file)
    requested = f"/resolutions/{resolution}"
    if requested in groups:
        return file + "::" + requested
    if "/" in groups and cooler.Cooler(file).binsize == resolution:
        return file
    raise ValueError(f"{file}: missing {resolution}-bp resolution; request corrected upstream output")


def reference(c):
    ref = load_structured(c["REFERENCE_REGISTRY"])
    require(ref["reference_id"] == c["REFERENCE_ID"], "Reference registry identity mismatch")
    base = Path(c["REFERENCE_REGISTRY"]).parent
    for key in ("chromsizes", "views", "blacklist", "phasing_track"):
        if key == "phasing_track" and not c["RUN_COMPARTMENTS"]:
            continue
        if ref.get(key):
            ref[key] = resolve(ref[key], base)
            require(Path(ref[key]).is_file(), f"Missing reference {key}: {ref[key]}")
            require(ref.get(key + "_sha256") == file_hash(ref[key]), f"Reference fingerprint mismatch: {key}")
    require(ref.get("chromsizes") and ref.get("views"), "Reference requires chromsizes and views")
    sizes = pd.read_csv(ref["chromsizes"], sep="\t", header=None, names=["chrom", "length"])
    require(sizes["chrom"].is_unique and (sizes["length"] > 0).all(), "Invalid chromosome sizes")
    views = pd.read_csv(ref["views"], sep="\t", header=None, names=["chrom", "start", "end", "name"])
    require(views["name"].is_unique and (views["end"] > views["start"]).all(), "Invalid/duplicate view names")
    for row in views.itertuples():
        require(row.chrom in set(sizes.chrom) and 0 <= row.start < row.end <= int(sizes.set_index("chrom").loc[row.chrom, "length"]), "View outside reference bounds")
    for chrom, group in views.groupby("chrom", sort=False):
        ordered = group.sort_values("start")
        require(all(a <= b for a, b in zip(ordered.end.iloc[:-1], ordered.start.iloc[1:])), f"Overlapping views on {chrom}")
    selection = chromosome_selection(c, ref, sizes.chrom, views.chrom)
    ref["effective_scope"] = selection
    views = views.loc[views.chrom.isin(selection["chromosomes"])].reset_index(drop=True)
    return ref, sizes, views


def weights(clr, c, provenance, resolution, needed=True):
    column = c["WEIGHT_COLUMN"]
    bins = clr.bins()[:]
    if not needed:
        return np.ones(len(bins), dtype=float)
    require(column in bins, f"{clr.uri}: missing {column}; request upstream balancing, no downstream repair")
    values = bins[column].to_numpy(dtype=float, copy=True)
    metadata = provenance.get("balancing", {}).get(str(resolution), {})
    require(metadata.get("weight_column", column) == column, "Weight-column provenance mismatch")
    require(metadata.get("convention") == c["WEIGHT_CONVENTION"], "Weight convention must be declared and compatible")
    if metadata.get("converged") is not True:
        require(metadata.get("converged") is None and c["LEGACY_WEIGHT_VALIDATION"] == "finite_positive",
                f"{clr.uri}: missing/failed convergence evidence; corrected upstream output or explicit legacy validation needed")
    require(not np.isinf(values).any() and not (values[np.isfinite(values)] <= 0).any(), "Infinite/nonpositive weights are invalid")
    require(np.isfinite(values).any(), f"{clr.uri}: no usable weights")
    if c["WEIGHT_CONVENTION"] == "divisive":
        values = 1.0 / values
    return values


def view_labels(bins, views):
    labels = np.full(len(bins), -1, dtype=int)
    for i, view in enumerate(views.itertuples()):
        mask = (bins.chrom.to_numpy().astype(str) == view.chrom) & (bins.start.to_numpy() >= view.start) & (bins.end.to_numpy() <= view.end)
        labels[mask] = i
    return labels


def analysis_bins(sizes, views, resolution):
    """Full bin grids on explicitly requested view chromosomes, in reference order."""
    selected = sizes.loc[sizes.chrom.isin(views.chrom)].set_index("chrom")["length"]
    return cooler.binnify(selected, resolution).reset_index(drop=True)


def tool_views(clr, views):
    """cooltools requires views in each file's own chromosome/coordinate order."""
    order = sorted(range(len(views)), key=lambda i: clr.extent(tuple(views.iloc[i][["chrom", "start", "end"]]))[0])
    return views.iloc[order].reset_index(drop=True)


def native_indices(clr, bins):
    """Map registered coordinates to native IDs; never assume file chromosome order."""
    native = clr.bins()[:][["chrom", "start", "end"]].copy()
    native["chrom"] = native.chrom.astype(str)
    wanted = bins.copy()
    wanted["chrom"] = wanted.chrom.astype(str)
    for chrom, length in wanted.groupby("chrom", sort=False).end.max().items():
        require(chrom in clr.chromsizes and int(clr.chromsizes[chrom]) == int(length),
                f"{clr.uri}: selected chromosome {chrom} is missing or its length differs from reference")
    require(not native.duplicated().any(), f"{clr.uri}: duplicate bin coordinates")
    # Checking the complete selected chromosome grid also checks its terminal length.
    selected = native.loc[native.chrom.isin(wanted.chrom)]
    require(len(selected) == len(wanted), f"{clr.uri}: selected chromosome bin grids differ")
    index = pd.MultiIndex.from_frame(native)
    ids = index.get_indexer(pd.MultiIndex.from_frame(wanted))
    require((ids >= 0).all(), f"{clr.uri}: missing chromosome or incompatible selected bin coordinates/lengths")
    return ids


def native_domain(clr, bins, common, labels, values):
    """Expand the common coordinate mask/weights onto one file's native axis."""
    ids = native_indices(clr, bins)
    native = clr.bins()[:][["chrom", "start", "end"]].copy()
    native["chrom"] = native.chrom.astype(str)
    mask = np.zeros(len(native), dtype=bool)
    native_labels = np.full(len(native), -1, dtype=int)
    native_values = np.full(len(native), np.nan)
    mask[ids] = common
    native_labels[ids] = labels
    native_values[ids] = values
    native_values[~mask] = np.nan
    return native, mask, native_labels, native_values


def domain(c, samples, resolution):
    ref, sizes, views = reference(c)
    bins = analysis_bins(sizes, views, resolution)
    labels = view_labels(bins, views)
    for view in views.itertuples():
        require(view.start%resolution==0 and (view.end%resolution==0 or view.end==int(sizes.set_index("chrom").loc[view.chrom,"length"])),
                f"View {view.name} is not aligned to requested {resolution}-bp bins")
    common = labels >= 0
    recipes, sample_weights = [], {}
    from .analyses import needs_loop_counts
    balanced_needed = c["RUN_COMPARTMENTS"] or c["RUN_INSULATION"] or c["RUN_LOOPS"] or c["RUN_AGGREGATES"] or c["RUN_QC"] or c.get("RUN_STRIPES", False) or needs_loop_counts(c)
    for s in samples:
        clr = cooler.Cooler(matrix_uri(sample_matrix(c, s, resolution), resolution))
        require(clr.binsize == resolution, "Variable-width bins are unsupported")
        ids = native_indices(clr, bins)
        require(clr.info.get("storage-mode") in ("symmetric-upper", "square"), "Unsupported storage mode")
        require(np.issubdtype(clr.pixels()["count"][0:0].dtype, np.integer), "Count field must store raw integers")
        prov = read_json(s["provenance_path"])
        require(prov.get("reference_id") == c["REFERENCE_ID"] and prov.get("count_units") == "unique_contacts", "Reference/count semantics require explicit compatible provenance")
        values = weights(clr, c, prov, resolution, needed=balanced_needed)[ids]
        sample_weights[s["sample_id"]] = values
        common &= np.isfinite(values)
        require(prov.get("filters"), "Missing upstream filters")
        # Actual chromosome scope is audited separately. All contact filters and
        # numerical balancing settings must still agree; no missing data is padded.
        filters = {k: v for k, v in prov["filters"].items() if k != "chromosome_policy"}
        recipes.append({"filters": filters, "balancing": prov.get("balancing", {}).get(str(resolution), {}).get("recipe", {})})
    require(len({digest(recipe) for recipe in recipes}) == 1, "Incompatible upstream filter/balancing recipes; correct upstream outputs")
    if ref.get("blacklist"):
        blacklist = pd.read_csv(ref["blacklist"], sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
        for b in blacklist.itertuples():
            common &= ~((bins.chrom.astype(str).to_numpy() == b.chrom) & (bins.start.to_numpy() < b.end) & (bins.end.to_numpy() > b.start))
    for key in sample_weights:
        sample_weights[key][~common] = np.nan
    require(common.any(), "No common usable bins remain in the registered analysis views")
    return bins, common, labels, views, sample_weights


def scope_record(c, samples, resolution, bins, common):
    selected = list(dict.fromkeys(bins.chrom.astype(str)))
    entries = []
    for sample in samples:
        clr = cooler.Cooler(matrix_uri(sample_matrix(c, sample, resolution), resolution))
        prov = read_json(sample["provenance_path"])
        entries.append({"sample_id": sample["sample_id"],
                        "native_chromosomes": list(clr.chromnames),
                        "excluded_chromosomes": [x for x in clr.chromnames if x not in selected],
                        "upstream_chromosome_policy": prov["filters"].get("chromosome_policy"),
                        "balancing_recipe": prov.get("balancing", {}).get(str(resolution), {}).get("recipe", {})})
    different = len({digest(e["native_chromosomes"]) for e in entries}) > 1
    return {"resolution": resolution, "analysis_chromosomes": selected,
            "chromosome_selection": reference(c)[0]["effective_scope"],
            "common_eligible_bins": int(common.sum()), "analysis_grid_bins": len(bins),
            "chromosome_coverage": [{"chrom": chrom, "grid_bins": int((bins.chrom.astype(str) == chrom).sum()),
                                    "common_eligible_bins": int((common & (bins.chrom.astype(str).to_numpy() == chrom)).sum())}
                                   for chrom in selected],
            "native_chromosome_sets_or_orders_differ": different, "samples": entries,
            "weights": "inherited_without_rebalancing",
            "interpretation": "Common coordinates and masks do not equalize upstream balancing scopes. Assess scope and sex-chromosome coverage in the real-data pilot."}


class BinSelector(RangeSelector1D):
    """Cooler-compatible lazy bin queries with analysis weights in memory only."""
    def __init__(self, clr, values, fields=None, **kwargs):
        self.clr, self.values, self.bin_kwargs = clr, values, kwargs

        def slicer(selected_fields, lo, hi):
            frame = clr.bins(**kwargs)[lo:hi].copy()
            frame["weight"] = values[frame.index.to_numpy()]
            return frame if selected_fields is None else frame[selected_fields]

        # cooler.annotate checks this selector type to choose end-exclusive
        # row slicing. Retain Cooler semantics without loading a whole bin table.
        super().__init__(fields, slicer, clr.extent, clr.info["nbins"])

    def __getitem__(self, key):
        if isinstance(key, (str, list)):
            return BinSelector(self.clr, self.values, key, **self.bin_kwargs)
        return super().__getitem__(key)


class MaskedCooler(cooler.Cooler):
    """Virtual canonical multiplicative weights; never writes a bin column."""
    def __init__(self, uri, values):
        super().__init__(uri)
        self.analysis_weights = np.asarray(values, dtype=float)
        self.raw_cooler = cooler.Cooler(uri)

    def bins(self, **kwargs):
        return BinSelector(self.raw_cooler, self.analysis_weights, **kwargs)

    def matrix(self, field=None, balance=True, sparse=False, as_pixels=False, **kwargs):
        if not balance:
            return self.raw_cooler.matrix(field=field, balance=False, sparse=sparse, as_pixels=as_pixels, **kwargs)
        require(not as_pixels, "Balanced pixel-selector adapter is unsupported; use bounded matrix queries")
        require(balance is True or balance == "weight", "Adapter exposes canonical multiplicative weight only")
        return WeightedSelector(self.raw_cooler, self.analysis_weights, field, sparse)


class WeightedSelector:
    def __init__(self, clr, values, field, sparse):
        self.clr, self.values, self.field, self.sparse = clr, values, field, sparse

    def __getitem__(self, key):
        from scipy.sparse import diags
        rows, cols = key
        array = self.clr.matrix(field=self.field, balance=False, sparse=self.sparse)[rows, cols]
        i = np.arange(len(self.values))[rows]
        j = np.arange(len(self.values))[cols]
        if self.sparse:
            return diags(self.values[i]) @ array @ diags(self.values[j])
        return array * self.values[i, None] * self.values[None, j]

    def fetch(self, region1, region2=None):
        lo, hi = self.clr.extent(region1)
        lo2, hi2 = self.clr.extent(region1 if region2 is None else region2)
        return self[lo:hi, lo2:hi2]


def pixel_chunks(clr, chunksize):
    for start in range(0, clr.info["nnz"], chunksize):
        frame = clr.pixels()[start:start+chunksize]
        require((frame["count"] >= 0).all(), f"Negative raw count in {clr.uri}")
        # Count one representative of unordered pixels under either storage mode.
        yield frame.loc[frame.bin1_id <= frame.bin2_id]


def aligned_pixel_chunks(clr, bins, chunksize):
    """Chunked raw counts on the selected coordinate axis; output may be unordered.

    Used only for explicitly requested NEW pools/thinned maps. Consumers must
    use Cooler's external merge/sort rather than assuming native index order.
    """
    ids = native_indices(clr, bins)
    mapping = np.full(clr.info["nbins"], -1, dtype=np.int64)
    mapping[ids] = np.arange(len(bins))
    for frame in pixel_chunks(clr, chunksize):
        i = mapping[frame.bin1_id.to_numpy()]
        j = mapping[frame.bin2_id.to_numpy()]
        keep = (i >= 0) & (j >= 0)
        yield pd.DataFrame({"bin1_id": np.minimum(i[keep], j[keep]),
                            "bin2_id": np.maximum(i[keep], j[keep]),
                            "count": frame["count"].to_numpy()[keep]}).sort_values(["bin1_id", "bin2_id"])


def unordered_total(clr,chunksize):
    if clr.info.get("storage-mode")=="symmetric-upper":return int(clr.info["sum"])
    return sum(int(chunk["count"].sum()) for chunk in pixel_chunks(clr,chunksize))


def raw_summary(clr, bins, mask, labels, c):
    total = cis = trans = exposure = selected_cis = distance_cis = selected_total = 0
    chrom = bins.chrom.astype(str).to_numpy()
    start = bins.start.to_numpy()
    for frame in pixel_chunks(clr, c["PIXEL_CHUNKSIZE"]):
        i, j, count = [frame[k].to_numpy() for k in ("bin1_id", "bin2_id", "count")]
        total += int(count.sum())
        same = chrom[i] == chrom[j]
        cis += int(count[same].sum()); trans += int(count[~same].sum())
        distance = start[j] - start[i]
        selected = (labels[i] >= 0) & (labels[j] >= 0)
        same_view = selected & (labels[i] == labels[j])
        distance_ok = same_view & (distance >= c["LOOP_MIN_DISTANCE_BP"]) & (distance <= c["LOOP_MAX_DISTANCE_BP"])
        selected_total += int(count[selected].sum())
        selected_cis += int(count[same_view].sum())
        distance_cis += int(count[distance_ok].sum())
        eligible = mask[i] & mask[j] & (labels[i] >= 0) & (labels[i] == labels[j]) & (distance >= c["LOOP_MIN_DISTANCE_BP"]) & (distance <= c["LOOP_MAX_DISTANCE_BP"])
        exposure += int(count[eligible].sum())
    return {"raw_total": total, "raw_cis": cis, "raw_trans": trans, "eligible_cis_exposure": exposure,
            "selected_view_total": selected_total, "selected_view_cis": selected_cis,
            "distance_eligible_cis_before_mask": distance_cis,
            "mask_excluded_distance_contacts": distance_cis-exposure,
            "excluded_bin_fraction": float(1-mask[labels >= 0].mean()) if (labels >= 0).any() else 1.0,
            "excluded_native_bin_fraction": float(1-mask.mean()),
            "raw_total_scope": "entire_native_matrix", "eligible_exposure_scope": "common_registered_views"}


def expected_table(clr, bins, mask, labels, views, values, c):
    """Streaming sums + FFT eligible-pair denominators; includes observable zeros."""
    from scipy.signal import fftconvolve
    frames, accumulators = [], {}
    for vid, view in enumerate(views.itertuples()):
        ids = np.flatnonzero(labels == vid)
        if not len(ids): continue
        valid = mask[ids].astype(float)
        denom = np.rint(fftconvolve(valid, valid[::-1], mode="full")[len(valid)-1:]).astype(np.int64)
        accumulators[vid] = {"offset": ids[0], "denom": denom, "raw": np.zeros(len(ids)), "balanced": np.zeros(len(ids))}
    for pixels in pixel_chunks(clr, c["PIXEL_CHUNKSIZE"]):
        i, j, count = [pixels[k].to_numpy() for k in ("bin1_id", "bin2_id", "count")]
        valid = mask[i] & mask[j] & (labels[i] >= 0) & (labels[i] == labels[j])
        for vid in np.unique(labels[i[valid]]):
            sel = valid & (labels[i] == vid)
            d = j[sel]-i[sel]
            a = accumulators[vid]
            a["raw"] += np.bincount(d, weights=count[sel], minlength=len(a["denom"]))
            a["balanced"] += np.bincount(d, weights=count[sel]*values[i[sel]]*values[j[sel]], minlength=len(a["denom"]))
    cutoff = c["INSULATION_IGNORE_DIAGS"]
    for vid, a in accumulators.items():
        name = views.iloc[vid]["name"]
        distance = np.arange(len(a["denom"]), dtype=np.int64)
        usable = (a["denom"] > 0) & (distance >= cutoff)
        frames.append(pd.DataFrame({
            "region1": name, "region2": name, "dist": distance,
            "distance_bp": distance*clr.binsize, "n_valid": a["denom"],
            "count.sum": a["raw"], "balanced.sum": a["balanced"],
            "count.avg": np.divide(a["raw"], a["denom"], out=np.full(len(distance), np.nan), where=usable),
            "balanced.avg": np.divide(a["balanced"], a["denom"], out=np.full(len(distance), np.nan), where=usable),
        }))
    columns = ["region1", "region2", "dist", "distance_bp", "n_valid", "count.sum", "balanced.sum", "count.avg", "balanced.avg"]
    return pd.concat(frames, ignore_index=True)[columns] if frames else pd.DataFrame(columns=columns)


def bounded_raw(clr, region1, region2, c):
    lo, hi = clr.extent(region1); lo2, hi2 = clr.extent(region2)
    require(max(hi-lo, hi2-lo2) <= c["MAX_DENSE_BINS"], "Requested dense region exceeds MAX_DENSE_BINS")
    return clr.matrix(balance=False).fetch(region1, region2)
