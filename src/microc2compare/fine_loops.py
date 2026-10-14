"""Auditable multi-resolution candidates and read-only fine-loop tool adapters."""
from __future__ import annotations
import importlib.metadata
import math
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from .contracts import atomic_json, digest, file_hash, require
from .features import LOOP_COLUMNS, canonical_loops, feature_id
from .loop_profiles import mustache_resource_plan


def _mustache_sparse_band(c, clr, lo, hi, mask, values, distance, resolution):
    """Read each upper-band pixel once, without fetching a full chromosome."""
    parts, known_bytes = [], 0
    for start in range(lo, hi, 512):
        end = min(start+512, hi)
        raw = clr.matrix(balance=False, sparse=True)[start:end, start:min(end+distance, hi)].tocoo()
        x, y = raw.row+start, raw.col+start
        valid = (y >= x) & (y-x <= distance) & mask[x] & mask[y]
        x, y = x[valid].astype(np.int64), y[valid].astype(np.int64)
        v = raw.data[valid].astype(float)*values[x]*values[y]
        known_bytes += x.nbytes+y.nbytes+v.nbytes
        # Concatenation briefly duplicates accumulated arrays; reject before
        # growing the global band beyond the configured reservation.
        mustache_resource_plan(c, resolution, 2*known_bytes)
        parts.append((x-lo, y-lo, v))
    if not parts:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([], dtype=float)
    return tuple(np.concatenate([p[k] for p in parts]) for k in range(3))


def cross_resolution_catalogue(c, sources):
    records = []
    for source, resolution, frame in sources:
        for row in canonical_loops(frame, c["REFERENCE_ID"], resolution).to_dict("records"):
            row.update(source=source, resolution=resolution)
            row["caller"] = row.get("caller", "unspecified")
            row["source_call_id"] = digest([source, resolution, row["caller"], row["feature_id"]])[:20]
            records.append(row)
    records.sort(key=lambda r: (r["chrom1"], r["chrom2"], r["start1"], r["start2"], r["resolution"], r["source"], r["caller"]))
    tolerance = c["CROSS_RESOLUTION_TOLERANCE_BINS"]
    groups, active = [], []
    maximum = max([r["resolution"] for r in records] or [1]) * tolerance
    def agrees(a, b):
        bound = max(a["resolution"], b["resolution"]) * tolerance
        return a["chrom1"] == b["chrom1"] and a["chrom2"] == b["chrom2"] and all(
            abs((a[f"start{k}"]+a[f"end{k}"])-(b[f"start{k}"]+b[f"end{k}"])) <= 2*bound for k in (1,2))
    for row in records:
        active = [g for g in active if g[0]["chrom1"] == row["chrom1"] and g[0]["chrom2"] == row["chrom2"] and row["start1"]-max(x["end1"] for x in g) <= maximum]
        matches = [g for g in active if all(agrees(row, other) for other in g)]
        if matches:
            # Stable first complete-linkage group; report ambiguities instead of transitive union.
            row["matching_cluster_count"] = len(matches)
            matches[0].append(row)
        else:
            row["matching_cluster_count"] = 1
            group = [row]; groups.append(group); active.append(group)
    candidates, members = [], []
    for group in groups:
        finest = min(r["resolution"] for r in group)
        eligible = [r for r in group if r["resolution"] == finest]
        rep = min(eligible, key=lambda a: (sum(abs(a["start1"]-b["start1"])+abs(a["start2"]-b["start2"]) for b in group), a["start1"], a["start2"], a["source"]))
        identity = [c["REFERENCE_ID"], *[rep[k] for k in LOOP_COLUMNS[:6]], finest]
        cid = "candidate_" + digest(identity)[:20]
        resolutions = sorted({r["resolution"] for r in group})
        candidates.append({**{k: rep[k] for k in LOOP_COLUMNS[:6]}, "feature_id": cid,
                           "representative_resolution": finest, "supporting_resolutions": ";".join(map(str,resolutions)),
                           "source_count": len(group), "biological_support": len({r["source"] for r in group if not r["source"].startswith("pool__")}),
                           "supporting_callers": ";".join(sorted({r["caller"] for r in group})),
                           "resolution_status": "multiple_resolutions" if len(resolutions)>1 else "single_resolution",
                           "ambiguous_correspondence": any(r["matching_cluster_count"]>1 for r in group),
                           "inference_status": "descriptive_same_data"})
        members.extend(dict(r, candidate_id=cid) for r in group)
    return (pd.DataFrame(candidates, columns=LOOP_COLUMNS+["representative_resolution","supporting_resolutions","source_count","biological_support","supporting_callers","resolution_status","ambiguous_correspondence","inference_status"]),
            pd.DataFrame(members, columns=list(dict.fromkeys(LOOP_COLUMNS+["source","resolution","source_call_id","candidate_id","matching_cluster_count"]+[k for r in records for k in r]))))


def map_candidates(candidates, reference_id, resolution):
    """One measurement per pixel; every original candidate retains its alias."""
    pixels, aliases = {}, []
    for row in candidates.to_dict("records"):
        a = (row["chrom1"], ((int(row["start1"])+int(row["end1"])-1)//2)//resolution*resolution)
        b = (row["chrom2"], ((int(row["start2"])+int(row["end2"])-1)//2)//resolution*resolution)
        if b<a: a,b=b,a
        coordinates = (a[0],a[1],a[1]+resolution,b[0],b[1],b[1]+resolution)
        mid = feature_id("loop", reference_id, coordinates, resolution)
        aliases.append({"candidate_id": row["feature_id"], "measurement_id": mid, "resolution": resolution,
                        "original_resolution": row.get("representative_resolution", row.get("resolution", resolution)),
                        "original_start1": row["start1"], "original_start2": row["start2"],
                        "original_end1": row["end1"], "original_end2": row["end2"],
                        "original_distance_bp": int(row["start2"])-int(row["start1"]) if row["chrom1"]==row["chrom2"] else np.nan,
                        "measured_distance_bp": b[1]-a[1] if a[0]==b[0] else np.nan,
                        "same_bin": a==b, "mapping_policy":"anchor_midpoint_lower_bin_boundary_tie",
                        "coarse_to_fine_projection":resolution < int(row.get("representative_resolution", row.get("resolution", resolution)))})
        pixels[mid] = dict(zip(LOOP_COLUMNS, (*coordinates, mid)))
    counts = pd.Series([r["measurement_id"] for r in aliases]).value_counts()
    for row in aliases: row["collapsed_candidate_count"] = int(counts[row["measurement_id"]])
    return pd.DataFrame(list(pixels.values()),columns=LOOP_COLUMNS), pd.DataFrame(aliases,columns=["candidate_id","measurement_id","resolution","original_resolution","original_start1","original_end1","original_start2","original_end2","original_distance_bp","measured_distance_bp","same_bin","collapsed_candidate_count","mapping_policy","coarse_to_fine_projection"])


def mustache_calls(c, clr, bins, mask, labels, views, values, resolution):
    from mustache import mustache as engine
    require(importlib.metadata.version("mustache-hic") == "1.3.3", "Mustache adapter requires tested mustache-hic 1.3.3")
    require(c["MUSTACHE_ITERATIONS"] == 10, "Mustache 1.3.3 hardcodes ten scale iterations; other values are unsupported")
    plan = mustache_resource_plan(c, resolution)
    distance, chunk_size = plan["maximum_distance_bins"], plan["tile_bins"]
    outputs, resource_plans = [], []
    octaves = [c["MUSTACHE_SIGMA_ZERO"]*2**i for i in range(c["MUSTACHE_OCTAVES"])]
    for vid, view in enumerate(views.itertuples()):
        lo, hi = clr.extent((view.chrom,int(view.start),int(view.end)))
        x, y, v = _mustache_sparse_band(c, clr, lo, hi, mask, values, distance, resolution)
        if not len(v): continue
        # Identical library distance transform and caller; no iterative matrix balancing.
        pweights = engine.normalize_sparse(x,y,v,resolution,distance)
        n = hi-lo
        starts = list(range(0,max(n-chunk_size,0)+1,max(chunk_size-distance,1)))
        if not starts or starts[-1]+chunk_size<n: starts.append(max(0,n-chunk_size))
        starts = sorted(set(starts))
        def call(start):
            end=min(start+chunk_size,n)
            selected=(x>=start)&(x<end)&(y>=start)&(y<end)
            tile=np.zeros((end-start,end-start),dtype=float)
            tile[x[selected]-start,y[selected]-start]=v[selected]
            return engine.mustache(tile,view.chrom,view.chrom,resolution,pweights,start,end,-1,distance,octaves,c["MUSTACHE_SPARSITY_THRESHOLD"],c["MUSTACHE_P_THRESHOLD"])
        plan = mustache_resource_plan(c, resolution, x.nbytes+y.nbytes+v.nbytes+np.asarray(pweights).nbytes)
        resource_plans.append(dict(plan, view=str(view.name), tile_count=len(starts)))
        with ThreadPoolExecutor(max_workers=plan["workers"]) as workers:
            for result in workers.map(call,starts):
                for i,j,qvalue,scale in result:
                    i,j=int(i)+lo,int(j)+lo
                    if not (mask[i] and mask[j] and labels[i]==labels[j]==vid): continue
                    a,b=bins.iloc[i],bins.iloc[j]
                    if a.end-a.start != resolution or b.end-b.start != resolution: continue
                    if not c["LOOP_MIN_DISTANCE_BP"]<=b.start-a.start<=c["LOOP_MAX_DISTANCE_BP"]: continue
                    outputs.append(dict(chrom1=str(a.chrom),start1=int(a.start),end1=int(a.end),chrom2=str(b.chrom),start2=int(b.start),end2=int(b.end),caller_qvalue=float(qvalue),caller_scale_bins=float(scale)))
    frame=pd.DataFrame(outputs,columns=LOOP_COLUMNS[:6]+["caller_qvalue","caller_scale_bins"])
    if len(frame): frame=frame.sort_values("caller_qvalue").drop_duplicates(LOOP_COLUMNS[:6])
    result = canonical_loops(frame,c["REFERENCE_ID"],resolution)
    result.attrs["resource_plans"] = resource_plans
    return result


def _expected_for(c, scope, resolution, clr, bins, mask, labels, views, values):
    """Reuse only a hash-verified completed QC output for this run and mask."""
    from .science import expected_table, expected_lookup
    from .state import complete, record_path
    from .contracts import read_json
    path = Path(c["run_dir"])/f"01_qc/ps/{scope}.{resolution}.tsv"
    if path.is_file() and complete(c, "qc", scope):
        record = read_json(record_path(c, "qc", scope))
        relative = str(path.relative_to(c["run_dir"]))
        if any(r["path"] == relative and not r.get("discarded_intermediate") for r in record["outputs"]):
            return expected_lookup(pd.read_csv(path, sep="\t")), "completed_qc"
    return expected_lookup(expected_table(clr, bins, mask, labels, views, values, c)), "computed_current_context"


def _reader_factory(clr):
    """Give each thread its own lazy read-only Cooler object and HDF5 handles."""
    if not getattr(clr, "uri", None):
        return lambda: clr  # Small adapter fixtures; caller must remain serial.
    import cooler
    local = threading.local()
    def reader():
        if not hasattr(local, "clr"):
            # oe_window requests raw counts and applies the inherited array.
            local.clr = cooler.Cooler(clr.uri)
        return local.clr
    return reader


def _chromosight_workers(c, per_worker_bytes, shared_bytes=0):
    threads = int(c.get("THREADS_CHROMOSIGHT", 4))
    budget = int(c.get("MEMORY_MB_CHROMOSIGHT", 16000))*1000000
    workers = min(threads, int(c.get("TOTAL_CPUS", threads)),
                  (int(budget*.8)-int(shared_bytes)) // max(int(per_worker_bytes), 1))
    require(workers >= 1, "Chromosight windows exceed MEMORY_MB_CHROMOSIGHT; lower maximum distance or increase its memory reservation")
    return workers


def chromosight_scores(c, scope):
    from importlib.resources import files
    from scipy.sparse import csr_matrix
    from chromosight.utils.detection import normxcorr2
    from .science import context, oe_window, save
    from .analyses import has_loop_discovery
    from .resolutions import quantification_resolutions
    version=importlib.metadata.version("chromosight")
    require(version=="1.6.3", "Chromosight adapter requires tested chromosight 1.6.3")
    kernel_path=files("chromosight.kernels").joinpath("artificial_template_loops_type1.txt")
    kernel=np.loadtxt(str(kernel_path)); half=kernel.shape[0]//2
    require(kernel.shape[0]==kernel.shape[1] and kernel.shape[0]%2==1, "Invalid loop kernel")
    root=Path(c["run_dir"]); outputs=[]
    for resolution in quantification_resolutions(c):
        clr,bins,mask,labels,views,values=context(c,scope,resolution)
        expected, expected_source = _expected_for(c,scope,resolution,clr,bins,mask,labels,views,values)
        directory=root/f"04_quantification/chromosight/{resolution}"
        source=root/f"04_quantification/discovered_loops/{resolution}/counts.tsv"
        if not has_loop_discovery(c): source=root/f"04_quantification/independent_loops/{resolution}/counts.tsv"
        features=pd.read_csv(source,sep="\t")
        reader = _reader_factory(clr)
        per_worker = (kernel.shape[0]+2)**2*8*32
        workers = _chromosight_workers(c, per_worker) if getattr(clr, "uri", None) else 1
        def score_feature(feature):
            raw = reader()
            row={"feature_id":feature["feature_id"],"sample_id":scope,"resolution":resolution,"score":np.nan,
                 "window_bins":kernel.shape[0],"missing_fraction":np.nan,"zero_fraction":np.nan,"status":"unavailable",
                 "PValue":np.nan,"FDR":np.nan}
            if not np.isfinite(feature[scope]):
                row["status"]="ineligible_candidate";return row
            span=(half+1)*resolution
            r1=(feature["chrom1"],int(feature["start1"])-span,int(feature["end1"])+span)
            r2=(feature["chrom2"],int(feature["start2"])-span,int(feature["end2"])+span)
            if r1[1]<0 or r2[1]<0 or r1[2]>raw.chromsizes[r1[0]] or r2[2]>raw.chromsizes[r2[0]]:
                row["status"]="edge_window";return row
            window=oe_window(raw,r1,r2,values,labels,views,expected,c)
            central=window[1:-1,1:-1];good=np.isfinite(central)
            missing=1-float(good.mean());row["missing_fraction"]=missing
            row["zero_fraction"]=float((central[good]==0).mean()) if good.any() else np.nan
            if missing>c["CHROMOSIGHT_MAX_MISSING_FRACTION"] or not good.any(): row["status"]="masked_window"
            elif np.std(central[good])==0: row["status"]="flat_window"
            else:
                missing_mask=csr_matrix(~np.isfinite(window),dtype=bool)
                signal=csr_matrix(np.nan_to_num(window,nan=0.0,posinf=0.0,neginf=0.0))
                corr,_=normxcorr2(signal,kernel,full=True,missing_mask=missing_mask,missing_tol=c["CHROMOSIGHT_MAX_MISSING_FRACTION"],pval=False)
                centre=window.shape[0]//2
                score=float(corr[centre,centre])
                row["score"]=score if np.isfinite(score) else np.nan
                row["status"]="measured" if np.isfinite(score) else "undefined"
            return row
        with ThreadPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(score_feature, features.to_dict("records")))
        outputs.append(save(pd.DataFrame(rows,columns=["feature_id","sample_id","resolution","score","window_bins","missing_fraction","zero_fraction","status","PValue","FDR"]),directory/f"{scope}.tsv"))
        provenance=directory/f"{scope}.json"
        atomic_json(provenance,{"software":"chromosight","version":version,"kernel_sha256":file_hash(str(kernel_path)),"normalization":"inherited_weights_and_common_view_observed_expected","adapter":"bounded_fixed_window_normxcorr2; adapted, not exact historical CLI reproduction","biological_significance":"not_computed","resolution":resolution,"sample_id":scope,"workers":workers,"expected_source":expected_source,"max_missing_fraction":c["CHROMOSIGHT_MAX_MISSING_FRACTION"],"input_mutation":False})
        outputs.append(provenance)
    return outputs


def chromosight_discovery(c, scope):
    """Independent CPU discovery on registered, inherited-weight matrices.

    Uses Chromosight's normalized correlation and focus picker in bounded
    rectangular tiles. It does not invoke the CLI's automatic normalization,
    alter source weights, or seed candidates from another caller. Tiled focus
    grouping and normalization are explicit adaptations of historical methods.
    """
    from importlib.resources import files
    from scipy.sparse import csr_matrix
    from chromosight.utils.detection import normxcorr2, pick_foci
    from .science import context, oe_window, save
    from .resolutions import chromosight_discovery_resolutions
    version = importlib.metadata.version("chromosight")
    require(version == "1.6.3", "Chromosight adapter requires tested chromosight 1.6.3")
    kernel_path = files("chromosight.kernels").joinpath("artificial_template_loops_type1.txt")
    kernel = np.loadtxt(str(kernel_path))
    require(kernel.ndim == 2 and kernel.shape[0] == kernel.shape[1] and kernel.shape[0]%2 == 1, "Invalid loop kernel")
    half, halo = kernel.shape[0]//2, kernel.shape[0]//2+1
    minimum = int(c.get("CHROMOSIGHT_DISCOVERY_MIN_DISTANCE_BP", 0) or c["LOOP_MIN_DISTANCE_BP"])
    maximum = int(c.get("CHROMOSIGHT_DISCOVERY_MAX_DISTANCE_BP", 0) or c["LOOP_MAX_DISTANCE_BP"])
    pearson = float(c.get("CHROMOSIGHT_DISCOVERY_MIN_PEARSON", .3))
    missing_tol = float(c["CHROMOSIGHT_MAX_MISSING_FRACTION"])
    # Native loops.json permits at most 10% observed zeros. Keep this caller
    # filter explicit; it is distinct from masked/undetected-bin tolerance.
    zero_tol = .1
    require(0 <= minimum < maximum and -1 <= pearson <= 1 and 0 <= missing_tol <= 1, "Invalid Chromosight discovery distance/correlation/missing settings")
    root, outputs = Path(c["run_dir"]), []
    for resolution in chromosight_discovery_resolutions(c):
        clr, bins, mask, labels, views, values = context(c, scope, resolution)
        expected, expected_source = _expected_for(c, scope, resolution, clr, bins, mask, labels, views, values)
        distance = maximum//resolution
        width = min(512, int(c["MAX_DENSE_BINS"])-distance-2*halo)
        require(width >= 1, "Chromosight distance band exceeds MAX_DENSE_BINS; lower maximum distance or use a coarser grid")
        per_worker = (width+2*halo)*(width+distance+2*halo)*8*32
        workers = _chromosight_workers(c, per_worker, values.nbytes+mask.nbytes+labels.nbytes)
        if not getattr(clr, "uri", None):
            workers = 1
        reader = _reader_factory(clr)
        tasks = []
        for vid, view in enumerate(views.itertuples()):
            lo, hi = clr.extent((view.chrom, int(view.start), int(view.end)))
            tasks.extend((vid, str(view.chrom), lo, hi, start, min(start+width, hi)) for start in range(lo, hi, width))
        def detect_tile(task):
            vid, chrom, lo, hi, start, end = task
            a, b = max(lo, start-halo), min(hi, end+halo)
            left, right = max(lo, start-halo), min(hi, end+distance+halo)
            if min(b-a, right-left) <= kernel.shape[0]:
                return []
            r1 = (chrom, int(bins.iloc[a].start), int(bins.iloc[b-1].end))
            r2 = (chrom, int(bins.iloc[left].start), int(bins.iloc[right-1].end))
            window = oe_window(reader(), r1, r2, values, labels, views, expected, c)
            rows = np.arange(a, b)[:, None]
            cols = np.arange(left, right)[None, :]
            valid = (np.isfinite(window) & (cols >= rows) & mask[a:b, None] & mask[None, left:right]
                     & (labels[a:b, None] == vid) & (labels[None, left:right] == vid))
            signal = csr_matrix(np.where(valid, window, 0.))
            missing = csr_matrix(~valid, dtype=bool)
            corr, logp = normxcorr2(signal, kernel, full=True, missing_mask=missing,
                                   missing_tol=missing_tol, pval=True)
            coords, _ = pick_foci(corr, pearson, min_size=2)
            if coords is None:
                return []
            found = []
            for i, j in coords:
                x, y = a+int(i), left+int(j)
                if not start <= x < end or not (mask[x] and mask[y] and labels[x] == labels[y] == vid):
                    continue
                anchor1, anchor2 = bins.iloc[x], bins.iloc[y]
                delta = int(anchor2.start)-int(anchor1.start)
                if y <= x or not minimum <= delta <= maximum or anchor1.end-anchor1.start != resolution or anchor2.end-anchor2.start != resolution:
                    continue
                # The entire kernel must fit the current registered view/tile.
                if i-half < 0 or j-half < 0 or i+half >= window.shape[0] or j+half >= window.shape[1]:
                    continue
                patch = window[i-half:i+half+1, j-half:j+half+1]
                good = valid[i-half:i+half+1, j-half:j+half+1]
                missing_fraction = 1-float(good.mean())
                zero_fraction = float((patch[good] == 0).mean()) if good.any() else np.nan
                # Match the native validator's strict fraction inequalities.
                if not good.any() or missing_fraction >= missing_tol or zero_fraction >= zero_tol:
                    continue
                score = float(corr[i, j])
                if not np.isfinite(score) or score < pearson:
                    continue
                pvalue = float(10.**float(logp[i, j]))
                found.append(dict(chrom1=chrom, start1=int(anchor1.start), end1=int(anchor1.end),
                                  chrom2=chrom, start2=int(anchor2.start), end2=int(anchor2.end),
                                  caller="chromosight", caller_score=score, caller_pvalue=pvalue,
                                  caller_qvalue=np.nan, caller_window_bins=kernel.shape[0],
                                  caller_missing_fraction=missing_fraction, caller_zero_fraction=zero_fraction,
                                  distance_bp=delta, resolution=resolution, inference_status="descriptive_same_data"))
            return found
        records = []
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for calls in executor.map(detect_tile, tasks):
                records.extend(calls)
        columns = LOOP_COLUMNS[:6]+["caller", "caller_score", "caller_pvalue", "caller_qvalue", "caller_window_bins",
                                    "caller_missing_fraction", "caller_zero_fraction", "distance_bp", "resolution", "inference_status"]
        frame = pd.DataFrame(records, columns=columns)
        if len(frame):
            frame = frame.sort_values(["caller_score"]+LOOP_COLUMNS[:6], ascending=[False]+[True]*6).drop_duplicates(LOOP_COLUMNS[:6])
        frame = canonical_loops(frame, c["REFERENCE_ID"], resolution)
        target = root/f"03_features/loops/{scope}.{resolution}.chromosight.tsv"
        outputs.append(save(frame, target))
        provenance = target.with_suffix(".json")
        atomic_json(provenance, {"software": "chromosight", "version": version, "sample_id": scope,
                    "reference_id": c["REFERENCE_ID"], "resolution": resolution,
                    "source_uri": getattr(clr, "uri", None), "kernel_sha256": file_hash(str(kernel_path)),
                    "normalization": "inherited_weights_and_common_view_observed_expected",
                    "expected_source": expected_source, "adapter": "bounded_rectangular_normxcorr2_and_native_pick_foci",
                    "exact_publication_reproduction": False, "independent_discovery": True,
                    "parameters": {"min_distance_bp": minimum, "max_distance_bp": maximum,
                                   "min_pearson": pearson, "max_missing_fraction": missing_tol,
                                   "max_zero_fraction": zero_tol, "focus_min_pixels": 2,
                                   "kernel_bins": list(kernel.shape), "centre_distance_filter_only": True,
                                   "neighbour_separation_bins": 1},
                    "geometry": {"tile_core_rows": width, "kernel_halo_bins": halo,
                                 "candidate_ownership": "first_anchor_core_row", "focus_grouping": "within_overlapping_tiles",
                                 "expected_ignored_diagonals": c["INSULATION_IGNORE_DIAGS"]},
                    "resources": {"workers": workers, "estimated_dense_bytes_per_worker": per_worker,
                                  "estimate_kind": "conservative_array_estimate_not_measured_peak_RSS",
                                  "memory_mb": c.get("MEMORY_MB_CHROMOSIGHT", 16000)},
                    "caller_statistics": "window correlation p-values are exploratory; caller_qvalue and biological significance are not computed",
                    "inference_status": "descriptive_same_data", "input_mutation": False})
        outputs.append(provenance)
    return outputs
