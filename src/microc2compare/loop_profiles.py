"""Assembly-independent loop presets and conservative caller admission plans.

Presets describe intended analyses; they never certify input depth, manufacture
fine pixels, approve resolution eligibility, or reproduce a historical result.
"""
from __future__ import annotations

import math
from .contracts import require

HSIEH_GRIDS = (400, 600, 800, 1000, 2000, 4000, 10000, 20000)
SHORT_RANGE_GRIDS = (400, 600, 800, 1000, 2000, 4000)


def profile_defaults(name):
    """Return a fresh config dict; sample identity and assembly stay user-owned."""
    aliases = {"hsieh": "hsieh_8grid",
               "short_range": "short_range_4kb", "short_range_profiles": "short_range_4kb"}
    name = aliases.get(name, name)
    require(name in ("hsieh_8grid", "short_range_4kb"), f"Unknown loop preset: {name}")
    grids = HSIEH_GRIDS if name == "hsieh_8grid" else SHORT_RANGE_GRIDS
    return {
        "ANALYSIS_PROFILE": name,
        "LOOP_CALLER": "mustache", "LOOP_RESOLUTION": 400,
        "LOOP_DISCOVERY_RESOLUTIONS": ",".join(map(str, grids)),
        "LOOP_SENSITIVITY_RESOLUTIONS": "",
        "LOOP_QUANTIFICATION_RESOLUTIONS": ",".join(map(str, grids)),
        "LOOP_MIN_DISTANCE_BP": 1000, "LOOP_MAX_DISTANCE_BP": 2000000 if name == "hsieh_8grid" else 500000,
        "MUSTACHE_P_THRESHOLD": .1, "MUSTACHE_SPARSITY_THRESHOLD": .88,
        "MUSTACHE_SIGMA_ZERO": 1.6, "MUSTACHE_OCTAVES": 2, "MUSTACHE_ITERATIONS": 10,
        "RUN_CHROMOSIGHT_DISCOVERY": True, "CHROMOSIGHT_DISCOVERY_RESOLUTIONS": ",".join(map(str, grids)),
        "CHROMOSIGHT_DISCOVERY_MIN_DISTANCE_BP": 0, "CHROMOSIGHT_DISCOVERY_MAX_DISTANCE_BP": 0,
        "CHROMOSIGHT_DISCOVERY_MIN_PEARSON": .3,
        "THREADS_CHROMOSIGHT": 4, "MEMORY_MB_CHROMOSIGHT": 16000,
        "RUN_CHROMOSIGHT": True,
        # 120000 is exactly divisible by every preset grid, including 600 bp.
        "APA_FLANK_BP": 120000,
    }


def mustache_geometry(c, resolution):
    """Best-case occupancy floor for the tested Mustache 1.3.3 code.

    The first assigned scale is sigma_zero * 2**(2/10).  Occupancy uses
    original nonzero pixels at diagonal >= 4, in radius ceil(scale) and twice
    that radius. This is a necessary best-case bound, never a detection promise.
    """
    require(int(resolution) > 0, "Resolution must be positive")
    require(c["MUSTACHE_ITERATIONS"] == 10, "Mustache 1.3.3 has ten fixed scale iterations")
    sigma = float(c["MUSTACHE_SIGMA_ZERO"])
    require(sigma > 0 and int(c["MUSTACHE_OCTAVES"]) > 0, "Mustache scales must be positive")
    threshold = float(c["MUSTACHE_SPARSITY_THRESHOLD"])
    require(0 <= threshold <= 1, "Mustache sparsity threshold must be in [0,1]")
    radius = math.ceil(sigma * 2 ** (2 / 10))
    def occupancy(distance, r):
        count = sum(distance + j - i >= 4 for i in range(-r, r+1) for j in range(-r, r+1))
        return count / (2*r+1)**2
    # c at diagonals <=4 is replaced by the constant 2 before the final
    # enrichment check, so offset 4 cannot survive that check either.
    minimum = next(d for d in range(5, 4*radius+5)
                   if occupancy(d, radius) >= threshold and occupancy(d, 2*radius) >= .6)
    configured = math.ceil(c["LOOP_MIN_DISTANCE_BP"] / resolution)
    return {"caller": "mustache", "version": "1.3.3", "resolution": int(resolution),
            "first_assigned_scale_bins": sigma * 2 ** (2 / 10), "inner_radius_bins": radius,
            "nonzero_diagonal_min_bins": 4, "inner_occupancy_threshold": threshold,
            "signal_candidate_min_distance_bins": 5,
            "outer_occupancy_threshold": .6, "geometric_min_distance_bins": minimum,
            "geometric_min_distance_bp": minimum * resolution,
            "effective_best_case_min_distance_bp": max(minimum, configured) * resolution,
            "occupancy_at_five_bins": {"inner": occupancy(5, radius), "outer": occupancy(5, 2*radius)},
            "occupancy_at_geometric_min": {"inner": occupancy(minimum, radius), "outer": occupancy(minimum, 2*radius)},
            "qualification": "ideal unmasked interior occupancy only; larger scales, zeros, masks, edges, enrichment and significance filters can raise the observed minimum"}


def mustache_resource_plan(c, resolution, sparse_bytes=0):
    """Admit fewer workers per grid before allocating dense Mustache tiles.

    MB denotes decimal MB, matching workflow reservations. Twenty float64
    arrays per dense worker is a conservative estimate, not a measured RSS cap.
    Reserve 20% for Python, sparse queries and caller overhead; subtract known
    global sparse arrays separately when supplied by the adapter.
    """
    distance = int(c["LOOP_MAX_DISTANCE_BP"]) // int(resolution)
    tile = max(2*distance, 2000)
    require(tile <= int(c["MAX_DENSE_BINS"]), "Mustache tile exceeds MAX_DENSE_BINS; lower maximum distance or use a coarser resolution")
    budget = int(c["MEMORY_MB_LOOPS"]) * 1000000
    per_worker = tile**2 * 8 * 20
    available = int(budget*.8) - int(sparse_bytes)
    workers = min(int(c["THREADS_LOOPS"]), int(c.get("TOTAL_CPUS", c["THREADS_LOOPS"])), available // per_worker)
    require(workers >= 1, f"Mustache {resolution} bp requires at least {math.ceil((per_worker+int(sparse_bytes))/.8/1000000)} MB; increase MEMORY_MB_LOOPS or lower maximum distance")
    return {"resolution": int(resolution), "tile_bins": tile, "maximum_distance_bins": distance,
            "requested_workers": int(c["THREADS_LOOPS"]), "workers": workers,
            "estimated_dense_bytes_per_worker": per_worker, "known_sparse_bytes": int(sparse_bytes),
            "reserved_overhead_fraction": .2, "memory_budget_bytes": budget,
            "estimate_kind": "conservative_array_estimate_not_measured_peak_RSS"}
