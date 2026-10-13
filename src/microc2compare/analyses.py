"""Named analysis selection; shared by normalization, execution and documentation."""
from .contracts import read_table, require

CATALOGUE = {
    "qc": "Contact totals and raw/balanced distance decay",
    "hicrep": "Raw-count pairwise contact-map concordance",
    "compartments": "Oriented compartment vectors and browser tracks",
    "insulation": "Insulation scores and boundary calls",
    "tads": "Per-replicate domains and descriptive TAD changes",
    "loops": "Loop discovery and consensus",
    "quantify": "Raw loop counts and eligible contact exposures",
    "chromosight": "Fixed-candidate pattern scores in biological replicates",
    "chromosight_discovery": "Independent CPU loop-pattern discovery on configured grids",
    "stripes": "CPU stripe detection and descriptive condition support",
    "browser_export": "Unified UCSC tracks and optional HiGlass conversions",
    "annotation": "Gene context and cCRE overlaps at both loop anchors",
    "apa": "Aggregate loop contact enrichment",
    "ada": "Aggregate domains and per-replicate domain strength",
    "saddle": "Compartment contact segregation",
    "comparisons": "Configured loop/insulation condition comparisons",
    "loop_comparison": "Only configured loop comparisons",
    "insulation_comparison": "Only configured insulation comparisons",
    "pool": "Explicit optional pooled discovery inputs",
    "depth_sensitivity": "Explicit optional sequencing-depth sensitivity",
    "report": "Report with configured analysis prerequisites",
}


def aggregate_modules(c):
    return set(c.get("AGGREGATE_MODULES", "apa,ada,saddle").split(",")) if c["RUN_AGGREGATES"] else set()


def requires_saddle(c):
    """The full fine-loop profiles include saddles; explicit module runs stay selective."""
    selected = set(c.get("_requested_analyses", []))
    return c.get("ANALYSIS_PROFILE") in ("practical_2_4kb", "hsieh_inspired", "hsieh_8grid", "short_range_4kb") and not (selected - {"report"})


def comparison_modules(c):
    if not c["RUN_DIFFERENTIAL"]:
        return set()
    rows = read_table(c["COMPARISONS"], ("modules",))
    allowed = set(c.get("_comparison_modules", ["loops", "insulation"]))
    return {m for row in rows for m in row["modules"].split(";")} & allowed


def configure_selection(c, request):
    """Set scientific switches before module-specific validation is performed."""
    if not request:
        return c
    names = set(request.split(","))
    require(names <= set(CATALOGUE) and "" not in names,
            "Unknown analysis selection; use --list-analyses")
    c["_requested_analyses"] = sorted(names)
    # A report-only request retains the configured science and refreshes/reuses
    # its prerequisites. Other selections produce a report unless disabled.
    if names == {"report"}:
        c["RUN_FINAL_REPORT"] = True
        return c
    original_pool, original_depth = c["RUN_POOL"], c["RUN_DEPTH_SENSITIVITY"]
    for key in ("RUN_QC", "RUN_HICREP", "RUN_POOL", "RUN_COMPARTMENTS", "RUN_INSULATION",
                "RUN_TADS", "RUN_LOOPS", "RUN_CHROMOSIGHT", "RUN_CHROMOSIGHT_DISCOVERY", "RUN_STRIPES", "RUN_BROWSER_EXPORT", "RUN_ANNOTATION", "RUN_AGGREGATES", "RUN_DIFFERENTIAL", "RUN_DEPTH_SENSITIVITY"):
        c[key] = False
    c["RUN_QC"] = "qc" in names
    c["RUN_HICREP"] = "hicrep" in names
    c["RUN_COMPARTMENTS"] = bool(names & {"compartments", "saddle"})
    c["RUN_INSULATION"] = bool(names & {"insulation", "tads", "ada", "insulation_comparison"})
    c["RUN_TADS"] = "tads" in names
    c["RUN_LOOPS"] = "loops" in names or bool(names & {"quantify", "apa", "chromosight", "annotation"}) and not c["INDEPENDENT_FEATURES_REGISTRY"] and "chromosight_discovery" not in names
    c["RUN_CHROMOSIGHT"] = "chromosight" in names
    c["RUN_CHROMOSIGHT_DISCOVERY"] = "chromosight_discovery" in names
    c["RUN_STRIPES"] = "stripes" in names
    c["RUN_BROWSER_EXPORT"] = "browser_export" in names
    c["RUN_ANNOTATION"] = "annotation" in names
    selected_aggregates = names & {"apa", "ada", "saddle"}
    c["RUN_AGGREGATES"] = bool(selected_aggregates)
    c["AGGREGATE_MODULES"] = ",".join(sorted(selected_aggregates)) if selected_aggregates else "apa,ada,saddle"
    if names & {"comparisons", "loop_comparison", "insulation_comparison"}:
        c["RUN_DIFFERENTIAL"] = True
        c["_comparison_modules"] = ["loops", "insulation"] if "comparisons" in names else sorted(
            {m for n, m in (("loop_comparison", "loops"), ("insulation_comparison", "insulation")) if n in names})
        modules = comparison_modules(c)
        require(modules, "No comparisons match the selected analysis")
        c["RUN_INSULATION"] |= "insulation" in modules
        rows = read_table(c["COMPARISONS"])
        # Formal loops quantify an independent universe; descriptive loops need
        # discovery unless a frozen external set was explicitly supplied.
        c["RUN_LOOPS"] |= "loops" in modules and not c["INDEPENDENT_FEATURES_REGISTRY"] and not c["RUN_CHROMOSIGHT_DISCOVERY"] and any(
            r["inference"] == "descriptive" and "loops" in r["modules"].split(";") for r in rows)
    if has_loop_discovery(c):
        c["RUN_QC"] = True
    if "pool" in names:
        require(original_pool, "--analysis pool requires RUN_POOL=true and explicit pool configuration")
        c["RUN_POOL"] = True
    if "depth_sensitivity" in names:
        require(original_depth, "--analysis depth_sensitivity requires RUN_DEPTH_SENSITIVITY=true")
        c["RUN_DEPTH_SENSITIVITY"] = True
        c["RUN_QC"] = True
        if c["DEPTH_SENSITIVITY_MODE"] == "balanced_features":
            c["RUN_INSULATION"] = c["RUN_LOOPS"] = True
    if "report" in names:
        c["RUN_FINAL_REPORT"] = True
    return c


def needs_loop_counts(c):
    return has_loop_discovery(c) or c.get("RUN_CHROMOSIGHT", False) or c.get("RUN_ANNOTATION", False) or "quantify" in c.get("_requested_analyses", []) or "loops" in comparison_modules(c)


def has_loop_discovery(c):
    return c["RUN_LOOPS"] or c.get("RUN_CHROMOSIGHT_DISCOVERY", False)


def resolution_set(c):
    """Required existing grids without scientific imports or contact-map I/O.

    Bootstrap state checks and scientific registration share this policy.
    """
    result = {c["PS_RESOLUTION"]} if c["RUN_QC"] or c["RUN_DEPTH_SENSITIVITY"] else set()
    for flag, key in (("RUN_COMPARTMENTS", "COMPARTMENT_RESOLUTION"), ("RUN_INSULATION", "INSULATION_RESOLUTION"), ("RUN_LOOPS", "LOOP_RESOLUTION"), ("RUN_HICREP", "HICREP_RESOLUTION")):
        if c[flag]: result.add(c[key])
    from .resolutions import loop_resolutions, quantification_resolutions
    if needs_loop_counts(c) or "apa" in aggregate_modules(c): result.update(quantification_resolutions(c))
    if "ada" in aggregate_modules(c): result.add(c.get("ADA_RESOLUTION", c["INSULATION_RESOLUTION"]))
    if c["RUN_LOOPS"]: result.update(loop_resolutions(c))
    if c.get("RUN_CHROMOSIGHT_DISCOVERY"):
        from .resolutions import chromosight_discovery_resolutions
        result.update(chromosight_discovery_resolutions(c))
    if c.get("RUN_STRIPES"): result.update(int(x) for x in c["STRIPE_RESOLUTIONS"].split(",") if x)
    return sorted(result or {c["PS_RESOLUTION"]})
