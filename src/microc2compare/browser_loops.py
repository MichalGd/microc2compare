"""Export saved catalogue and native loop geometry without calling a caller."""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import math
from pathlib import Path

from .browser_common import require, token


COORDS = ("chrom1", "start1", "end1", "chrom2", "start2", "end2")
COLORS = {"promoter_enhancer_like": "230,159,0", "promoter_promoter": "0,114,178",
          "enhancer_like_enhancer_like": "0,158,115", "other_or_unresolved": "128,128,128", "unannotated": "128,128,128"}


def indexed(rows, field, what):
    result = {}
    for row in rows:
        key = row[field]
        require(key and key not in result, f"Empty/duplicate {what}: {key!r}")
        token(key)
        result[key] = row
    return result


def anchor_id(assembly, coordinates):
    return "export_anchor_" + hashlib.sha256(repr((assembly, *coordinates)).encode()).hexdigest()[:20]


def _records(ctx, rows, annotations=None):
    records = []
    for row in rows:
        record = dict(row)
        ann = (annotations or {}).get(row["feature_id"], {})
        category = ann.get("candidate_category", "unannotated")
        record.update(category=category, color=COLORS.get(category, "128,128,128"))
        for n in (1, 2):
            coords = ctx.coords(tuple(row[f"{field}{n}"] for field in ("chrom", "start", "end")))
            for field, value in zip(("chrom", "start", "end"), coords):
                record[f"{field}{n}"] = value
                require(f"anchor{n}_{field}" not in ann or str(ann[f"anchor{n}_{field}"]) == str(value), "Annotation/catalogue coordinate mismatch")
            record[f"anchor{n}_id"] = ann.get(f"anchor{n}_id", anchor_id(ctx.assembly, coords))
        records.append(record)
    return records


def _track_set(ctx, relative, records, scope, label, source, description, resolution=""):
    ctx.interactions(relative, records, description, "loops", scope, label, source, resolution)
    anchors = {}
    for row in records:
        for n in (1, 2):
            aid = row[f"anchor{n}_id"]
            coords = tuple(row[f"{field}{n}"] for field in ("chrom", "start", "end"))
            require(aid not in anchors or anchors[aid] == coords, "Anchor ID reused for different coordinates")
            anchors[aid] = coords
    ctx.bed(str(Path(relative).with_name(Path(relative).name.removesuffix(".custom.txt") + ".anchors.custom.txt")),
            [(*coords, aid) for aid, coords in sorted(anchors.items())],
            description + "; exact unique anchors", "loop_anchors", scope, label, source, resolution)


def _effects(ctx, candidates, annotations, comparisons):
    for comparison in comparisons:
        comp_id = str(comparison["comparison_id"])
        require("/" not in comp_id and "\\" not in comp_id and comp_id not in (".", ".."), "Comparison ID cannot be a path")
        directory = ctx.run / "06_comparisons" / comp_id
        for effect_path in sorted(directory.glob("*/loops.tsv")):
            resolution = int(effect_path.parent.name)
            params = ctx.json(effect_path.parent / "parameters.json")
            if params.get("loop_universe") != "discovered_loops":
                ctx.add_missing(f"{comp_id}/{resolution} loop effects", "Independent measurement universe; not remapped to discovery catalogue")
                continue
            require(int(params.get("measurement_resolution", resolution)) == resolution, "Effect measurement resolution mismatch")
            effects = indexed(ctx.table(effect_path, ("feature_id",)), "feature_id", "measurement ID")
            source = f"04_quantification/discovered_loops/{resolution}/candidate_mapping.tsv"
            mapping = ctx.table(source, ("candidate_id", "measurement_id", "resolution", "original_resolution", "original_start1", "original_start2", "collapsed_candidate_count"))
            aliases = indexed(mapping, "candidate_id", "candidate alias")
            require(set(aliases) == set(candidates), "Effect candidate mapping does not cover catalogue")
            frequencies = Counter(row["measurement_id"] for row in mapping)
            displayed, details = [], []
            for cid, alias in aliases.items():
                original = candidates[cid]
                require(int(alias["resolution"]) == resolution and int(alias["original_resolution"]) == int(original["representative_resolution"]), "Effect alias resolution mismatch")
                require(all(int(alias[f"original_start{n}"]) == int(original[f"start{n}"]) for n in (1, 2)), "Effect alias original coordinates mismatch")
                require(int(alias["collapsed_candidate_count"]) == frequencies[alias["measurement_id"]], "Effect alias collision count mismatch")
                require(alias["measurement_id"] in effects, "Missing effect measurement")
                measured = effects[alias["measurement_id"]]
                field = "log2_allocation_ratio" if "log2_allocation_ratio" in measured else "logFC"
                require(field in measured, "No recorded effect column")
                try:
                    value = float(measured[field])
                except ValueError:
                    value = math.nan
                details.append({"candidate_id": cid, "comparison": comp_id, "resolution": resolution,
                                "effect_column": field, "effect": measured[field], "mapping": alias, "measurement": measured})
                if math.isfinite(value):
                    record = _records(ctx, [original], annotations)[0]
                    record.update(value=abs(value), color="213,94,0" if value > 0 else "0,114,178" if value < 0 else "128,128,128",
                                  category=f"{comp_id}|{field}={value:.8g}|measurement={alias['measurement_id']}|aliases={alias['collapsed_candidate_count']}")
                    displayed.append(record)
            ctx.write_table(f"loops/comparisons/{token(comp_id)}/effects_{resolution}.tsv", details,
                            ("candidate_id", "comparison", "resolution", "effect_column", "effect", "mapping", "measurement"))
            ctx.interactions(f"effects_{token(comp_id)}_{resolution}.custom.txt", displayed,
                             f"{comp_id}: descriptive effect; repeated aliases are one measurement; no significance encoded",
                             "loops", "comparison", comp_id, str(effect_path), resolution)


def export_loops(ctx, comparisons=()):
    catalogue_source = "03_features/loops/catalogue.tsv"
    sample_ids = {row["sample_id"] for row in ctx.samples}
    directory = ctx.run / "03_features/loops"
    native_files = sorted(directory.glob("*.bedpe.tsv"))
    # Native per-grid calls survive even when no consensus catalogue exists yet.
    for sample in ctx.samples:
        sid = sample["sample_id"]
        for path in native_files:
            if not path.name.startswith(sid + "."):
                continue
            suffix = path.name[len(sid) + 1:].removesuffix(".bedpe.tsv")
            if not suffix.isdigit():
                continue
            resolution = int(suffix)
            require(resolution > 0, "Native call resolution must be positive")
            native = ctx.table(path, (*COORDS, "feature_id"))
            indexed(native, "feature_id", "native loop ID")
            records = _records(ctx, native)
            for row in native:
                require("resolution" not in row or int(row["resolution"]) == resolution, "Native call resolution mismatch")
            _track_set(ctx, f"native_{resolution}.custom.txt", records, "individual", sid, str(path),
                       f"{sid}: native {resolution}-bp caller coordinates; no catalogue-coordinate substitution", resolution)
            ctx.copy_source(path, f"loops/individual/{token(sid)}/native_{resolution}.tsv")
            if path.with_name(f"{sid}.{resolution}.json").is_file():
                ctx.copy_source(path.with_name(f"{sid}.{resolution}.json"), f"loops/individual/{token(sid)}/native_{resolution}.json")
    if not ctx.exists(catalogue_source):
        ctx.add_missing("Loop catalogue", "No saved catalogue; available native calls exported without reconciliation")
        return
    catalogue = ctx.table(catalogue_source, (*COORDS, "feature_id", "representative_resolution"))
    candidates = indexed(catalogue, "feature_id", "catalogue candidate ID")
    annotation_source = "03_features/annotation/loops.tsv"
    annotations = {}
    if ctx.exists(annotation_source):
        annotations = indexed(ctx.table(annotation_source, ("feature_id", "anchor1_id", "anchor2_id", "candidate_category")), "feature_id", "annotation ID")
        require(set(annotations) == set(candidates), "Catalogue and annotation IDs differ")
        ctx.copy_source(annotation_source, "loops/annotations.tsv")
    else:
        ctx.add_missing("Loop functional annotation", "Annotation unavailable; geometry retained with unannotated styling and export-derived anchor IDs")
    records = _records(ctx, catalogue, annotations)
    ctx.copy_source(catalogue_source, "loops/catalogue.tsv")
    _track_set(ctx, "all_loops.custom.txt", records, "shared", "catalogue", catalogue_source, "Original catalogue representatives; descriptive annotations")
    for category in sorted({row["category"] for row in records}):
        ctx.interactions(f"category_{token(category)}.custom.txt", [r for r in records if r["category"] == category],
                         f"Catalogue {category}; candidate annotation, not regulatory causality", "loops", "shared", "catalogue", catalogue_source)
    membership_source = "03_features/loops/source_membership.tsv"
    if ctx.exists(membership_source):
        membership = ctx.table(membership_source, ("candidate_id", "source", "resolution", *COORDS))
        support = defaultdict(set)
        pool_sources = set()
        for member in membership:
            require(member["candidate_id"] in candidates, "Membership references unknown candidate")
            require(int(member["resolution"]) > 0, "Invalid membership resolution")
            for n in (1, 2):
                ctx.coords(tuple(member[f"{key}{n}"] for key in ("chrom", "start", "end")))
            if member["source"] in sample_ids:
                support[member["candidate_id"]].add(member["source"])
            else:
                require(member["source"].startswith("pool__"), "Membership references unknown biological observation")
                pool_sources.add(member["source"])
        ctx.copy_source(membership_source, "loops/source_membership.tsv")
        summaries = [{"candidate_id": cid, "sample_ids": ";".join(sorted(support[cid])),
                      **{f"support:{condition}": len(support[cid] & set(members)) for condition, members in ctx.conditions.items()}}
                     for cid in candidates]
        ctx.write_table("loops/shared/detection_support.tsv", summaries,
                        ("candidate_id", "sample_ids", *(f"support:{condition}" for condition in ctx.conditions)))
        for sample in ctx.samples:
            sid = sample["sample_id"]
            _track_set(ctx, "detected.custom.txt", [r for r in records if sid in support[r["feature_id"]]], "individual", sid, membership_source,
                       f"{sid}: detected catalogue representatives; distinct observation support, not native localization")
        for condition, members in ctx.conditions.items():
            for minimum, stem in ((1, "union"), (2, "support_ge2")):
                if minimum > len(members):
                    continue
                selected = [r for r in records if len(support[r["feature_id"]] & set(members)) >= minimum]
                _track_set(ctx, stem + ".custom.txt", selected, "condition", condition, membership_source,
                           f"{condition}: detection support >= {minimum}/{len(members)}; not specificity or statistical significance")
    else:
        ctx.add_missing("Loop detection support", "Membership unavailable; missing observations are not treated as zero support")
    _effects(ctx, candidates, annotations, comparisons)
    ctx.notes.append("Native caller/grid loop coordinates and catalogue representatives are separate tracks. Biological support deduplicates sample IDs across grids and excludes pools. Geometry is never averaged; caller detection is not differential significance.")
