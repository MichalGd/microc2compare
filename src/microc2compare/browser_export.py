"""Generic browser export from frozen result tables; never opens contact matrices."""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import struct
import subprocess

from .browser_common import ExportContext, require, sha256, token
from .contracts import load_structured, native_path, resolve
from .browser_domains import export_domains
from .browser_loops import export_loops
from .browser_signals import export_signals


EXPORT_KEYS = ("EXPORT_UCSC", "EXPORT_HIGLASS", "BROWSER_OUTPUT_DIR", "HIGLASS_CHROMSIZES")


def _bool(value):
    if isinstance(value, bool):
        return value
    require(str(value).lower() in ("true", "false"), f"Expected boolean export option, found {value!r}")
    return str(value).lower() == "true"


def _json(ctx, rel, data):
    path = ctx.destination(rel)
    path.write_text(json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return path


def _sizes(path):
    values = {}
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        require(len(fields) == 2 and fields[0] not in values and fields[1].isdigit() and int(fields[1]) > 0, "Invalid/duplicate chromosome sizes")
        values[fields[0]] = int(fields[1])
    require(values, "Empty chromosome-size file")
    return values


def _write_sizes(ctx, relative, values):
    path = ctx.destination(relative)
    path.write_text("".join(f"{chrom}\t{length}\n" for chrom, length in values.items()), encoding="utf-8")
    return path


def _reference(config, run):
    registry_path = native_path(resolve(config["REFERENCE_REGISTRY"], config.get("config_dir", run)))
    require(registry_path.is_file(), f"Missing registered reference: {registry_path}")
    registry = load_structured(registry_path)
    require(registry.get("reference_id") == config["REFERENCE_ID"], "Reference registry identity mismatch")
    sizes_path = native_path(resolve(registry.get("chromsizes", ""), registry_path.parent))
    require(sizes_path.is_file(), "Registered chromosome-size file unavailable")
    require(registry.get("chromsizes_sha256") == sha256(sizes_path), "Reference chromosome-size fingerprint mismatch")
    return registry_path, sizes_path, _sizes(sizes_path)


def _profiles(ctx, full_sizes):
    profiles, assignments = {}, []
    def add(order, kind, sample_id="", resolution=""):
        require(order and len(order) == len(set(order)), "Empty/duplicate chromosome axis")
        require(set(order) <= set(full_sizes), "Native chromosome axis lacks registered lengths")
        values = {chrom: full_sizes[chrom] for chrom in order}
        identity = [[chrom, size] for chrom, size in values.items()]
        profile = token(ctx.assembly) + "_" + hashlib.sha256(json.dumps(identity, separators=(",", ":")).encode()).hexdigest()[:16]
        if profile not in profiles:
            path = _write_sizes(ctx, f"reference/profiles/{profile}.chromsizes.tsv", values)
            profiles[profile] = {"profile_id": profile, "assembly": ctx.assembly, "coord_system": profile,
                                 "chromosomes": identity, "path": path.relative_to(ctx.out).as_posix(), "sha256": sha256(path), "kind": kind}
        assignments.append({"profile_id": profile, "sample_id": sample_id, "resolution": resolution, "kind": kind})
        return profile
    analysis_profile = add(list(ctx.sizes), "analysis_grid")
    for resolution, order in getattr(ctx, "registered_orders", {}).items():
        add(order, "analysis_grid", resolution=resolution)
    for path in sorted((ctx.run / "00_metadata/masks").glob("*.scope.json")):
        scope = ctx.json(path)
        for sample in scope.get("samples", []):
            require(sample["sample_id"] in {s["sample_id"] for s in ctx.samples}, "Unknown sample in native-axis metadata")
            add(sample["native_chromosomes"], "matrix_native", sample["sample_id"], scope.get("resolution", ""))
    override = str(ctx.config.get("HIGLASS_CHROMSIZES", ""))
    if override:
        path = ctx.record(native_path(resolve(override, ctx.config.get("config_dir", ctx.run))))
        selected = _sizes(path)
        require(all(full_sizes.get(chrom) == size for chrom, size in selected.items()), "HiGlass chromosome profile length/reference mismatch")
        require(set(ctx.sizes) <= set(selected), "HiGlass chromosome profile omits exported chromosomes")
        analysis_profile = add(list(selected), "explicit_higlass")
    ctx.write_table("reference/profile_assignments.tsv", assignments, ("profile_id", "sample_id", "resolution", "kind"))
    _json(ctx, "reference/coordinate_profiles.json", {"assembly": ctx.assembly, "default_profile": analysis_profile,
                                                    "profiles": list(profiles.values()),
                                                    "coordinates": "zero-based half-open", "matrix_headers_read": False})
    return profiles, analysis_profile


def _stripes(ctx):
    directory = ctx.run / "03_features/stripes"
    files = sorted(directory.glob("*.tsv"))
    count, stripe_ids, stripe_counts = 0, set(), {}
    for sample in ctx.samples:
        sid = sample["sample_id"]
        for path in files:
            if not path.name.startswith(sid + "."):
                continue
            ending = path.name[len(sid) + 1:].removesuffix(".tsv")
            if not ending.isdigit():
                continue
            resolution = int(ending)
            rows = ctx.table(path, ("stripe_id", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "direction", "matrix_location"))
            stripe_counts[(sid, resolution)] = len(rows)
            require(len({r["stripe_id"] for r in rows}) == len(rows) and all(r["stripe_id"] for r in rows), "Empty/duplicate stripe IDs")
            rectangles, anchors, bodies = [], [], []
            for row in rows:
                stripe_ids.add((sid, resolution, row["stripe_id"]))
                require(row["direction"] in ("upstream", "downstream") and row["matrix_location"] in ("LT", "UT"), "Invalid stripe direction/triangle")
                require(row.get("sample_id", sid) == sid and row.get("reference_id", ctx.assembly) == ctx.assembly, "Stripe sample/reference mismatch")
                require(row.get("condition", sample["condition"]) == sample["condition"], "Stripe condition mismatch")
                require(int(row.get("resolution", resolution)) == resolution, "Stripe resolution mismatch")
                a = ctx.coords(tuple(row[f"{field}1"] for field in ("chrom", "start", "end")))
                b = ctx.coords(tuple(row[f"{field}2"] for field in ("chrom", "start", "end")))
                require(a[0] == b[0], "Stripe geometry must be cis")
                # Native StripePy rectangles include the diagonal seed in the body.
                # The body therefore may overlap its narrow anchor; do not trim it.
                seed = int(row.get("seed_bp", a[1]))
                require(a[1] <= seed < a[2], "Stripe seed outside narrow anchor")
                require((b[1] < seed and b[2] <= a[2]) if row["direction"] == "upstream"
                        else (b[2] > seed + resolution and b[1] >= a[1]), "Stripe body direction disagrees with anchor")
                for prefix, coords in (("anchor", a), ("body", b)):
                    for field, value in zip(("chrom", "start", "end"), coords):
                        require(f"{prefix}_{field}" not in row or str(row[f"{prefix}_{field}"]) == str(value), "Stripe coordinate aliases disagree")
                rectangles.append((*a, *b, row["stripe_id"]))
                anchors.append((*a, row["stripe_id"] + ":anchor"))
                bodies.append((*b, row["stripe_id"] + ":body"))
            output = Path("stripes") / "individual" / token(sid)
            paired = ctx.destination(output / "higlass" / f"rectangles_{resolution}.bedpe")
            paired.write_text("".join("\t".join(map(str, row)) + "\n" for row in rectangles), encoding="utf-8")
            ctx.add_track(None, paired, "stripes", "individual", sid, "bedpe", str(path),
                          "Exact narrow anchor/body rectangle; direction is genomic, not gene strand", len(rectangles), resolution, bedpe_path=paired)
            ctx.bed(f"anchor_{resolution}.bed", anchors, f"{sid}: stripe narrow anchors, {resolution} bp; descriptive calls", "stripes", "individual", sid, str(path), resolution)
            ctx.bed(f"body_{resolution}.bed", bodies, f"{sid}: stripe extended bodies, {resolution} bp; preserve triangle/direction in sidecars", "stripes", "individual", sid, str(path), resolution)
            ctx.copy_source(path, f"stripes/individual/{token(sid)}/{resolution}.tsv")
            count += 1
    support_files = [directory / name for name in ("summary.tsv", "matches.tsv", "summary.json")]
    if any(path.is_file() for path in support_files):
        require(all(path.is_file() for path in support_files), "Incomplete saved stripe support products")
        samples = {s["sample_id"]: s for s in ctx.samples}
        summaries = ctx.table(support_files[0], ("sample_id", "condition", "reference_id", "resolution", "stripe_count", "biological_replicates_in_condition"))
        for row in summaries:
            require(row["sample_id"] in samples, "Unknown sample in stripe support summary")
            sample = samples[row["sample_id"]]
            require(row["condition"] == sample["condition"] and row["reference_id"] == ctx.assembly, "Stripe support sample/reference mismatch")
            require(int(row["biological_replicates_in_condition"]) == len(ctx.conditions[row["condition"]]), "Stripe support denominator differs from frozen cohort")
            require(stripe_counts.get((row["sample_id"], int(row["resolution"]))) == int(row["stripe_count"]), "Stripe summary count differs from saved calls")
        matches = ctx.table(support_files[1], ("resolution", "sample1", "sample2", "stripe1", "stripe2", "direction"))
        for row in matches:
            resolution = int(row["resolution"])
            require(all((row[f"sample{i}"], resolution, row[f"stripe{i}"]) in stripe_ids for i in (1, 2)), "Stripe support link references an unknown saved call")
            require(row["direction"] in ("upstream", "downstream"), "Invalid stripe support direction")
        summary = ctx.json(support_files[2])
        require({k: set(v) for k, v in summary.get("condition_members", {}).items()} == {k: set(v) for k, v in ctx.conditions.items()}, "Stripe support cohort differs from frozen metadata")
        for path in support_files:
            ctx.copy_source(path, "stripes/support/" + path.name)
        ctx.notes.append("Stripe support sidecars contain pairwise many-to-many links and observation counts. No merged stripe geometry or differential significance is inferred.")
    if not count:
        ctx.add_missing("Stripes", "No saved optional stripe calls; no stripe discovery is performed by export")


def _binary_check(path, filetype, raw=None, sizes=None):
    if filetype == "bigwig":
        with path.open("rb") as stream:
            require(stream.read(4) in (b"\x26\xfc\x8f\x88", b"\x88\x8f\xfc\x26"), "Converter output is not bigWig")
        try:
            import pyBigWig
        except ImportError:
            return "magic_checked; decoding/round-trip and rendering not verified"
        def normalized(records):
            result = []
            for start, end, value in records:
                if result and result[-1][1] == start and result[-1][2] == value:
                    result[-1] = (result[-1][0], end, value)
                else:
                    result.append((start, end, value))
            return result
        expected = {}
        for line in raw.read_text(encoding="utf-8").splitlines():
            chrom, start, end, value = line.split("\t")
            try:
                # bigWig stores float32 values; compare their defined rounding.
                value = struct.unpack("<f", struct.pack("<f", float(value)))[0]
            except OverflowError as error:
                raise ValueError("Signal cannot be represented as bigWig float32") from error
            expected.setdefault(chrom, []).append((int(start), int(end), value))
        reader = pyBigWig.open(str(path))
        require(reader is not None, "Cannot decode converted bigWig")
        try:
            require(reader.isBigWig(), "Converter output is not a decoded bigWig")
            observed_sizes = reader.chroms()
            require(all(sizes.get(chrom) == size for chrom, size in observed_sizes.items()) and set(expected) <= set(observed_sizes), "Converted bigWig chromosome header differs from reference")
            for chrom in observed_sizes:
                require(normalized(list(reader.intervals(chrom) or ())) == normalized(expected.get(chrom, [])), "Converted bigWig intervals/float32 values differ from source")
        finally:
            reader.close()
        return "header_and_float32_interval_roundtrip_checked; rendering not verified"
    with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as connection:
        require(connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok", "Clodius SQLite integrity failure")
        require(connection.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall(), "Clodius database has no tables")
    return "sqlite_integrity_checked; base-record round-trip and rendering not verified"


def _higlass(ctx, profiles, default_profile):
    if not ctx.config["EXPORT_HIGLASS"]:
        return
    bigwig = shutil.which("bedGraphToBigWig")
    clodius = shutil.which("clodius")
    commands, receipts = [], []
    for profile in profiles.values():
        relative = profile["path"]
        commands.append(["higlass-manage", "ingest", relative, "--filetype", "chromsizes-tsv", "--datatype", "chromsizes", "--assembly", profile["coord_system"]])
    for track in ctx.tracks:
        paired = track.get("bedpe_path", "")
        filetype = "bigwig" if track["format"] == "bedGraph" else "bed2ddb" if paired else "beddb"
        tool = bigwig if filetype == "bigwig" else clodius
        track["higlass_filetype"] = filetype
        track["higlass_datatype"] = "vector" if filetype == "bigwig" else "2d-rectangle-domains" if paired else "bedlike"
        if not track["rows"]:
            track["higlass_status"] = "empty; conversion not required"
            continue
        if tool is None:
            track["higlass_status"] = "converter_missing: " + ("bedGraphToBigWig" if filetype == "bigwig" else "clodius")
            continue
        raw = Path(paired or track["raw_path"])
        chromosomes = set()
        for line in raw.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            chromosomes.add(fields[0])
            if paired:
                chromosomes.add(fields[3])
        compatible = [pid for pid, profile in profiles.items() if chromosomes <= {r[0] for r in profile["chromosomes"]}]
        require(default_profile in compatible, "Default browser profile omits track chromosomes")
        track["higlass_profiles"] = ";".join(compatible)
        # bigWig chromosomes are named internally; interval DB offsets depend on order.
        track_profiles = [default_profile] if filetype == "bigwig" else compatible
        outputs = []
        for profile_id in track_profiles:
            profile = profiles[profile_id]
            coords = ctx.out / profile["path"]
            directory = Path(track["raw_path"]).parent
            if directory.name == "ucsc":
                directory = directory.parent / "higlass"
            output = ctx.destination(directory.relative_to(ctx.out) / profile_id / (raw.stem + "." + ("bw" if filetype == "bigwig" else filetype)))
            if filetype == "bigwig":
                command = [tool, str(raw), str(coords), str(output)]
            elif paired:
                command = [tool, "aggregate", "bedpe", "--chromsizes-filename", str(coords), "--chr1-col", "1", "--from1-col", "2", "--to1-col", "3",
                           "--chr2-col", "4", "--from2-col", "5", "--to2-col", "6", "--output-file", str(output), str(raw)]
            else:
                command = [tool, "aggregate", "bedfile", "--chromsizes-filename", str(coords), "--delimiter", "\t", "--output-file", str(output), str(raw)]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            if not output.is_file() and Path(str(output) + ".db").is_file():
                Path(str(output) + ".db").rename(output)
            require(output.is_file() and output.stat().st_size > 0, "Converter produced no nonempty output")
            validation = _binary_check(output, filetype, raw, dict(profile["chromosomes"]))
            log = ctx.destination(str(output.relative_to(ctx.out)) + ".conversion.log")
            log.write_text("Exit code: 0\nCommand argv: " + json.dumps(command) + "\n\n" +
                           completed.stdout + completed.stderr, encoding="utf-8")
            relative = output.relative_to(ctx.out).as_posix()
            outputs.append(relative)
            registration_profiles = compatible if filetype == "bigwig" else [profile_id]
            for registration in registration_profiles:
                commands.append(["higlass-manage", "ingest", relative, "--filetype", filetype, "--datatype", track["higlass_datatype"], "--assembly", profiles[registration]["coord_system"]])
            receipts.append({"path": relative, "profile_id": profile_id, "filetype": filetype,
                             "command": command, "validation": validation, "ingested": False, "render_tested": False})
        track["higlass_paths"] = ";".join(outputs)
        track["higlass_status"] = "converted; not ingested or render-tested"
    _json(ctx, "higlass/conversion_receipts.json", receipts)
    # A reviewable, structured recipe: no server writes and no shell string interpolation.
    _json(ctx, "higlass/ingestion_commands.json", {"working_directory": ".", "commands": commands,
                                                "note": "Choose the profile matching each displayed matrix. These commands were not executed."})


def export_browser(c, scope="run"):
    """Return exported artifact Paths. Any new cohort uses its own frozen metadata.

The default destination is a new UTC-named child of run/browser_tracks. An
explicit BROWSER_OUTPUT_DIR denotes a new bundle directory, never a reuse target.
"""
    require(scope == "run", "Browser export has one run scope")
    run = native_path(c["run_dir"]).resolve()
    config_path, samples_path = run / "00_metadata/config.json", run / "00_metadata/samples.json"
    require(config_path.is_file() and samples_path.is_file(), "Browser export requires frozen run config and sample metadata")
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    require(config.get("REFERENCE_ID"), "Frozen run reference ID is missing")
    frozen_registry = native_path(resolve(config["REFERENCE_REGISTRY"], config.get("config_dir", run))).resolve()
    frozen_inputs_path = run / "00_metadata/inputs.json"
    frozen_inputs = json.loads(frozen_inputs_path.read_text(encoding="utf-8-sig")) if frozen_inputs_path.is_file() else {}
    require(not c.get("REFERENCE_ID") or c["REFERENCE_ID"] == config["REFERENCE_ID"], "Active/frozen reference mismatch")
    for key in EXPORT_KEYS:
        if key in c:
            config[key] = c[key]
    config["EXPORT_UCSC"] = _bool(config.get("EXPORT_UCSC", True))
    config["EXPORT_HIGLASS"] = _bool(config.get("EXPORT_HIGLASS", False))
    require(config["EXPORT_UCSC"] or config["EXPORT_HIGLASS"], "Enable at least one browser export format")
    # Permit relocation of the reference registry without changing scientific identity.
    if c.get("REFERENCE_REGISTRY"):
        config["REFERENCE_REGISTRY"] = c["REFERENCE_REGISTRY"]
    config["config_dir"] = c.get("config_dir", config.get("config_dir", str(run)))
    metadata = json.loads(samples_path.read_text(encoding="utf-8-sig"))
    registry, sizes_path, full_sizes = _reference(config, run)
    expected_registry = next((value for key, value in frozen_inputs.items() if native_path(key.split("::", 1)[0]).resolve() == frozen_registry), None)
    if frozen_registry.is_file():
        original = load_structured(frozen_registry)
        require(original.get("reference_id") == config["REFERENCE_ID"] and original.get("chromsizes_sha256") == sha256(sizes_path), "Reference chromosome sizes differ from frozen registry")
        if expected_registry:
            require(expected_registry.get("sha256") == sha256(frozen_registry), "Frozen reference registry fingerprint changed")
    else:
        # A copied run may no longer reach the server's original asset paths.
        # Bind relocated content to saved preflight fingerprints, never merely
        # to the replacement registry's own newly declared hash.
        expected_hashes = {value.get("sha256") for value in frozen_inputs.values() if isinstance(value, dict)}
        require(sha256(sizes_path) in expected_hashes, "Relocated reference needs matching frozen input fingerprints")
    selected, selected_order, registered_orders = set(), [], {}
    registered = sorted((run / "00_metadata/masks").glob("*.bins.tsv"))
    for path in registered:
        import csv
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            require(reader.fieldnames and set(("chrom", "start", "end")) <= set(reader.fieldnames), "Registered grid columns missing")
            prior, seen = {}, set()
            order = []
            for row in reader:
                require(row["chrom"] in full_sizes and str(row["start"]).isdigit() and str(row["end"]).isdigit(), "Invalid registered chromosome grid")
                start, end = int(row["start"]), int(row["end"])
                require(0 <= start < end <= full_sizes[row["chrom"]], "Registered bins exceed reference bounds")
                require((row["chrom"], start, end) not in seen and prior.get(row["chrom"], -1) <= start, "Duplicate/overlapping registered bins")
                seen.add((row["chrom"], start, end)); prior[row["chrom"]] = end
                selected.add(row["chrom"])
                if row["chrom"] not in order:
                    order.append(row["chrom"])
                if row["chrom"] not in selected_order:
                    selected_order.append(row["chrom"])
            require(order, "Registered chromosome grid is empty")
            registered_orders[int(path.name.split(".", 1)[0])] = order
    sizes = {chrom: full_sizes[chrom] for chrom in (selected_order or list(full_sizes))}
    value = str(config.get("BROWSER_OUTPUT_DIR", ""))
    out = native_path(resolve(value, config["config_dir"])) if value else run / "browser_tracks" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    out = out.resolve()
    require(not out.exists(), f"Browser output must be new: {out}")
    require(out != run and out not in run.parents, "Browser output cannot replace or contain scientific run")
    # Validate metadata before creating output; failed export never receives a passed audit.
    ctx = ExportContext(run, out, config, metadata.get("samples", []), sizes)
    ctx.registered_orders = registered_orders
    out.mkdir(parents=True, exist_ok=False)
    try:
        for path in (config_path, samples_path, registry, sizes_path, *registered):
            ctx.record(path)
        if frozen_inputs_path.is_file():
            ctx.record(frozen_inputs_path)
        if frozen_registry.is_file():
            ctx.record(frozen_registry)
        _write_sizes(ctx, "reference/chromsizes.tsv", full_sizes)
        _write_sizes(ctx, "reference/analysis.chromsizes.tsv", sizes)
        profiles, default_profile = _profiles(ctx, full_sizes)
        export_loops(ctx, metadata.get("comparisons", []))
        export_signals(ctx)
        export_domains(ctx)
        _stripes(ctx)
        _higlass(ctx, profiles, default_profile)
        ctx.verify_inputs()
        inventory = []
        for row in ctx.tracks:
            row = dict(row)
            for key in ("path", "raw_path", "bedpe_path"):
                row[key] = Path(row[key]).relative_to(out).as_posix() if row.get(key) else ""
            row["coordinate_profile"] = default_profile
            inventory.append(row)
        fields = ("path", "raw_path", "bedpe_path", "analysis", "scope", "label", "format", "assembly", "coordinate_profile", "resolution", "window", "rows",
                  "source", "description", "ucsc_status", "higlass_status", "higlass_filetype", "higlass_datatype", "higlass_profiles", "higlass_paths")
        ctx.write_table("inventory.tsv", inventory, fields)
        ctx.write_table("samples.tsv", ctx.samples)
        ctx.write_table("missing_features.tsv", ctx.missing, ("feature", "reason"))
        ctx.write_table("source_files.tsv", [{"source": path, "sha256": digest} for path, digest in sorted(ctx.inputs.items())], ("source", "sha256"))
        counts = dict(Counter(row["analysis"] for row in inventory))
        description = [f"# Browser tracks: {config.get('RUN_ID', run.name)}\n", f"Assembly: **{ctx.assembly}**. Frozen run: `{run}`.\n",
                       "All browser coordinates are zero-based half-open. No matrices were opened, scientific stages rerun, weights changed, or files published.\n",
                       "Upload nonempty `ucsc/*.custom.txt` files using UCSC Custom Tracks for the matching assembly. Raw BED/BEDPE and TSV/JSON sidecars preserve geometry and values.\n",
                       "`inventory.tsv` lists observations, conditions, windows, source identities and per-format status. `reference/coordinate_profiles.json` binds ordered chromosome names/lengths; HiGlass offsets must match the displayed matrices.\n",
                       "New compartment/insulation/ADA condition means give each original observation equal weight and require complete finite coverage. Counts and available-case summaries remain separate. TAD summaries retain saved means across supporting called domains.\n",
                       "Condition detection is not biological absence or differential significance. Raw compartment eigenvectors retain arbitrary signs; selected oriented means are descriptive, not pooled eigenvectors.\n",
                       "HiGlass conversion is optional. Converter availability, conversion, ingestion and actual rendering are distinct. JSON ingestion commands are reviewable recipes and were not executed.\n",
                       "Track counts: " + json.dumps(counts, sort_keys=True) + ".\n",
                       "\n".join(ctx.notes), "\nUnavailable/qualified outputs:\n" + "\n".join(f"- {row['feature']}: {row['reason']}" for row in ctx.missing)]
        ctx.destination("README.md").write_text("\n".join(description) + "\n", encoding="utf-8")
        audit = {"status": "passed", "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(), "run": str(run),
                 "assembly": ctx.assembly, "condition_members": dict(ctx.conditions), "tracks": len(inventory), "track_counts": counts,
                 "input_sha256": ctx.inputs, "copies": ctx.copies, "coordinate_profile": default_profile,
                 "contact_matrices_opened": False, "scientific_analyses_rerun": False, "published": False,
                 "higlass_ingested": False, "higlass_render_tested": False}
        _json(ctx, "export_audit.json", audit)
        artifacts = [{"path": p.relative_to(out).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in sorted(out.rglob("*")) if p.is_file()]
        ctx.write_table("artifact_checksums.tsv", artifacts, ("path", "bytes", "sha256"))
        ctx.verify_inputs()
        return sorted(p for p in out.rglob("*") if p.is_file())
    except Exception as error:
        # Preserve diagnostics; do not imply completion after a partially written bundle.
        passed = out / "export_audit.json"
        if passed.is_file():
            passed.unlink()
        (out / "export_failure.json").write_text(json.dumps({"status": "failed", "error": str(error), "run": str(run)}, indent=2) + "\n", encoding="utf-8")
        raise
