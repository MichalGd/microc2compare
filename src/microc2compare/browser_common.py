"""Table-only browser export primitives shared by every reference and cohort.

Signal/domain logic is adapted from this repository's 11 October 2026 exporters.
No scientific package or contact-matrix reader is imported here.
"""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def token(value):
    value = str(value)
    require(value and not any(ord(x) < 32 for x in value), "Empty/control-character identifier")
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._") or "label"
    reserved = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)), *(f"LPT{i}" for i in range(1, 10))}
    if cleaned != value or cleaned.split(".", 1)[0].upper() in reserved:
        cleaned = cleaned[:80] + "_" + hashlib.sha256(value.encode()).hexdigest()[:12]
    return cleaned


def text_value(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return format(value, ".17g") if math.isfinite(value) else "NA"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, ensure_ascii=False)
    return str(value)


class ExportContext:
    def __init__(self, run, out, config, samples, sizes):
        self.run, self.out = Path(run).resolve(), Path(out).resolve()
        self.config, self.samples = config, samples
        self.assembly = str(config["REFERENCE_ID"])
        self.sizes = dict(sizes)
        self.conditions = defaultdict(list)
        ids = set()
        for sample in samples:
            sid, condition = str(sample["sample_id"]), str(sample["condition"])
            require(sid not in ids and not sid.startswith("pool__"), "Duplicate/pooled biological observation")
            require(sample.get("reference_id", self.assembly) == self.assembly, "Sample reference mismatch")
            token(sid); token(condition)
            require("/" not in sid and "\\" not in sid and sid not in (".", ".."), "Sample ID cannot be a path")
            ids.add(sid)
            self.conditions[condition].append(sid)
        require(samples, "No biological observations in frozen metadata")
        self.notes, self.missing, self.tracks = [], [], []
        self.inputs, self.copies = {}, []

    def source(self, rel):
        path = Path(rel)
        return path if path.is_absolute() else self.run / path

    def exists(self, rel):
        return self.source(rel).is_file()

    def record(self, rel):
        path = self.source(rel).resolve()
        require(path.is_file(), f"Missing source file: {path}")
        require(self.out != path and self.out not in path.parents, "Export products cannot be scientific inputs")
        current = sha256(path)
        old = self.inputs.setdefault(str(path), current)
        require(current == old, f"Input changed during export: {path}")
        return path

    def json(self, rel):
        return json.loads(self.record(rel).read_text(encoding="utf-8-sig"))

    def table(self, rel, required=()):
        path = self.record(rel)
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            require(reader.fieldnames is not None and set(required) <= set(reader.fieldnames), f"Missing headers {required}: {path}")
            require(len(reader.fieldnames) == len(set(reader.fieldnames)), f"Duplicate headers: {path}")
            rows = list(reader)
        require(all(None not in row and all(v is not None for v in row.values()) for row in rows), f"Malformed TSV: {path}")
        return rows

    def destination(self, rel):
        path = (self.out / rel).resolve()
        require(self.out in path.parents, f"Output path escapes export directory: {rel}")
        require(not path.exists(), f"Refusing to overwrite export product: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def write_table(self, rel, rows, fields=None):
        rows = list(rows)
        fields = fields or list(dict.fromkeys(k for row in rows for k in row))
        require(fields, f"Empty table needs columns: {rel}")
        path = self.destination(rel)
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: text_value(row.get(key)) for key in fields})
        return path

    def copy_source(self, src_rel, dest_rel):
        source = self.record(src_rel)
        # Retain source tables separately from ready-to-load browser artifacts.
        target = self.destination(Path("source_tables") / dest_rel)
        shutil.copyfile(source, target)
        require(sha256(target) == self.inputs[str(source)], "Source copy verification failed")
        self.copies.append({"source": str(source), "path": str(target.relative_to(self.out)), "method": "copy"})
        return target

    def coords(self, row):
        chrom, start, end = (row[k] for k in ("chrom", "start", "end")) if isinstance(row, dict) else row[:3]
        require(str(start).isdigit() and str(end).isdigit(), "Browser coordinates must be integer, zero-based half-open")
        start, end = int(start), int(end)
        require(chrom in self.sizes and 0 <= start < end <= self.sizes[chrom], f"Invalid/out-of-scope interval: {chrom}:{start}-{end}")
        return str(chrom), start, end

    def add_missing(self, feature, reason):
        self.missing.append({"feature": feature, "reason": reason})

    def location(self, rel, analysis, scope, label):
        # Existing helper paths remain useful identifiers; the public layout is uniform.
        parts = Path(rel).parts
        scope_dir = "individual" if scope in ("individual", "sample") else "conditions" if scope == "condition" else "shared"
        prefix = Path(analysis) / scope_dir
        if scope_dir != "shared":
            prefix /= token(label)
        leaf = token(parts[-1])
        windows = [p for p in parts[:-1] if p.startswith("window_")]
        if windows and not leaf.startswith(windows[-1]):
            leaf = windows[-1] + "." + leaf
        return prefix / "ucsc" / leaf

    def add_track(self, path, raw_path, analysis, scope, label, format, source, description, rows, resolution="", window="", bedpe_path=""):
        self.tracks.append({"path": str(path or ""), "raw_path": str(raw_path), "analysis": analysis,
                            "scope": scope, "label": label, "format": format, "source": str(source),
                            "description": description, "rows": rows, "resolution": resolution, "window": window,
                            "assembly": self.assembly, "bedpe_path": str(bedpe_path),
                            "ucsc_status": "ready" if path else "disabled", "higlass_status": "disabled"})

    def _track(self, rel, rows, description, analysis, scope, label, source, resolution="", window="", kind="bedGraph"):
        rows = list(rows)
        rows.sort(key=lambda row: (*self.coords(row), str(row[3]) if len(row) > 3 else ""))
        previous = None
        for row in rows:
            chrom, start, end = self.coords(row)
            if kind == "bedGraph":
                require(len(row) == 4 and math.isfinite(float(row[3])), "bedGraph requires finite numeric values")
                require(previous is None or previous[0] != chrom or previous[2] <= start, f"Overlapping bedGraph intervals: {rel}")
                previous = chrom, start, end
            else:
                require(3 <= len(row) <= 12, "Unsupported BED column count")
                if len(row) >= 5:
                    require(str(row[4]).isdigit() and 0 <= int(row[4]) <= 1000, "BED score must be integer 0..1000")
                if len(row) >= 6:
                    require(row[5] in (".", "+", "-"), "Invalid BED strand")
                if len(row) >= 8:
                    require(start <= int(row[6]) <= int(row[7]) <= end, "Invalid BED thick coordinates")
                if len(row) >= 9:
                    require(re.fullmatch(r"\d{1,3},\d{1,3},\d{1,3}", str(row[8])) and all(0 <= int(v) <= 255 for v in str(row[8]).split(",")), "Invalid BED RGB")
            require(all(not any(x in str(v) for x in ("\n", "\r", "\t")) for v in row), "Control character in browser field")
        relative = self.location(rel, analysis, scope, label)
        stem = relative.name.removesuffix(".custom.txt")
        if stem.endswith((".bed", ".bedgraph")):
            stem = stem.rsplit(".", 1)[0]
        raw = self.destination(relative.with_name(stem + (".bedGraph" if kind == "bedGraph" else ".bed")))
        data = "".join("\t".join(text_value(v) for v in row) + "\n" for row in rows)
        raw.write_text(data, encoding="utf-8")
        custom = None
        if self.config.get("EXPORT_UCSC", True):
            custom = self.destination(relative.with_name(stem + ".custom.txt"))
            ident = "m2c_" + token(f"{analysis}_{scope}_{label}")[:140] + "_" + hashlib.sha256(str(relative).encode()).hexdigest()[:10]
            options = "visibility=full autoScale=on alwaysZero=on color=213,94,0 altColor=0,114,178" if kind == "bedGraph" else "visibility=pack itemRgb=On useScore=0"
            custom.write_text(f"track type={kind} name={json.dumps(ident)} description={json.dumps(description)} {options}\n" + data, encoding="utf-8")
        self.add_track(custom, raw, analysis, scope, label, kind, source, description, len(rows), resolution, window)
        if kind == "bed" and analysis in ("tads", "ada"):
            # A domain square is the same exact interval on both matrix axes.
            # Preserve the interval track independently; no averaging/flattening.
            paired = self.destination(relative.parent.parent / "higlass" / (stem + ".squares.bedpe"))
            paired.write_text("".join("\t".join(map(str, (*row[:3], *row[:3], row[3]))) + "\n" for row in rows), encoding="utf-8")
            self.add_track(None, paired, analysis, scope, label, "bedpe", source,
                           description + "; exact domain squares on contact matrix", len(rows), resolution, window, bedpe_path=paired)
        return custom or raw

    def bedgraph(self, rel, rows, description, analysis, scope, label, source, resolution="", window=""):
        return self._track(rel, rows, description, analysis, scope, label, source, resolution, window, "bedGraph")

    def bed(self, rel, rows, description, analysis, scope, label, source, resolution="", window=""):
        return self._track(rel, rows, description, analysis, scope, label, source, resolution, window, "bed")

    def interactions(self, rel, records, description, analysis, scope, label, source, resolution=""):
        records = list(records)
        rows, paired = [], []
        for record in records:
            axes = [self.coords((record[f"chrom{i}"], record[f"start{i}"], record[f"end{i}"])) for i in (1, 2)]
            fid = str(record["feature_id"])
            token(fid)
            color = record.get("color", "128,128,128")
            require(re.fullmatch(r"\d{1,3},\d{1,3},\d{1,3}", color) and all(0 <= int(v) <= 255 for v in color.split(",")), "Invalid interaction RGB")
            paired.append((*axes[0], *axes[1], fid))
            names = [record.get(f"anchor{i}_id", fid + f":anchor{i}") for i in (1, 2)]
            displays = [(axes[0][0], min(axes[0][1], axes[1][1]), max(axes[0][2], axes[1][2]))] if axes[0][0] == axes[1][0] else axes
            for display in displays:
                rows.append((*display, fid, 0, record.get("value", 0), record.get("category", "unannotated"), color,
                             *axes[0], names[0], ".", *axes[1], names[1], "."))
        rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
        paired.sort(key=lambda r: (r[0], r[1], r[3], r[4], r[6]))
        require(all(not any(c in str(v) for c in ("\n", "\r", "\t")) for row in rows for v in row), "Control character in interaction field")
        relative = self.location(rel, analysis, scope, label)
        stem = relative.name.removesuffix(".custom.txt").removesuffix(".interact.bed")
        raw = self.destination(relative.with_name(stem + ".interact.bed"))
        raw.write_text("".join("\t".join(text_value(v) for v in r) + "\n" for r in rows), encoding="utf-8")
        bedpe = self.destination(relative.with_name(stem + ".bedpe"))
        bedpe.write_text("".join("\t".join(text_value(v) for v in r) + "\n" for r in paired), encoding="utf-8")
        custom = None
        if self.config.get("EXPORT_UCSC", True):
            custom = self.destination(relative.with_name(stem + ".custom.txt"))
            ident = "m2c_" + hashlib.sha256(str(relative).encode()).hexdigest()[:16]
            custom.write_text(f"track type=interact name={json.dumps(ident)} description={json.dumps(description)} visibility=full\n" + raw.read_text(encoding="utf-8"), encoding="utf-8")
        self.add_track(custom, raw, analysis, scope, label, "interact", source, description, len(rows), resolution, bedpe_path=bedpe)
        return custom or raw

    def verify_inputs(self):
        for path, recorded in self.inputs.items():
            require(sha256(path) == recorded, f"Input changed during export: {path}")
