#!/usr/bin/env python3
"""Export any microc2compare run's saved tables to UCSC and optional HiGlass files."""
import argparse
import json
from pathlib import Path
import sys
import subprocess

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from microc2compare.browser_export import export_browser


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--config-json", type=Path, help="Resolved workflow config containing run_dir")
    selection.add_argument("--run", type=Path, help="Run containing frozen 00_metadata config and samples")
    parser.add_argument("--out", type=Path, help="New bundle directory; default run/browser_tracks/<UTC>")
    parser.add_argument("--reference-registry", type=Path, help="Relocated registered reference with unchanged ID and chromosome-size digest")
    parser.add_argument("--higlass-chromsizes", type=Path, help="Explicit ordered chromosome profile matching displayed matrices")
    parser.add_argument("--higlass", action="store_true", help="Opt-in local format conversion; no ingestion or publication")
    parser.add_argument("--no-ucsc", action="store_true", help="Omit UCSC custom headers; raw geometry is retained")
    args = parser.parse_args()
    config = json.loads(args.config_json.read_text(encoding="utf-8-sig")) if args.config_json else {"run_dir": str(args.run.resolve())}
    if args.out:
        config["BROWSER_OUTPUT_DIR"] = str(args.out.resolve())
    if args.reference_registry:
        config["REFERENCE_REGISTRY"] = str(args.reference_registry.resolve())
    if args.higlass_chromsizes:
        config["HIGLASS_CHROMSIZES"] = str(args.higlass_chromsizes.resolve())
    if args.higlass:
        config["EXPORT_HIGLASS"] = True
    if args.no_ucsc:
        config["EXPORT_UCSC"] = False
    outputs = export_browser(config)
    audit = next(p for p in outputs if p.name == "export_audit.json")
    print(json.dumps({"status": "passed", "output_directory": str(audit.parent), "artifacts": len(outputs)}))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
