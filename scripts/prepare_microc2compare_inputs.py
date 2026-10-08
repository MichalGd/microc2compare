#!/usr/bin/env python3
"""Preview microc2tracks outputs and optionally export reviewable input drafts."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from microc2compare.handoff import discover, summary, write_drafts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, help="Existing microc2tracks CSV/TSV")
    parser.add_argument("--upstream-results", required=True, help="Existing microc2tracks OUTDIR")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--preview", action="store_true", help="Read-only inventory (also the default)")
    mode.add_argument("--write-drafts", metavar="NEW_DIRECTORY", help="Export drafts outside upstream results; never overwrite")
    parser.add_argument("--project", default="project", help="Downstream project identifier")
    parser.add_argument("--resolutions", default="10000,100000", help="Comma-separated resolutions to inspect, in bp")
    parser.add_argument("--reference-id", help="Select one canonical reference/alias from a mixed run")
    parser.add_argument("--legacy-reference-id", help="Explicit identity ONLY for rows lacking sheet/run reference metadata")
    parser.add_argument("--json", action="store_true", help="Print the full machine-readable inventory")
    args = parser.parse_args(argv)
    try:
        resolutions = tuple(int(value.strip()) for value in args.resolutions.split(","))
        report = discover(args.samplesheet, args.upstream_results, args.project, resolutions,
                          args.reference_id, args.legacy_reference_id)
        if args.write_drafts:
            output = write_drafts(report, args.write_drafts)
            print(f"Drafts written to {output}; review READ_BEFORE_RUNNING.md before preflight.", file=sys.stderr)
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) if args.json else summary(report))
        return 0
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
