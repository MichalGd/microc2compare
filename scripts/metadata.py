#!/usr/bin/env python3
"""Standard-library CLI used before scientific environments are entered."""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from microc2compare.contracts import ContractError, STAGES, enabled, normalize, validate_metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("normalize", "plan", "samples", "get"))
    parser.add_argument("--config-dir")
    parser.add_argument("--config")
    parser.add_argument("--key")
    args = parser.parse_args()
    if args.command == "normalize":
        result = normalize(os.environ, args.config_dir)
        profile = os.environ.get("M2C_LOADED_SOFTWARE_PROFILE")
        if profile:
            path = Path(profile).resolve()
            checksum = hashlib.sha256(path.read_bytes()).hexdigest()
            if checksum != os.environ.get("M2C_LOADED_SOFTWARE_PROFILE_SHA256"):
                raise ContractError("Software profile changed while loading configuration")
            result["_software_profile"] = {"path": str(path), "sha256": checksum,
                                           "selection": os.environ["M2C_SOFTWARE_PROFILE_ORIGIN"]}
        print(json.dumps(result, sort_keys=True))
        return
    c = json.loads(Path(args.config).read_text())
    if args.command == "get":
        value = c[args.key]
        print(str(value).lower() if isinstance(value, bool) else value)
        return
    metadata = validate_metadata(c)
    if args.command == "samples":
        print("\n".join(s["sample_id"] for s in metadata["samples"]))
        return
    print(f"Run: {c['run_dir']}\nImported matrices: read-only counts + stored weights")
    print("Analysis selection: " + ",".join(c.get("_requested_analyses", ["configured full pipeline"])))
    print("Chromosome scope: "+c['CHROMOSOME_SCOPE']+(' '+c['SEX_CHROMOSOMES'] if c['SEX_CHROMOSOMES'] else ''))
    from microc2compare.resolutions import loop_resolutions, quantification_resolutions, chromosight_discovery_resolutions
    from microc2compare.analyses import has_loop_discovery, needs_loop_counts
    if c['RUN_LOOPS']:print('Loop discovery: '+c['LOOP_CALLER']+'; grids '+','.join(map(str,loop_resolutions(c)))+' bp')
    if c['RUN_CHROMOSIGHT_DISCOVERY']:print('Independent CPU loop discovery: chromosight; grids '+','.join(map(str,chromosight_discovery_resolutions(c)))+' bp')
    if has_loop_discovery(c) or needs_loop_counts(c):print('Loop measurements: '+','.join(map(str,quantification_resolutions(c)))+' bp')
    if c['RUN_STRIPES']:print('Stripe discovery: StripePy CPU; grids '+c['STRIPE_RESOLUTIONS']+' bp; common eligible segments')
    if c['RUN_BROWSER_EXPORT']:print('Browser exports: '+','.join(name for key,name in [('EXPORT_UCSC','UCSC'),('EXPORT_HIGLASS','HiGlass conversion')] if c[key]))
    if c["RUN_AGGREGATES"]: print("Aggregate modules: " + c["AGGREGATE_MODULES"])
    for i, stage in enumerate(STAGES, 1):
        mode = "enabled" if enabled(c, stage) else "skipped_config"
        print(f"{i:02d} {stage}: {mode}")
    print(f"Resource budget: {c['TOTAL_CPUS']} CPUs; {c['MAX_MEMORY_MB']} decimal MB")
    print("New persistent matrices: " + ("explicit optional derived products" if c["RUN_POOL"] and c["POOL_MODE"] == "create" else "none"))
    print("Depth sensitivity: " + (c["DEPTH_SENSITIVITY_MODE"] if c["RUN_DEPTH_SENSITIVITY"] else "disabled"))
    print("Imported-map balancing: forbidden")


if __name__ == "__main__":
    try:
        main()
    except (ContractError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
