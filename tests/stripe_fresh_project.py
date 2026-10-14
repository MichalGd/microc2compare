"""Generate a fresh uneven-replication stripe project for actual Linux CI.

This original synthetic fixture tests software integration, not biological
calibration. Its two untreated cultures share a donor and remain descriptive.
"""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))

import cooler
import h5py
import numpy as np
import pandas as pd

from microc2compare.contracts import atomic_json, file_hash, SAMPLE_FIELDS, write_table


def generate(destination, core_prefix, stripe_prefix, conda_exe):
    root = Path(destination).resolve()
    if root.exists() and any(root.iterdir()):
        raise ValueError("Fresh stripe fixture requires an empty destination")
    root.mkdir(parents=True, exist_ok=True)
    bins = cooler.binnify(pd.Series({"freshChr": 640000}), 4000)
    i, j = np.triu_indices(len(bins))
    counts = np.maximum(1, 200//(j-i+1)).astype(np.int64)
    counts[(i == 25) & (j > 25)] += 100
    samples, membership = [], []
    for index, (sid, condition, subject) in enumerate([("newA", "untreated", "sharedDonor"),
            ("newB", "untreated", "sharedDonor"), ("newC", "perturbed", "otherDonor")]):
        matrix, provenance = root/f"{sid}.cool", root/f"{sid}.provenance.json"
        cooler.create_cooler(str(matrix), bins, pd.DataFrame(dict(bin1_id=i, bin2_id=j, count=counts)),
                             ordered=True, assembly="freshAssembly")
        with h5py.File(matrix, "r+") as handle:
            weights = np.ones(len(bins)); weights[80] = np.nan
            handle["bins"].create_dataset("weight", data=weights)
        atomic_json(provenance, dict(reference_id="freshAssembly", count_units="unique_contacts",
            filters=dict(mapq=30, min_distance_bp=0, chromosome_policy="synthetic_freshChr"),
            balancing={"4000": dict(weight_column="weight", convention="multiplicative", converged=True,
                                    recipe={"kind": "declared_synthetic_weights"})},
            dedup_scope="joint_library", library_ids=[sid+"_library"], source_ids=[sid+"_source"]))
        samples.append(dict(sample_id=sid, project="freshStripe", condition=condition, assay="microc",
            biological_replicate=str(index+1), subject_id=subject, experimental_unit_id=sid+"_culture",
            reference_id="freshAssembly", matrix_path=str(matrix), provenance_path=str(provenance)))
        membership.append(dict(sample_id=sid, library_id=sid+"_library", source_id=sid+"_source",
            source_kind="independent_library", dedup_group_id=sid+"_dedup", dedup_scope="joint_library",
            source_checksum=file_hash(matrix)))
    write_table(root/"samples.tsv", samples, list(SAMPLE_FIELDS)+["experimental_unit_id"])
    write_table(root/"technical_membership.tsv", membership, list(membership[0]))
    (root/"chromsizes.tsv").write_text("freshChr\t640000\n", encoding="utf-8")
    (root/"views.bed").write_text("freshChr\t0\t640000\tselected\n", encoding="utf-8")
    atomic_json(root/"reference.json", dict(reference_id="freshAssembly", chromsizes="chromsizes.tsv",
        chromsizes_sha256=file_hash(root/"chromsizes.tsv"), views="views.bed", views_sha256=file_hash(root/"views.bed"),
        chromosome_classes={"autosomes": ["freshChr"]}))
    settings = dict(PROJECT="freshStripe", ANALYSIS_SET_ID="newSet", RUN_ID="newRun", REFERENCE_ID="freshAssembly",
        SAMPLESHEET="samples.tsv", TECHNICAL_MEMBERSHIP="technical_membership.tsv", REFERENCE_REGISTRY="reference.json",
        OUTDIR="results", SCRATCH_ROOT="scratch", ENVIRONMENT_MODE="reuse", CONDA_EXE=conda_exe,
        PYTHON_BOOTSTRAP=sys.executable, CORE_ENV_PREFIX=core_prefix, STRIPE_ENV_PREFIX=stripe_prefix,
        RUN_STRIPES="true", STRIPE_RESOLUTIONS="4000", THREADS_STRIPES="1", MEMORY_MB_STRIPES="1000",
        TOTAL_CPUS="2", MAX_MEMORY_MB="3000", MIN_FREE_DISK_MB="100", RUN_FINAL_REPORT="true",
        THREADS_PREFLIGHT="1", MEMORY_MB_PREFLIGHT="1000", THREADS_REPORT="1", MEMORY_MB_REPORT="1000",
        THREADS_BROWSER_EXPORT="1", MEMORY_MB_BROWSER_EXPORT="1000")
    path = root/"project.conf"
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for key, value in settings.items():
            if any(char in str(value) for char in ('"', "\n", "$", "`")):
                raise ValueError(f"Unsafe shell configuration value for {key}")
            handle.write(f'{key}="{value}"\n')
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("--core-prefix", required=True)
    parser.add_argument("--stripe-prefix", required=True)
    parser.add_argument("--conda-exe", required=True)
    args = parser.parse_args()
    print(generate(args.out, args.core_prefix, args.stripe_prefix, args.conda_exe))
