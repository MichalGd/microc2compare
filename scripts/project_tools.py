#!/usr/bin/env python3
"""Project actions behind the shared command; software and runtime stay internal."""
import argparse
import csv
import importlib.util
import os
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from microc2compare.contracts import (atomic_json, file_hash, load_structured,
                                    read_json, read_table, require, resolve, write_table)


def script(name):
    spec = importlib.util.spec_from_file_location('project_' + name, ROOT / 'scripts' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def token():
    return datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S') + '-' + uuid4().hex[:8]


def initialize(path, preset, reference_id='mm39'):
    from microc2compare.contracts import IDENTIFIER
    require(preset in ('existing', 'practical_2_4kb', 'hsieh_8grid', 'short_range_4kb'), 'Unknown project preset')
    require(IDENTIFIER.fullmatch(reference_id), 'Supply a safe reference assembly identifier')
    path = Path(path).resolve()
    require(not path.exists(), 'Project configuration exists; refusing overwrite')
    # No fictitious samples, biological metadata, reference identities or hashes.
    companions = [path.parent / name for name in ('upstream_inputs.tsv', 'comparisons.tsv', 'references.json')]
    require(not any(p.exists() for p in companions), 'Companion files exist; choose a new configuration directory')
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f'''# Edit experiment-specific values. Software and runtime are supplied internally.
# This preset line supplies defaults; assignments below can override them.
use_analysis_profile {preset}
PROJECT="experiment"
ANALYSIS_SET_ID="cell_type_comparison"
RUN_ID="run001"
REFERENCE_ID="{reference_id}"
SAMPLESHEET="samples.csv"
TECHNICAL_MEMBERSHIP="technical_membership.tsv"
COMPARISONS="comparisons.tsv"
REFERENCE_REGISTRY="references.json"
# Only needed when importing new microc2tracks outputs:
UPSTREAM_INPUTS="upstream_inputs.tsv"
OUTDIR="../results"
SCRATCH_ROOT="../scratch"
TOTAL_CPUS="24"
MAX_MEMORY_MB="160000"
MAX_PARALLEL_JOBS="8"
MAX_PARALLEL_SAMPLES="8"
MAX_PARALLEL_MATRIX="8"
RUN_BROWSER_EXPORT="true"
# Optional annotations need genuine assembly/source-release information:
RUN_ANNOTATION="false"
ANNOTATION_REGISTRY=""
'''
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(text)
    write_table(companions[0], [], ['samplesheet', 'upstream_results', 'individual_matrices'])
    write_table(companions[1], [], ['comparison_id', 'case_condition', 'reference_condition', 'design', 'inference', 'modules'])
    atomic_json(companions[2], {'reference_id': reference_id, 'chromsizes': '', 'phasing_track': ''})
    print(f'Created {path}\nEdit the project identity, input paths and reference assets before preparation.')
    print(f'Next: microc2compare -c {shlex.quote(str(path))} --prepare-inputs')


def destination(c, action, out):
    # Every action writes to a NEW project-owned directory; no general shared data area.
    return Path(out).resolve() if out else Path(c['OUTDIR']) / c['PROJECT'] / 'preparation' / (action + '-' + token())


def variant(config, changes, label):
    config = Path(config).resolve()
    target = config.with_name(config.stem + '.' + label + '-' + token() + config.suffix)
    text = config.read_text(encoding='utf-8-sig') + '\n# Prepared inputs; original configuration preserved.\n'
    for key, value in changes.items():
        require('\n' not in str(value) and '\r' not in str(value), 'Configuration values must be single-line literals')
        text += key + '=' + shlex.quote(str(value)) + '\n'
    with target.open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(text)
    print(f'Prepared configuration: {target}')
    print(f'Next: microc2compare -c {shlex.quote(str(target))} --dry-run')
    return target


def prepare_inputs(c, config, out):
    from microc2compare.handoff import discover, summary, write_drafts
    from microc2compare.analyses import resolution_set
    require(c['UPSTREAM_INPUTS'], 'Set UPSTREAM_INPUTS to a table of upstream samplesheets/results, or use your completed SAMPLESHEET directly')
    table = Path(c['UPSTREAM_INPUTS'])
    rows = read_table(table, ('samplesheet', 'upstream_results'))
    require(rows, 'UPSTREAM_INPUTS has no rows; supply your experiment-specific input paths')
    out = Path(out)
    require(not out.exists(), 'Choose a NEW preparation directory')
    reports = []
    for row in rows:
        mode = row.get('individual_matrices', '').lower() or 'false'
        require(mode in ('true', 'false'), 'individual_matrices must be true or false')
        report = discover(resolve(row['samplesheet'], table.parent), resolve(row['upstream_results'], table.parent),
                          c['PROJECT'], tuple(resolution_set(c)), c['REFERENCE_ID'], c['REFERENCE_ID'], mode == 'true')
        require(report['units'], 'No biological matrices discovered; inspect the upstream sample descriptions')
        reports.append(report)
    # Preserve each upstream inventory and its evidence gaps; do not invent provenance.
    samples, members = [], []
    for index, report in enumerate(reports, 1):
        base = out / f'upstream-{index:03d}'
        write_drafts(report, base)
        print(summary(report))
        units = read_table(base / 'samples.csv')
        for sample in units:
            for key in ('matrix_path', 'provenance_path'):
                sample[key] = resolve(sample[key], base)
        samples.extend(units)
        members.extend(read_table(base / 'technical_membership.tsv'))
    require(len({r['sample_id'] for r in samples}) == len(samples), 'Repeated sample IDs across upstream inputs; correct the source descriptions')
    require(len({r['matrix_path'] for r in samples}) == len(samples), 'A matrix was assigned to more than one biological observation')
    require(len({r['source_id'] for r in members}) == len(members), 'Repeated technical source IDs across upstream inputs')
    with (out / 'samples.csv').open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(samples[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(samples)
    write_table(out / 'technical_membership.tsv', members, list(members[0]))
    (out / 'READ_BEFORE_RUNNING.md').write_text(
        '# Draft inputs: review required\n\nEach upstream directory contains its inventory and READ_BEFORE_RUNNING.md. '
        'Resolve missing library, deduplication, filtering and processing evidence in the draft provenance files. '
        'Biological and technical identities come from the samplesheets; this helper does not establish independence. '
        'Historical reference fields missing upstream use the explicit REFERENCE_ID from project.conf and retain evidence gaps. '
        'Do not remove draft markers without completing the recorded evidence. Preflight remains mandatory.\n', encoding='utf-8')
    print(f'Drafts require review: {out / "READ_BEFORE_RUNNING.md"}')
    return variant(config, {'SAMPLESHEET': out / 'samples.csv', 'TECHNICAL_MEMBERSHIP': out / 'technical_membership.tsv'}, 'inputs')


def register_resources(c, config, out):
    from microc2compare.matrix import reference
    out = Path(out)
    require(not out.exists(), 'Choose a NEW resource registration directory')
    registry = Path(c['REFERENCE_REGISTRY'])
    ref = load_structured(registry)
    require(ref.get('reference_id') == c['REFERENCE_ID'], 'Reference registry identity mismatch')
    require(ref.get('chromsizes'), 'Supply the actual chromsizes path in REFERENCE_REGISTRY')
    changes = {}

    def asset(data, key, base):
        if not data.get(key):
            return
        path = resolve(data[key], base)
        require(Path(path).is_file(), f'Missing resource {key}: {path}')
        sha = file_hash(path)
        expected = data.get(key + '_sha256')
        require(not expected or expected == sha, f'Existing fingerprint mismatch: {key}; inspect changed assets before registering')
        data[key], data[key + '_sha256'] = path, sha

    for key in ('chromsizes', 'views', 'blacklist', 'phasing_track'):
        asset(ref, key, registry.parent)
    annotation = None
    if c['ANNOTATION_REGISTRY']:
        annotation_path = Path(c['ANNOTATION_REGISTRY'])
        annotation = load_structured(annotation_path)
        require(annotation.get('reference_id') == c['REFERENCE_ID'], 'Annotation assembly differs from the project')
        require(annotation.get('gtf') and annotation.get('gtf_release'), 'Annotation needs the genuine GTF path and release')
        require(not annotation.get('ccre') or annotation.get('ccre_source'), 'cCRE annotation needs its source/release and any conversion history')
        for key in ('gtf', 'ccre'):
            asset(annotation, key, annotation_path.parent)
    out.mkdir(parents=True)
    if not ref.get('views'):
        sizes = read_sizes(ref['chromsizes'])
        views = out / 'whole_chromosome.views.bed'
        views.write_text(''.join(f'{chrom}\t0\t{length}\t{chrom}\n' for chrom, length in sizes), encoding='utf-8')
        ref['views'] = str(views.resolve()); ref['views_sha256'] = file_hash(views)
    registered = out / 'references.json'
    atomic_json(registered, ref)
    # Validate coordinates/scope without claiming the scientific phasing gate passed.
    reference(dict(c, REFERENCE_REGISTRY=str(registered), RUN_COMPARTMENTS=False))
    changes['REFERENCE_REGISTRY'] = registered
    if annotation is not None:
        path = out / 'annotation.json'; atomic_json(path, annotation)
        changes['ANNOTATION_REGISTRY'] = path
    atomic_json(out / 'registration.json', {'reference_id': c['REFERENCE_ID'], 'assets_modified': False,
                                          'statistical_calibration': False, 'phasing_track_present': bool(ref.get('phasing_track'))})
    print('Resource paths and SHA256 identities registered; source assets were not copied or modified.')
    if not ref.get('phasing_track'):
        print('Full compartments/saddle analysis still needs a genuine common phasing track; QC can run without it.')
    return variant(config, changes, 'resources')


def read_sizes(path):
    rows = []
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split()
        require(len(parts) == 2 and parts[1].isdigit() and int(parts[1]) > 0, 'Invalid chromosome sizes')
        rows.append((parts[0], int(parts[1])))
    require(rows and len({r[0] for r in rows}) == len(rows), 'Missing/duplicate chromosomes')
    return rows


def run_action(action, c, config, out=None, run_id=None, overrides=()):
    if action == 'export-browser':
        from microc2compare.browser_export import export_browser
        require(not (Path(c['run_dir']) / '.run-lock').exists(), 'Finish the active run before exporting its saved results')
        settings = dict(c)
        if out: settings['BROWSER_OUTPUT_DIR'] = str(Path(out).resolve())
        products = export_browser(settings)
        print('Browser exports: ' + str(next(p.parent for p in products if p.name == 'README.md')))
        return 0
    if action == 'input-capabilities':
        from microc2compare.contracts import validate_metadata
        from microc2compare.resolutions import input_capabilities
        out = destination(c, action, out)
        out.mkdir(parents=True, exist_ok=False)
        report = input_capabilities(c, validate_metadata(c)['samples'])
        atomic_json(out / 'capabilities.json', report)
        print(f'Input grid capabilities: {out / "capabilities.json"}')
        return 0
    if action == 'variant':
        from microc2compare.contracts import IDENTIFIER
        require(run_id and IDENTIFIER.fullmatch(run_id) and run_id != c['RUN_ID'], 'Supply a new safe run ID')
        allowed = read_json(ROOT / 'schemas/config.schema.json')['properties']
        changes = {'RUN_ID': run_id}
        for setting in overrides:
            key, separator, value = setting.partition('=')
            require(separator and key in allowed and key != 'RUN_ID', 'Use --set KNOWN_PARAMETER=value; the run ID comes from --variant')
            changes[key] = value
        variant(config, changes, run_id)
        return 0
    if action == 'validate':
        from microc2compare.analyses import configure_selection
        saved = read_json(Path(c['run_dir']) / '00_metadata/config.json')
        # Validate the modules actually run, with TODAY'S project values and source.
        configure_selection(c, ','.join(saved.get('_requested_analyses', [])))
        print(f'Checking recorded run: {c["run_dir"]}')
        return script('validate_outputs').validate(c)
    out = destination(c, action, out)
    if action == 'prepare-inputs':
        prepare_inputs(c, config, out)
    elif action == 'register-resources':
        register_resources(c, config, out)
    elif action == 'audit-resolutions':
        from microc2compare.analyses import resolution_set
        script('audit_resolutions').audit(c['SAMPLESHEET'], c['REFERENCE_REGISTRY'], out, resolution_set(c),
                                         c['MATRIX_REGISTRY'], c['CHROMOSOME_SCOPE'], c['SEX_CHROMOSOMES'], c['ANALYSIS_CHROMOSOMES'])
        evidence = read_json(out/'audit.json')
        fine = [r for r in resolution_set(c) if r <= 1000]
        if fine:
            atomic_json(out/'eligibility.template.json', {
                'reference_id': c['REFERENCE_ID'], 'scope_digest': evidence['scope']['digest'],
                'audit_path': 'audit.json', 'audit_sha256': file_hash(out/'audit.json'),
                'resolutions': {str(r): {'approved': False, 'reviewer': '', 'retention_review': '',
                                         'reproducibility_review': ''} for r in fine},
                'status': 'template_only_not_scientific_admission'})
    elif action == 'prepare-resolutions':
        from microc2compare.analyses import resolution_set
        from microc2compare.contracts import allocation_limits
        limits = allocation_limits()
        threads = min(c['THREADS_BALANCE'], c['TOTAL_CPUS'], max(1, int(limits['effective_cpus'])))
        require(c['MEMORY_MB_BALANCE'] <= c['MAX_MEMORY_MB'], 'MEMORY_MB_BALANCE exceeds MAX_MEMORY_MB')
        require(limits['available_memory_bytes'] is None or c['MEMORY_MB_BALANCE'] * 1000000 <= limits['available_memory_bytes'], 'Insufficient available memory for the configured balancing reservation')
        script('prepare_resolutions').prepare(c['SAMPLESHEET'], out, resolution_set(c), threads=threads, matrix_registry=c['MATRIX_REGISTRY'])
        variant(config, {'SAMPLESHEET': out / 'samples.tsv', 'MATRIX_REGISTRY': out / 'matrix_registry.json'}, 'resolutions')
    elif action == 'prepare-from-pairs':
        from microc2compare.pairs_preparation import prepare_from_pairs
        from microc2compare.analyses import resolution_set
        from microc2compare.contracts import allocation_limits
        limits = allocation_limits()
        threads = min(c['THREADS_BALANCE'], c['TOTAL_CPUS'], max(1, int(limits['effective_cpus'])))
        require(c['MEMORY_MB_BALANCE'] <= c['MAX_MEMORY_MB'], 'MEMORY_MB_BALANCE exceeds MAX_MEMORY_MB')
        require(limits['available_memory_bytes'] is None or c['MEMORY_MB_BALANCE']*1000000 <= limits['available_memory_bytes'], 'Insufficient available memory for the configured balancing reservation')
        prepare_from_pairs(dict(c, THREADS_BALANCE=threads), out, resolution_set(c))
        variant(config, {'SAMPLESHEET': out / 'samples.tsv', 'MATRIX_REGISTRY': out / 'matrix_registry.json',
                         'WEIGHT_COLUMN': 'weight', 'WEIGHT_CONVENTION': 'multiplicative'}, 'pairs-grids')
    elif action == 'demo':
        conf = script('generate_demo').generate(out, raw_only=True)
        print('Generated six synthetic biological-sample maps; these are software fixtures, not real biological data.')
        print(f'Next: microc2compare -c {shlex.quote(str(conf))} --analysis qc,hicrep')
    else:
        raise ValueError('Unknown project action')
    print(f'Preparation/check output: {out}')
    return 0


def main(argv=None):
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--action', required=True, choices=['init', 'prepare-inputs', 'register-resources', 'audit-resolutions', 'prepare-resolutions', 'prepare-from-pairs', 'input-capabilities', 'export-browser', 'validate', 'demo', 'variant'])
    parser.add_argument('--config-json'); parser.add_argument('--project-config', required=True)
    parser.add_argument('--out'); parser.add_argument('--preset', default='practical_2_4kb')
    parser.add_argument('--reference', default='mm39')
    parser.add_argument('--run-id'); parser.add_argument('--set', action='append', default=[])
    args = parser.parse_args(argv)
    if args.action == 'init':
        initialize(args.project_config, args.preset, args.reference)
        return 0
    require(args.config_json, 'Normalized project configuration required')
    return run_action(args.action, read_json(args.config_json), args.project_config, args.out, args.run_id, args.set)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(2)
