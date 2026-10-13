"""Exercise convenience commands through an installed, profile-bound launcher."""
import csv
import json
import subprocess
from pathlib import Path

import pytest

from microc2compare.contracts import file_hash, read_json, read_table
from test_software_profiles import installation, BASH, shell_path, module, load
from test_handoff import upstream, matrix, write_sheet


def command(installation, *args):
    return subprocess.run([str(BASH), shell_path(installation[3]), *map(str, args)],
                          env=installation[5], capture_output=True, text=True, timeout=120)


def test_init_installed_preset_and_override_without_software_settings(installation, tmp_path):
    conf = tmp_path / 'my private project' / 'project.conf'
    result = command(installation, '--init', shell_path(conf))
    assert result.returncode == 0, result.stderr
    assert 'PYTHON_BOOTSTRAP=' not in conf.read_text() and 'export ' not in conf.read_text()
    assert len(read_table(conf.parent / 'upstream_inputs.tsv')) == 0
    with conf.open('a') as handle:
        handle.write('\nLOOP_MIN_DISTANCE_BP="6000"\n')
    effective = load(installation, conf)
    assert effective.returncode == 0, effective.stderr
    c = json.loads(effective.stdout)
    assert c['LOOP_MIN_DISTANCE_BP'] == 6000
    assert c['LOOP_CALLER'] == 'mustache' and c['LOOP_QUANTIFICATION_RESOLUTIONS'] == '2000,4000'
    assert c['AGGREGATE_MODULES'] == 'apa,ada,saddle' and c['CHROMOSOME_SCOPE'] == 'autosomes'
    assert c['PYTHON_BOOTSTRAP'] != 'python3' and c['_software_profile']['selection'] == 'launcher_default'
    before = file_hash(conf)
    assert command(installation, '--init', shell_path(conf)).returncode != 0
    assert file_hash(conf) == before


def test_variant_preserves_relative_inputs_and_uses_literal_parameters(installation, tmp_path):
    conf = tmp_path / 'project.conf'
    assert command(installation, '--init', shell_path(conf)).returncode == 0
    before = file_hash(conf)
    result = command(installation, '-c', shell_path(conf), '--variant', 'tads_05', '--set', 'INSULATION_BOUNDARY_THRESHOLD=0.5')
    assert result.returncode == 0, result.stderr
    updated = list(tmp_path.glob('project.tads_05-*.conf'))[0]
    c = json.loads(load(installation, updated).stdout)
    assert c['RUN_ID'] == 'tads_05' and c['INSULATION_BOUNDARY_THRESHOLD'] == '0.5'
    assert c['SAMPLESHEET'] == str(tmp_path / 'samples.csv') and file_hash(conf) == before
    bad = command(installation, '-c', shell_path(conf), '--variant', 'bad', '--set', 'TYPO=1')
    assert bad.returncode != 0 and not list(tmp_path.glob('project.bad-*.conf'))


@pytest.mark.parametrize('args', [('--audit-resolutions', '--prepare-inputs'),
                                  ('--audit-resolutions', '--dry-run'),
                                  ('--validate', '--analysis', 'qc')])
def test_project_actions_reject_ambiguous_execution(installation, args):
    result = command(installation, *args)
    assert result.returncode != 0 and ('Choose one' in result.stderr or 'cannot be combined' in result.stderr)


def test_register_assets_once_generates_views_and_keeps_sources(installation, tmp_path):
    conf = tmp_path / 'project.conf'
    assert command(installation, '--init', shell_path(conf)).returncode == 0
    sizes = tmp_path / 'mm39.chrom.sizes'; sizes.write_text(''.join(f'chr{i}\t300000\n' for i in range(1, 20)) + 'chrX\t100000\nchrY\t100000\n')
    ref = conf.parent / 'references.json'
    ref.write_text(json.dumps({'reference_id': 'mm39', 'chromsizes': sizes.name}))
    before = {p: file_hash(p) for p in (conf, sizes, ref)}
    out = tmp_path / 'registered'
    result = command(installation, '-c', shell_path(conf), '--register-resources', '--out', shell_path(out))
    assert result.returncode == 0, result.stderr
    registered = read_json(out / 'references.json')
    assert registered['chromsizes_sha256'] == file_hash(sizes)
    assert file_hash(registered['views']) == registered['views_sha256']
    assert all(file_hash(p) == sha for p, sha in before.items())
    updated = list(tmp_path.glob('project.resources-*.conf'))
    assert len(updated) == 1
    effective = json.loads(load(installation, updated[0]).stdout)
    assert effective['REFERENCE_REGISTRY'] == str(out / 'references.json')
    assert 'phasing track' in result.stdout


def test_register_existing_fingerprint_mismatch_never_blesses_changed_asset(installation, tmp_path):
    conf = tmp_path / 'project.conf'
    assert command(installation, '--init', shell_path(conf)).returncode == 0
    sizes = tmp_path / 'sizes'; sizes.write_text('chr1\t300000\n')
    (tmp_path / 'references.json').write_text(json.dumps({'reference_id': 'mm39', 'chromsizes': str(sizes), 'chromsizes_sha256': '0' * 64}))
    out = tmp_path / 'registered'
    result = command(installation, '-c', shell_path(conf), '--register-resources', '--out', shell_path(out))
    assert result.returncode != 0 and 'Existing fingerprint mismatch' in result.stderr
    assert not out.exists()


def test_prepare_multiple_upstreams_derives_project_reference_and_grids(installation, upstream, tmp_path):
    root, sheet, _, merged, single, _ = upstream
    second = tmp_path / 'second upstream'; second.mkdir()
    sheet2 = tmp_path / 'second.csv'
    write_sheet(sheet2, [{'sample': 'T1', 'reference_genome': 'mm39', 'assay': 'microc', 'condition': 'treated', 'biological_replicate': '1', 'technical_replicate': '1'}])
    target = matrix(second / 'T1/04_matrices/T1.norm.mcool')
    conf = tmp_path / 'config' / 'project.conf'
    assert command(installation, '--init', shell_path(conf), '--preset', 'existing').returncode == 0
    with (conf.parent / 'upstream_inputs.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['samplesheet', 'upstream_results', 'individual_matrices'])
        writer.writerow([sheet, root, 'false']); writer.writerow([sheet2, second, 'true'])
    before = {p: file_hash(p) for p in (merged, single, target, sheet, sheet2)}
    out = tmp_path / 'combined'
    result = command(installation, '-c', shell_path(conf), '--prepare-inputs', '--out', shell_path(out))
    assert result.returncode == 0, result.stderr
    samples = read_table(out / 'samples.csv')
    assert len(samples) == 3 and {s['project'] for s in samples} == {'experiment'}
    assert len(read_table(out / 'technical_membership.tsv')) == 4
    assert all(Path(s['provenance_path']).is_absolute() for s in samples)
    assert all(read_json(s['provenance_path'])['draft'] for s in samples)
    assert before == {p: file_hash(p) for p in before}
    assert len(list(conf.parent.glob('project.inputs-*.conf'))) == 1


def test_installed_audit_and_validation_need_only_project_config(installation, tmp_path):
    data = tmp_path / 'private demo'
    starter = tmp_path / 'starter' / 'project.conf'
    assert command(installation, '--init', shell_path(starter)).returncode == 0
    generated = command(installation, '-c', shell_path(starter), '--demo', '--out', shell_path(data))
    assert generated.returncode == 0 and 'synthetic' in generated.stdout, generated.stderr
    conf = data / 'project.conf'
    out = tmp_path / 'audit'
    before = {p: file_hash(p) for p in data.glob('*.mcool')}
    result = command(installation, '-c', shell_path(conf), '--audit-resolutions', '--out', shell_path(out))
    assert result.returncode == 0, result.stderr
    audit = read_json(out / 'audit.json')
    assert audit['resolutions'] == [10000] and len(audit['rows']) == 6 and audit['read_only']
    result = command(installation, '-c', shell_path(conf), '--validate')
    assert result.returncode != 0, 'An unrun project must not validate successfully'
    # Reduced fixture budgets stay project variables, not software settings.
    with conf.open('a') as handle:
        handle.write('\nTOTAL_CPUS="2"\nMAX_MEMORY_MB="2000"\nMEMORY_MB_PREFLIGHT="1000"\nMEMORY_MB_QC="1000"\nMEMORY_MB_REPORT="1000"\n')
    result = command(installation, '-c', shell_path(conf), '--analysis', 'qc')
    assert result.returncode == 0, result.stderr
    result = command(installation, '-c', shell_path(conf), '--validate')
    assert result.returncode == 0 and 'qc/C1: valid' in result.stdout, result.stderr
    # Changed scientific settings must be reported stale, even for a previously good run.
    with conf.open('a') as handle:
        handle.write('\nLOOP_MIN_DISTANCE_BP="30000"\n')
    result = command(installation, '-c', shell_path(conf), '--validate')
    assert result.returncode != 0 and 'missing_or_stale' in result.stdout
    assert before == {p: file_hash(p) for p in before}
    assert not list((data / 'results').rglob('.run-lock'))


def test_prepare_grids_reuses_registered_products_without_rebalancing(installation, tmp_path):
    import cooler
    import h5py
    import numpy as np
    import pandas as pd
    from microc2compare.contracts import atomic_json
    source = tmp_path / 'source.cool'
    bins = cooler.binnify(pd.Series({'chr1': 80000}), 1000)
    i, j = np.triu_indices(len(bins))
    cooler.create_cooler(str(source), bins, pd.DataFrame({'bin1_id': i, 'bin2_id': j, 'count': np.full(len(i), 10, dtype=np.int64)}), ordered=True, assembly='mm39')
    recipe = {'cis_only': False, 'ignore_diags': 0, 'mad_max': 5, 'min_nnz': 1, 'min_count': 0, 'tol': 1e-5}
    with h5py.File(source, 'r+') as handle:
        ds = handle['bins'].create_dataset('weight', data=np.ones(len(bins)))
        ds.attrs.update(recipe); ds.attrs['converged'] = True
    provenance = tmp_path / 'source.json'
    atomic_json(provenance, {'reference_id': 'mm39', 'balancing': {'1000': {'convention': 'multiplicative', 'converged': True, 'recipe': recipe}}})
    sheet = tmp_path / 'samples.tsv'
    pd.DataFrame([{'sample_id': 'one', 'reference_id': 'mm39', 'matrix_path': str(source), 'provenance_path': str(provenance)}]).to_csv(sheet, sep='\t', index=False)
    conf = tmp_path / 'project.conf'
    conf.write_text('REFERENCE_ID="mm39"\nSAMPLESHEET="samples.tsv"\nPS_RESOLUTION="2000"\nHICREP_RESOLUTION="2000"\n'
                    'RUN_COMPARTMENTS="false"\nRUN_INSULATION="false"\nRUN_LOOPS="false"\nRUN_AGGREGATES="false"\nRUN_DIFFERENTIAL="false"\n'
                    'THREADS_BALANCE="1"\nMEMORY_MB_BALANCE="1000"\n')
    before = file_hash(source)
    out = tmp_path / 'grids'
    result = command(installation, '-c', shell_path(conf), '--prepare-resolutions', '--out', shell_path(out))
    assert result.returncode == 0, result.stdout + result.stderr
    derived = out / 'one.2000.cool'
    assert cooler.Cooler(str(derived)).info['sum'] == cooler.Cooler(str(source)).info['sum']
    assert file_hash(source) == before
    updated = list(tmp_path.glob('project.resolutions-*.conf'))[0]
    original_derived = file_hash(derived)
    again = tmp_path / 'reused'
    result = command(installation, '-c', shell_path(updated), '--prepare-resolutions', '--out', shell_path(again))
    assert result.returncode == 0, result.stdout + result.stderr
    assert not list(again.glob('*.cool'))
    assert read_json(again / 'matrix_registry.json')['samples']['one']['2000'] == str(derived)
    assert file_hash(derived) == original_derived and file_hash(source) == before


def test_annotation_registration_requires_source_and_assembly_metadata(installation, tmp_path):
    conf = tmp_path / 'project.conf'
    assert command(installation, '--init', shell_path(conf)).returncode == 0
    sizes = tmp_path / 'sizes'; sizes.write_text(''.join(f'chr{i}\t300000\n' for i in range(1, 20)))
    (tmp_path / 'references.json').write_text(json.dumps({'reference_id': 'mm39', 'chromsizes': str(sizes)}))
    gtf = tmp_path / 'annotation.gtf'; gtf.write_text('# fixture only\n')
    ann = tmp_path / 'annotation.json'; ann.write_text(json.dumps({'reference_id': 'mm39', 'gtf': gtf.name}))
    with conf.open('a') as handle: handle.write('\nANNOTATION_REGISTRY="annotation.json"\n')
    out = tmp_path / 'registered'
    failed = command(installation, '-c', shell_path(conf), '--register-resources', '--out', shell_path(out))
    assert failed.returncode != 0 and 'GTF path and release' in failed.stderr
    ann.write_text(json.dumps({'reference_id': 'mm39', 'gtf': gtf.name, 'gtf_release': 'synthetic fixture'}))
    result = command(installation, '-c', shell_path(conf), '--register-resources', '--out', shell_path(out))
    assert result.returncode == 0, result.stderr
    assert read_json(out / 'annotation.json')['gtf_sha256'] == file_hash(gtf)
