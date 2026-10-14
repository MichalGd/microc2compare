"""Coordinate alignment against real raw Cooler maps, including X/Y and scaffolds."""
import itertools
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import ContractError, atomic_json, file_hash, read_json
from microc2compare.features import LOOP_COLUMNS, quantify
from microc2compare.matrix import (aligned_pixel_chunks, domain, expected_table,
                                 native_indices, raw_summary, tool_views)
from microc2compare import science, qc_parallel


@pytest.fixture
def aligned(tmp_path, config, monkeypatch):
    sizes = tmp_path / 'sizes.tsv'
    sizes.write_text('chr1\t80000\nchrX\t80000\nchrY\t80000\nscaffold\t20000\n')
    views = tmp_path / 'views.bed'
    views.write_text(''.join(f'{chrom}\t0\t80000\t{chrom}\n' for chrom in ('chr1', 'chrX', 'chrY')))
    reference = tmp_path / 'reference.json'
    atomic_json(reference, dict(reference_id='mm39', chromosome_classes={'autosomes':['chr1'],'sex':['chrX','chrY'],'other':['scaffold']}, chromsizes=str(sizes), chromsizes_sha256=file_hash(sizes),
                               views=str(views), views_sha256=file_hash(views)))
    samples = []
    for sid, chromosomes in [('A', ['scaffold', 'chrY', 'chr1', 'chrX']), ('B', ['chr1', 'chrX', 'chrY'])]:
        lengths = pd.Series({chrom: 20000 if chrom == 'scaffold' else 80000 for chrom in chromosomes})
        bins = cooler.binnify(lengths, 10000)
        rows = []
        for i, j in itertools.combinations_with_replacement(range(len(bins)), 2):
            a, b = bins.iloc[i], bins.iloc[j]
            # Identical selected raw counts by coordinates, independent of native IDs.
            count = 100000 if 'scaffold' in (a.chrom, b.chrom) else (int(a.start + b.start) // 10000 + 1)
            rows.append((i, j, count))
        path = tmp_path / f'{sid}.cool'
        cooler.create_cooler(str(path), bins, pd.DataFrame(rows, columns=['bin1_id','bin2_id','count']), ordered=True)
        values = (bins.start.to_numpy() / 10000 + 1) / 10
        values[(bins.chrom.astype(str) == 'chrY') & (bins.start == 70000)] = np.nan
        if sid == 'B': values[(bins.chrom.astype(str) == 'chrX') & (bins.start == 20000)] = np.nan
        with h5py.File(path, 'r+') as h:
            ds = h['bins'].create_dataset('weight', data=values); ds.attrs['converged'] = True
        prov = tmp_path / f'{sid}.json'
        atomic_json(prov, dict(reference_id='mm39', count_units='unique_contacts',
                               filters=dict(mapq=30, min_distance_bp=0, chromosome_policy='all' if sid == 'A' else 'canonical'),
                               balancing={'10000':dict(convention='multiplicative', converged=True, recipe={'ignore_diags':2})}))
        samples.append(dict(sample_id=sid, matrix_path=str(path), provenance_path=str(prov), condition=sid))
    config.update(CHROMOSOME_SCOPE='autosomes_plus_sex', SEX_CHROMOSOMES='chrX,chrY', REFERENCE_REGISTRY=str(reference), REFERENCE_ID='mm39', run_dir=str(tmp_path/'results'),
                  SCRATCH_ROOT=str(tmp_path/'scratch'), LOOP_MIN_DISTANCE_BP=10000, LOOP_MAX_DISTANCE_BP=70000,
                  PS_RESOLUTION=10000, INSULATION_RESOLUTION=10000, LOOP_RESOLUTION=10000,
                  COMPARTMENT_RESOLUTION=10000, HICREP_RESOLUTION=10000, RUN_COMPARTMENTS=False,
                  RUN_INSULATION=False, RUN_LOOPS=False, RUN_AGGREGATES=False, RUN_DIFFERENTIAL=False,
                  PIXEL_CHUNKSIZE=17, INSULATION_IGNORE_DIAGS=0)
    monkeypatch.setattr(science, 'validate_metadata', lambda c: {'samples':samples})
    monkeypatch.setattr(qc_parallel, 'validate_metadata', lambda c: {'samples':samples})
    science.register(config, 'run')
    return config, samples


def test_alignment_retains_xy_common_masks_exact_counts_and_weights(aligned, monkeypatch):
    c, samples = aligned
    before = [file_hash(s['matrix_path']) for s in samples]
    def forbidden(*a, **kw): raise AssertionError('Default comparison attempted matrix creation/balancing')
    for name in ('create_cooler', 'balance_cooler', 'merge_coolers', 'zoomify_cooler'):
        monkeypatch.setattr(cooler, name, forbidden)
    bins, common, labels, views, weights = domain(c, samples, 10000)
    assert list(dict.fromkeys(bins.chrom)) == ['chr1', 'chrX', 'chrY']
    assert len(bins) == 24 and common.sum() == 22
    np.testing.assert_allclose(weights['A'], weights['B'], equal_nan=True)
    features = pd.DataFrame([[chrom, 0, 10000, chrom, 30000, 40000, chrom] for chrom in ('chr1','chrX','chrY')], columns=LOOP_COLUMNS)
    expected = []; exposures = []; summaries = []
    for sample in samples:
        clr, native, mask, local_labels, local_views, values = science.context(c, sample['sample_id'], 10000)
        ids = native_indices(clr, bins)
        np.testing.assert_array_equal(mask[ids], common)
        np.testing.assert_allclose(values[ids], weights[sample['sample_id']], equal_nan=True)
        for chrom in ('chr1', 'chrX', 'chrY'):
            raw = clr.matrix(balance=False).fetch((chrom,0,40000))
            balanced = clr.matrix(balance=True).fetch((chrom,0,40000))
            local = clr.bins()['weight'].fetch((chrom,0,40000)).to_numpy()
            np.testing.assert_allclose(balanced, raw * local[:,None] * local[None,:], equal_nan=True)
        counts, exposure, reasons = quantify(clr, features, native, mask, local_labels, c)
        assert list(counts) == [4,4,4] and not reasons
        exposures.append(exposure)
        expected.append(expected_table(clr,native,mask,local_labels,local_views,values,c))
        summaries.append(raw_summary(clr,native,mask,local_labels,c))
        science.qc(c,sample['sample_id'])
    assert exposures[0] == exposures[1] > 0
    pd.testing.assert_frame_equal(expected[0], expected[1])
    assert summaries[0]['raw_total'] > summaries[1]['raw_total']  # honest native scope
    assert summaries[0]['excluded_bin_fraction'] == summaries[1]['excluded_bin_fraction']
    assert [file_hash(s['matrix_path']) for s in samples] == before
    assert not list(Path(c['run_dir']).rglob('*.cool'))
    audit = read_json(Path(c['run_dir'])/'00_metadata/masks/10000.scope.json')
    assert audit['samples'][0]['excluded_chromosomes'] == ['scaffold']
    assert audit['analysis_chromosomes'][-2:] == ['chrX','chrY']
    first=cooler.Cooler(samples[0]['matrix_path'])
    assert list(tool_views(first,views).name)==['chrY','chr1','chrX']


@pytest.mark.parametrize('change', ['missing_y','different_length','different_grid','failed_weight','mapq','recipe'])
def test_selected_mismatches_are_not_silently_accepted(aligned, change):
    c, samples = aligned
    if change in ('failed_weight','mapq','recipe'):
        path = samples[1]['provenance_path']; prov = read_json(path)
        if change == 'failed_weight': prov['balancing']['10000']['converged'] = False
        if change == 'mapq': prov['filters']['mapq'] = 10
        if change == 'recipe': prov['balancing']['10000']['recipe']['ignore_diags'] = 3
        atomic_json(path,prov)
    else:
        with h5py.File(samples[1]['matrix_path'],'r+') as h:
            if change == 'missing_y': h['chroms/name'][2] = b'other'
            if change == 'different_length':
                h['chroms/length'][2] = 79999; h['bins/end'][-1] = 79999
            if change == 'different_grid': h['bins/start'][17] += 1
    with pytest.raises(ContractError): domain(c,samples,10000)


def test_new_selected_pool_pixel_remapping_preserves_cis_and_trans(aligned, tmp_path):
    c,samples=aligned
    bins,*_=domain(c,samples,10000)
    results=[]
    for sample in samples:
        native=cooler.Cooler(sample['matrix_path'])
        path=tmp_path/(sample['sample_id']+'.selected.cool')
        cooler.create_cooler(str(path),bins,aligned_pixel_chunks(native,bins,17),ordered=False,mergebuf=17,dtypes={'count':'int64'})
        results.append(cooler.Cooler(str(path)).matrix(balance=False)[:,:])
    np.testing.assert_array_equal(*results)
    assert results[0][0,8] > 0 and results[0][0,16] > 0


def test_raw_thinning_uses_common_grid_including_xy(aligned, monkeypatch):
    c,samples=aligned
    c.update(RUN_DEPTH_SENSITIVITY=True, RUN_HICREP=False, DOWNSAMPLING_REPEATS=1,
             DEPTH_SENSITIVITY_MODE='raw_only', KEEP_DOWNSAMPLED_MAPS=True)
    def forbidden(*a, **kw): raise AssertionError('Raw-only thinning tried balancing')
    monkeypatch.setattr(cooler,'balance_cooler',forbidden)
    science.depth_sensitivity(c,'run')
    for sample in samples:
        path=Path(c['run_dir'])/f"01_qc/depth_sensitivity/seed_0/{sample['sample_id']}/10000.cool"
        assert cooler.Cooler(str(path)).chromnames == ['chr1','chrX','chrY']


def test_hicrep_shared_xy_views_with_different_native_indices(aligned):
    pytest.importorskip('hicrep')
    c,_=aligned
    c.update(HICREP_SMOOTHING_BINS=1,HICREP_MAX_DISTANCE_BP=70000)
    products=science.hicrep(c,'run')
    table=pd.read_csv(products[0],sep='\t')
    assert set(table.region) == {'chr1','chrX','chrY'}
    assert np.allclose(table.scc,1)


def test_selected_pool_stage_conserves_coordinate_counts(aligned, monkeypatch):
    c,samples=aligned
    c.update(RUN_POOL=True,POOL_MODE='create',BALANCE_NEW_DERIVED_MAPS=True)
    for sample in samples:
        prov=read_json(sample['provenance_path'])
        prov.update(library_ids=[sample['sample_id']+'_library'],source_ids=[sample['sample_id']+'_source'])
        atomic_json(sample['provenance_path'],prov)
    atomic_json(Path(c['run_dir'])/'00_metadata/inputs.json',
                {s['matrix_path']:{'sha256':file_hash(s['matrix_path'])} for s in samples})
    def accepted_owned_weights(c,path,blacklist):
        # Numerical balancing is tested separately. Here the oracle is raw sum
        # preservation on remapped coordinates, including native chromosome reversal.
        assert Path(path).is_relative_to(Path(c['run_dir'])/'02_pools')
        with h5py.File(path,'r+') as h:
            values=np.ones(len(h['bins/start']));values[blacklist]=np.nan
            h['bins'].create_dataset('weight',data=values)
    monkeypatch.setattr(science,'balance_new_map',accepted_owned_weights)
    science.pools(c,'run')
    bins,*_=domain(c,samples,10000)
    expected=cooler.Cooler(samples[1]['matrix_path']).matrix(balance=False)[:,:] * 2
    pooled=cooler.Cooler(str(Path(c['run_dir'])/'02_pools/all_samples/10000.cool'))
    np.testing.assert_array_equal(pooled.matrix(balance=False)[:,:],expected)
    clr,native,mask,labels,views,values=science.context(c,'pool__all_samples',10000)
    assert clr.chromnames==['chr1','chrX','chrY'] and mask.sum()==22


def test_feature_export_uses_registered_order_for_tad_grid(aligned):
    c,_=aligned
    clr,native,mask,labels,_,_=science.context(c,'A',10000)
    measurements=native.assign(score=np.arange(len(native)))
    projected=science.registered_table(c,10000,measurements)
    bins,*_=domain(c,science.validate_metadata(c)['samples'],10000)
    pd.testing.assert_frame_equal(projected[['chrom','start','end']],bins.astype({'chrom':str}),check_dtype=False)
    np.testing.assert_array_equal(projected.score,native_indices(clr,bins))


def test_human_chr22_and_xy_use_registered_coordinates_without_mouse_limit(aligned):
    # Small synthetic chromosome lengths test reference handling, not human biology.
    c,samples=aligned
    registry=read_json(c['REFERENCE_REGISTRY'])
    registry['reference_id']='hg38'
    registry['chromosome_classes']['autosomes']=['chr22']
    for key in ('chromsizes','views'):
        path=Path(registry[key]);path.write_text(path.read_text().replace('chr1\t','chr22\t').replace('\tchr1\n','\tchr22\n'))
        registry[key+'_sha256']=file_hash(path)
    atomic_json(c['REFERENCE_REGISTRY'],registry)
    c['REFERENCE_ID']='hg38'
    for sample in samples:
        raw=cooler.Cooler(sample['matrix_path'])
        native=raw.bins()[:][['chrom','start','end']].astype({'chrom':str})
        native['chrom']=native.chrom.replace({'chr1':'chr22'})
        path=Path(sample['matrix_path']).with_name('hg38-'+sample['sample_id']+'.cool')
        values=raw.bins()['weight'][:].to_numpy()
        cooler.create_cooler(str(path),native,raw.pixels()[:],ordered=True,assembly='hg38')
        with h5py.File(path,'r+') as h:
            ds=h['bins'].create_dataset('weight',data=values);ds.attrs['converged']=True
        sample['matrix_path']=str(path)
        prov=read_json(sample['provenance_path']);prov['reference_id']='hg38'
        atomic_json(sample['provenance_path'],prov)
    bins,common,*_=domain(c,samples,10000)
    assert list(dict.fromkeys(bins.chrom))==['chr22','chrX','chrY']
    assert common.sum()==22
