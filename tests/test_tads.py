import numpy as np
import pandas as pd
import pytest

from microc2compare.tads import (DOMAIN_COLUMNS, call_domains, condition_domains,
                                 match_domains, run)


def domains(intervals, strength=1.0):
    return pd.DataFrame([dict(feature_id=f"d_{a}_{b}", chrom="chr1", start=a, end=b, region="arm",
                              window_bp=100000, length_bp=b-a, left_boundary_strength=strength,
                              right_boundary_strength=strength, mean_log2_insulation=0., support=3,
                              n_samples=3, support_fraction=1., sample_ids="1;2;3", inference_status="descriptive")
                         for a, b in intervals], columns=DOMAIN_COLUMNS)


COMPARISON = dict(comparison_id="case_vs_reference", reference_condition="control", case_condition="treated")


@pytest.mark.parametrize("left,right,event", [
    ([(100000, 300000)], [(100000, 300000)], "stable"),
    ([(100000, 300000)], [(100000, 400000)], "growing"),
    ([(100000, 400000)], [(100000, 300000)], "shrinking"),
    ([(100000, 200000), (200000, 300000)], [(100000, 300000)], "merge"),
    ([(100000, 300000)], [(100000, 200000), (200000, 300000)], "split"),
    ([], [(100000, 300000)], "appearing"),
    ([(100000, 300000)], [], "disappearing"),
    ([(100000, 300000)], [(130000, 330000)], "shifted"),
    ([(100000, 220000), (220000, 340000)], [(100000, 280000), (280000, 340000)], "complex"),
])
def test_known_events(left, right, event, config):
    frame = match_domains(domains(left), domains(right), config, COMPARISON, 100000)
    assert len(frame) == 1 and frame.event.iloc[0] == event
    assert frame.PValue.isna().all() and frame.FDR.isna().all()
    assert frame.event_id.is_unique
    if event == "merge": assert frame.reference_ids.iloc[0].count(";") == 1


def test_strength_missing_coverage_and_disagreement(config):
    left = domains([(100000, 300000)])
    right = domains([(100000, 300000)], strength=2.)
    result = match_domains(left, right, config, COMPARISON, 100000)
    assert result.event.iloc[0] == "strength_change"
    assert result.boundary_strength_change.iloc[0] == 1
    assert match_domains(left, domains([]), config, COMPARISON, 100000,
                         evaluable=lambda *a: 0.).event.iloc[0] == "uncertain"
    assert match_domains(left, domains([]), config, COMPARISON, 100000,
                         ambiguous=lambda *a: True).reason.iloc[0] == "replicate_domain_disagreement"


def test_direction_reversal_and_resolution_tolerance(config):
    a, b = domains([(100000, 300000)]), domains([(100000, 400000)])
    forward = match_domains(a, b, config, COMPARISON, 100000).iloc[0]
    reverse = match_domains(b, a, config, COMPARISON, 100000).iloc[0]
    assert forward.length_change_bp == -reverse.length_change_bp
    assert forward.length_ratio * reverse.length_ratio == pytest.approx(1.)
    close = domains([(110000, 310000)])
    assert match_domains(a, close, config, COMPARISON, 100000).event.iloc[0] == "stable"


def score_table():
    frame = pd.DataFrame(dict(chrom=["chr1"]*100, start=np.arange(100)*10000,
                              end=(np.arange(100)+1)*10000, region=["arm"]*100))
    frame["log2_insulation_score_100000"] = 0.
    frame["valid_pixel_fraction_100000"] = 1.
    frame["boundary_strength_100000"] = np.nan
    frame.loc[[10, 40, 70], "boundary_strength_100000"] = [1., .3, 1.]
    return frame


def test_threshold_and_mask_never_bridge_exclusions(config):
    table = score_table(); mask = np.ones(len(table), dtype=bool)
    two = call_domains(table, mask, config, 100000, "s", threshold=0.)
    one = call_domains(table, mask, config, 100000, "s", threshold=.5)
    assert len(two) == 2 and len(one) == 1
    assert one.iloc[0].length_bp == 600000
    mask[30] = False
    assert len(call_domains(table, mask, config, 100000, "s", threshold=.5)) == 0
    table.loc[10, "boundary_strength_100000"] = np.nan
    assert len(call_domains(table, np.ones(len(table), bool), config, 100000, "s", threshold=.5)) == 0


def test_replicates_and_complete_span_clustering(config):
    frames = {"a": domains([(100000, 300000)]), "b": domains([(110000, 310000)]),
              "c": domains([(120000, 320000)])}
    accepted, rejected = condition_domains(frames, config, 100000)
    assert accepted.empty and len(rejected) == 2  # no transitive 20-kb cluster
    frames["c"] = domains([(110000, 310000)])
    accepted, rejected = condition_domains(frames, config, 100000)
    assert len(accepted) == 1 and not len(rejected)
    assert accepted.support.iloc[0] == 3


def test_tad_adapter_exports_events_sensitivity_and_browser_files(tmp_path, config, monkeypatch):
    import microc2compare.tads as module
    c = dict(config, run_dir=str(tmp_path), REFERENCE_ID="ref", INSULATION_WINDOWS_BP="100000", INSULATION_BOUNDARY_THRESHOLD="0.1")
    samples = [dict(sample_id=f"{cond}{i}", condition=cond) for cond in ("control", "treated") for i in range(3)]
    monkeypatch.setattr(module, "validate_metadata", lambda c: dict(samples=samples, comparisons=[COMPARISON]))
    d = tmp_path/"03_features/insulation"; d.mkdir(parents=True)
    table = score_table()
    for sample in samples:
        frame = table.copy()
        if sample["condition"] == "treated": frame.loc[40, "boundary_strength_100000"] = 0.
        frame.to_csv(d/f"{sample['sample_id']}.tsv", sep="\t", index=False)
    meta = tmp_path/"00_metadata/masks"; meta.mkdir(parents=True)
    table[["chrom", "start", "end"]].to_csv(meta/"10000.bins.tsv", sep="\t", index=False)
    np.savez_compressed(meta/"10000.npz", mask=np.ones(len(table), bool))
    products = run(c, "run")
    assert all(p.is_file() and p.stat().st_size > 0 for p in products)
    events = pd.read_csv(tmp_path/"06_comparisons/case_vs_reference/tad_events.100000.tsv", sep="\t")
    assert events.event.tolist() == ["merge"]
    assert len(list((tmp_path/"03_features/tads/samples").glob("*.bed"))) == 6
    sensitivity = pd.read_csv(tmp_path/"03_features/tads/sensitivity.tsv", sep="\t")
    assert set(sensitivity.boundary_threshold) == {0., .1, .2, .5}
    assert not list(tmp_path.rglob("*.cool"))
    # Test the report's TAD consumer with real bounded Cooler reads and known
    # score tables; cooltools insulation computation remains a separate gate.
    import cooler
    monkeypatch.setenv('MPLCONFIGDIR', str(tmp_path/'plot-cache'))
    from microc2compare import report
    from microc2compare.contracts import atomic_json, file_hash
    matrix = tmp_path/'original-report-fixture.cool'
    cooler.create_cooler(str(matrix), table[['chrom','start','end']],
                        pd.DataFrame(dict(bin1_id=np.arange(100), bin2_id=np.arange(100), count=np.ones(100,dtype=int))), ordered=True)
    before = file_hash(matrix)
    for sample in samples: sample.update(assay='microc', subject_id=sample['sample_id'])
    c.update(RUN_QC=False, RUN_HICREP=False, RUN_COMPARTMENTS=False, RUN_LOOPS=False,
             RUN_AGGREGATES=False, RUN_DIFFERENTIAL=False, RUN_POOL=False, RUN_DEPTH_SENSITIVITY=False,
             RUN_INSULATION=True, RUN_TADS=True)
    monkeypatch.setattr(report, 'validate_metadata', lambda c: dict(samples=samples))
    monkeypatch.setattr(report, 'complete', lambda c, stage, scope: stage=='tads')
    monkeypatch.setattr(report, 'scopes', lambda c, stage: [s['sample_id'] for s in samples] if stage=='insulation' else ['run'])
    clr = cooler.Cooler(str(matrix)); views = pd.DataFrame([dict(chrom='chr1',start=0,end=1000000,name='arm')])
    monkeypatch.setattr('microc2compare.science.context', lambda *a: (clr,table,np.ones(100,bool),np.zeros(100),views,np.ones(100)))
    atomic_json(tmp_path/'completion/run/tads.done.json', {'stage':'tads','scope':'run','outputs':[{'path':str(p.relative_to(tmp_path))} for p in products]})
    atomic_json(tmp_path/'00_metadata/inputs.json', {str(matrix):{'path':str(matrix),'sha256':before}})
    stale = tmp_path/'06_comparisons/unused/tad_events.100000.tsv'
    stale.parent.mkdir(parents=True); events.assign(comparison_id='STALE_COMPARISON_SHOULD_NOT_APPEAR').to_csv(stale,sep='\t',index=False)
    report_outputs = report.render(c, 'run')
    page = (tmp_path/'08_reports/index.html').read_text(encoding='utf-8')
    assert 'TAD changes' in page and 'merge' in page and 'No event-level biological P values' in page
    assert 'STALE_COMPARISON_SHOULD_NOT_APPEAR' not in page and 'unused/tad_events' not in page
    assert any(p.name=='tad_events.png' for p in report_outputs)
    assert len(list((tmp_path/'08_reports/figures').glob('locus.*.png')))==6
    assert file_hash(matrix)==before
