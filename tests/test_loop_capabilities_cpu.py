"""Real caller/kernel and Cooler fixtures for generic CPU loop preparation."""
from importlib.resources import files
import json
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import ContractError, atomic_json, file_hash
from microc2compare.fine_loops import chromosight_discovery, chromosight_scores, cross_resolution_catalogue, _mustache_sparse_band, _expected_for
from microc2compare.loop_profiles import HSIEH_GRIDS, profile_defaults, mustache_geometry, mustache_resource_plan
from microc2compare.resolutions import input_capabilities, chromosight_discovery_resolutions
from microc2compare import science


def make_map(path, resolution=400, n=180, assembly="custom", peak=(55,65)):
    bins = cooler.binnify(pd.Series({"contigA": n*resolution}), resolution)
    counts = np.full((n,n), 10000, dtype=np.int64)
    if peak:
        pytest.importorskip("chromosight.utils.detection")
        kernel = np.loadtxt(str(files("chromosight.kernels").joinpath("artificial_template_loops_type1.txt")))
        x, y = peak; r=kernel.shape[0]//2
        counts[x-r:x+r+1,y-r:y+r+1] += np.rint((kernel-1)*10000).astype(np.int64)
    i,j=np.triu_indices(n)
    cooler.create_cooler(str(path), bins, pd.DataFrame({"bin1_id":i,"bin2_id":j,"count":counts[i,j]}), ordered=True, assembly=assembly)
    with h5py.File(path,"r+") as handle:
        weight=handle["bins"].create_dataset("weight", data=np.linspace(.9,1.1,n))
        weight.attrs["converged"]=True
    return cooler.Cooler(str(path)), bins


def attach_context(monkeypatch, clr, bins, config):
    n=len(bins); mask=np.ones(n,dtype=bool); labels=np.zeros(n,dtype=np.int64)
    views=pd.DataFrame([["contigA",0,int(bins.end.iloc[-1]),"viewA"]],columns=["chrom","start","end","name"])
    values=clr.bins()[:].weight.to_numpy()
    monkeypatch.setattr(science,"context",lambda *args:(clr,bins,mask,labels,views,values))
    return bins,mask,labels,views,values


def test_profiles_have_no_cohort_or_reference_and_divisible_apa():
    full=profile_defaults("hsieh_8grid")
    assert tuple(map(int,full["LOOP_DISCOVERY_RESOLUTIONS"].split(",")))==HSIEH_GRIDS
    assert "REFERENCE_ID" not in full and "SAMPLESHEET" not in full
    assert all(full["APA_FLANK_BP"]%r==0 for r in HSIEH_GRIDS)
    short=profile_defaults("short_range")
    assert short["ANALYSIS_PROFILE"]=="short_range_4kb" and short["LOOP_MAX_DISTANCE_BP"]==500000
    full["LOOP_RESOLUTION"]=999
    assert profile_defaults("hsieh_8grid")["LOOP_RESOLUTION"]==400


def test_two_callers_preserve_identity_without_inflating_biological_support(config):
    first=pd.DataFrame([{"chrom1":"chr1","start1":4000,"end1":4400,"chrom2":"chr1","start2":8000,"end2":8400,"caller":"mustache"}])
    second=first.assign(caller="chromosight")
    candidates,members=cross_resolution_catalogue(config,[("sample-1",400,first),("sample-1",400,second)])
    assert len(candidates)==1 and members.source_call_id.nunique()==2
    assert candidates.biological_support.iloc[0]==1 and candidates.source_count.iloc[0]==2
    assert candidates.supporting_callers.iloc[0]=="chromosight;mustache"


def test_geometry_uses_real_first_assigned_scale_and_diagonal_occupancy(config):
    config.update(MUSTACHE_SIGMA_ZERO=1.6,MUSTACHE_OCTAVES=2,MUSTACHE_ITERATIONS=10,MUSTACHE_SPARSITY_THRESHOLD=.88,LOOP_MIN_DISTANCE_BP=1000)
    # Independent construction of ideal original nz mask in the caller's
    # actual diagonal convention verifies the 19/25 versus 22/25 distinction.
    nz=np.triu(np.ones((80,80),dtype=bool),4)
    assert nz[28:33,33:38].sum()==19
    assert nz[28:33,34:39].sum()==22
    for resolution in HSIEH_GRIDS:
        result=mustache_geometry(config,resolution)
        assert result["inner_radius_bins"]==2
        assert result["geometric_min_distance_bins"]==6
        assert result["geometric_min_distance_bp"]==6*resolution
        assert result["first_assigned_scale_bins"]==pytest.approx(1.6*2**.2)
    assert mustache_geometry(config,800)["geometric_min_distance_bp"]==4800
    assert mustache_geometry(config,400)["occupancy_at_five_bins"]["inner"]==.76


def test_adaptive_workers_admit_fine_grids_and_account_for_sparse_arrays(config):
    config.update(THREADS_LOOPS=16,TOTAL_CPUS=72,MEMORY_MB_LOOPS=80000,LOOP_MAX_DISTANCE_BP=2000000,MAX_DENSE_BINS=12000)
    assert [mustache_resource_plan(config,r)["workers"] for r in (400,600,800,1000)]==[4,9,16,16]
    assert mustache_resource_plan(config,400,1000)["workers"]==3
    config["MEMORY_MB_LOOPS"]=10000
    with pytest.raises(ContractError,match="requires at least"):
        mustache_resource_plan(config,400)


@pytest.mark.parametrize("assembly,sample_count",[("hg38",1),("mm39",2),("custom",5)])
def test_native_capabilities_are_generic_and_never_upsample(tmp_path,config,assembly,sample_count):
    config.update(REFERENCE_ID=assembly,MATRIX_REGISTRY="")
    samples=[]
    for i in range(sample_count):
        path=tmp_path/f"sample{i}.cool";make_map(path,1000,40,assembly,peak=None)
        samples.append({"sample_id":f"sample-{i}.v1","matrix_path":str(path),"reference_id":assembly})
    before={s["sample_id"]:file_hash(s["matrix_path"]) for s in samples}
    report=input_capabilities(config,samples,[400,600,800,1000,2000,4000])
    assert report["sample_count"]==sample_count
    for sid in before:
        rows={r["resolution"]:r for r in report["rows"] if r["sample_id"]==sid}
        assert [rows[r]["status"] for r in (400,600,800)]==["missing_native_grid"]*3
        assert rows[1000]["status"]=="available"
        assert rows[2000]["status"]==rows[4000]["status"]=="coarsenable"
        assert "exact-position pairs" in rows[400]["action"]
    assert before=={s["sample_id"]:file_hash(s["matrix_path"]) for s in samples}


def test_registry_and_assembly_mismatch_are_not_silently_admitted(tmp_path,config):
    good=tmp_path/"good.cool";bad=tmp_path/"bad.cool"
    make_map(good,1000,40,"custom",peak=None);make_map(bad,400,40,"other",peak=None)
    registry=tmp_path/"registry.json"
    atomic_json(registry,{"reference_id":"custom","samples":{"one":{"400":str(bad)}}})
    config.update(REFERENCE_ID="custom",MATRIX_REGISTRY=str(registry))
    report=input_capabilities(config,[{"sample_id":"one","matrix_path":str(good)}],[400])
    assert report["rows"][0]["status"]=="invalid_native_metadata"
    assert not any(r.get("assembly_matches") for r in report["native_inputs"]["one"] if r.get("resolution")==400)


def test_banded_mustache_reads_equal_full_fetch_without_source_mutation(tmp_path,config):
    path=tmp_path/"map.cool";clr,bins=make_map(path,400,650,peak=None)
    mask=np.ones(len(bins),dtype=bool);mask[520]=False
    values=clr.bins()[:].weight.to_numpy();before=file_hash(path)
    config.update(LOOP_MAX_DISTANCE_BP=20000,MAX_DENSE_BINS=12000,MEMORY_MB_LOOPS=4000,THREADS_LOOPS=2,TOTAL_CPUS=2)
    x,y,v=_mustache_sparse_band(config,clr,0,len(bins),mask,values,50,400)
    raw=clr.matrix(balance=False,sparse=True).fetch("contigA").tocoo()
    keep=(raw.col>=raw.row)&(raw.col-raw.row<=50)&mask[raw.row]&mask[raw.col]
    expected=pd.DataFrame({"x":raw.row[keep],"y":raw.col[keep],"v":raw.data[keep]*values[raw.row[keep]]*values[raw.col[keep]]}).sort_values(["x","y"]).reset_index(drop=True)
    actual=pd.DataFrame({"x":x,"y":y,"v":v}).sort_values(["x","y"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(actual,expected,check_dtype=False)
    assert file_hash(path)==before


def test_real_independent_chromosight_finds_four_kb_and_preserves_input(tmp_path,config,monkeypatch):
    path=tmp_path/"input.cool";clr,bins=make_map(path)
    attach_context(monkeypatch,clr,bins,config)
    config.update(run_dir=str(tmp_path/"run"),REFERENCE_ID="custom",LOOP_RESOLUTION=400,LOOP_DISCOVERY_RESOLUTIONS="400",CHROMOSIGHT_DISCOVERY_RESOLUTIONS="",
                  LOOP_MIN_DISTANCE_BP=1000,LOOP_MAX_DISTANCE_BP=20000,CHROMOSIGHT_DISCOVERY_MIN_DISTANCE_BP=0,CHROMOSIGHT_DISCOVERY_MAX_DISTANCE_BP=0,
                  CHROMOSIGHT_DISCOVERY_MIN_PEARSON=.3,CHROMOSIGHT_MAX_MISSING_FRACTION=.1,THREADS_CHROMOSIGHT=2,MEMORY_MB_CHROMOSIGHT=2000,TOTAL_CPUS=2,MAX_DENSE_BINS=512,INSULATION_IGNORE_DIAGS=2)
    before=file_hash(path);products=chromosight_discovery(config,"sample-1.v1")
    frame=pd.read_csv(products[0],sep="\t")
    assert len(frame)>0
    closest=frame.iloc[np.argmin(abs(frame.start1-22000)+abs(frame.start2-26000))]
    # O/E and partial near-diagonal kernel support can shift the local maximum
    # by one native bin; retain measured coordinates rather than snapping it.
    assert abs(closest.start1-22000)<=400 and abs(closest.start2-26000)<=400
    assert abs(closest.distance_bp-4000)<=400 and closest.caller_score>.8
    assert frame.caller_qvalue.isna().all()
    assert file_hash(path)==before
    assert not list(Path(config["run_dir"]).rglob("*.cool"))
    provenance=json.loads(products[1].read_text())
    assert provenance["independent_discovery"] and not provenance["input_mutation"]
    assert chromosight_discovery_resolutions(config)==[400]


def test_native_discovery_covers_tile_seam_and_is_thread_deterministic(tmp_path,config,monkeypatch):
    path=tmp_path/"input.cool";clr,bins=make_map(path,400,700,peak=(443,488))
    attach_context(monkeypatch,clr,bins,config)
    config.update(run_dir=str(tmp_path/"run"),REFERENCE_ID="custom",LOOP_RESOLUTION=400,LOOP_DISCOVERY_RESOLUTIONS="400",CHROMOSIGHT_DISCOVERY_RESOLUTIONS="400",
                  LOOP_MIN_DISTANCE_BP=1000,LOOP_MAX_DISTANCE_BP=20000,CHROMOSIGHT_DISCOVERY_MIN_DISTANCE_BP=0,CHROMOSIGHT_DISCOVERY_MAX_DISTANCE_BP=0,
                  CHROMOSIGHT_DISCOVERY_MIN_PEARSON=.3,CHROMOSIGHT_MAX_MISSING_FRACTION=.1,THREADS_CHROMOSIGHT=2,MEMORY_MB_CHROMOSIGHT=2000,TOTAL_CPUS=2,MAX_DENSE_BINS=512,INSULATION_IGNORE_DIAGS=2)
    products=chromosight_discovery(config,"seam");parallel=pd.read_csv(products[0],sep="\t")
    nearest=parallel.iloc[np.argmin(abs(parallel.start1-443*400)+abs(parallel.start2-488*400))]
    assert abs(nearest.start1-443*400)<=400 and abs(nearest.start2-488*400)<=400 and nearest.caller_score>.8
    assert parallel.feature_id.is_unique
    config["THREADS_CHROMOSIGHT"]=1
    serial=pd.read_csv(chromosight_discovery(config,"seam")[0],sep="\t")
    pd.testing.assert_frame_equal(parallel,serial)


def test_expected_reuse_requires_valid_qc_and_listed_output(tmp_path,config,monkeypatch):
    from microc2compare import state
    config["run_dir"]=str(tmp_path)
    expected=pd.DataFrame([{"region1":"view","region2":"view","dist":0,"balanced.avg":2.}])
    source=tmp_path/"01_qc/ps/sample.400.tsv";science.save(expected,source)
    calls=[]
    monkeypatch.setattr(science,"expected_table",lambda *args:calls.append(1) or expected)
    monkeypatch.setattr(state,"complete",lambda *args:False)
    _,origin=_expected_for(config,"sample",400,None,None,None,None,None,None)
    assert origin=="computed_current_context" and len(calls)==1
    record=state.record_path(config,"qc","sample")
    atomic_json(record,{"outputs":[{"path":str(source.relative_to(tmp_path))}]})
    monkeypatch.setattr(state,"complete",lambda *args:True)
    lookup,origin=_expected_for(config,"sample",400,None,None,None,None,None,None)
    assert origin=="completed_qc" and lookup["view"][0]==2 and len(calls)==1
    atomic_json(record,{"outputs":[]})
    _,origin=_expected_for(config,"sample",400,None,None,None,None,None,None)
    assert origin=="computed_current_context" and len(calls)==2


def test_fixed_scoring_handles_safe_punctuation_ids_in_parallel(tmp_path,config,monkeypatch):
    path=tmp_path/"input.cool";clr,bins=make_map(path,400,180,peak=(55,100))
    attach_context(monkeypatch,clr,bins,config)
    scope="sample-1.v1"
    config.update(run_dir=str(tmp_path/"run"),REFERENCE_ID="custom",LOOP_RESOLUTION=400,LOOP_QUANTIFICATION_RESOLUTIONS="400",RUN_LOOPS=True,
                  THREADS_CHROMOSIGHT=2,MEMORY_MB_CHROMOSIGHT=2000,TOTAL_CPUS=2,CHROMOSIGHT_MAX_MISSING_FRACTION=.1)
    candidates=pd.DataFrame([{ "chrom1":"contigA","start1":22000,"end1":22400,"chrom2":"contigA","start2":40000,"end2":40400,"feature_id":"yes",scope:100},
                             { "chrom1":"contigA","start1":22400,"end1":22800,"chrom2":"contigA","start2":40400,"end2":40800,"feature_id":"yes-two",scope:100},
                             { "chrom1":"contigA","start1":24000,"end1":24400,"chrom2":"contigA","start2":42000,"end2":42400,"feature_id":"excluded",scope:np.nan}])
    science.save(candidates,Path(config["run_dir"])/"04_quantification/discovered_loops/400/counts.tsv")
    before=file_hash(path);products=chromosight_scores(config,scope)
    frame=pd.read_csv(products[0],sep="\t")
    assert frame.status.tolist()==["measured","measured","ineligible_candidate"] and frame.score.iloc[0]>.8
    assert json.loads(products[1].read_text())["workers"]==2
    assert file_hash(path)==before
    config["THREADS_CHROMOSIGHT"]=1
    serial=pd.read_csv(chromosight_scores(config,scope)[0],sep="\t")
    pd.testing.assert_frame_equal(frame,serial)
