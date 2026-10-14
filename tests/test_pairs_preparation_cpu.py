"""Exact coordinate preparation, convergence, conservation and fail-closed provenance."""
from pathlib import Path
import json
import sys

import cooler
import jsonschema
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import ContractError, atomic_json, file_hash, write_table, read_json
from microc2compare.pairs_preparation import prepare_from_pairs


def setup_pairs(tmp_path,config,minimum=0,old_minimum=None):
    pairs=tmp_path/"retained.pairs"
    positions=np.arange(0,40000,200)+113
    i,j=np.triu_indices(len(positions))
    kept=positions[j]-positions[i]>=minimum
    i,j=i[kept],j[kept]
    with pairs.open("w",encoding="utf-8") as handle:
        handle.write("## pairs format v1.0\n#genome_assembly: custom\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n")
        for k,(a,b) in enumerate(zip(positions[i],positions[j])):
            handle.write(f"r{k}\tcontigA\t{a}\tcontigA\t{b}\t+\t-\tUU\n")
        handle.write("excluded\tcontigA\t113\tcontigA\t513\t+\t-\tDD\n")
        handle.write("unmapped\t!\t0\tcontigA\t513\t+\t-\tUU\n")
    sizes=tmp_path/"sizes.tsv";sizes.write_text("contigA\t40000\n")
    views=tmp_path/"views.bed";views.write_text("contigA\t0\t40000\tcontigA\n")
    ref=tmp_path/"reference.json";atomic_json(ref,{"reference_id":"custom","chromsizes":str(sizes),"chromsizes_sha256":file_hash(sizes),
        "views":str(views),"views_sha256":file_hash(views),"chromosome_classes":{"autosomes":["contigA"],"sex":[]}})
    provenance=tmp_path/"sample.json"
    atomic_json(provenance,{"reference_id":"custom","pairs_provenance":{
        "reference_id":"custom","source_sha256":file_hash(pairs),"complete":True,"deduplicated":True,
        "dedup_scope":"joint_library","coordinate_system":"1-based","records_total":len(i)+2,
        "upstream_filters":{"retained_pair_types":["UU"],"minimum_cis_distance_bp":minimum,"mapq":30,"chromosome_policy":"contigA"},
        "balancing_recipe":{"cis_only":False,"ignore_diags":0,"mad_max":5,"min_nnz":1,"min_count":0,"tol":1e-6,"max_iters":500}}})
    if old_minimum is not None:
        original=read_json(provenance)
        original.update(count_units="unique_contacts",filters={"mapq":50,"min_distance_bp":old_minimum,"chromosome_policy":"historical_subset"},
                        balancing={"10000":{"recipe":{"historical":True}}},dedup_scope="historical",library_ids=["old_lib"],source_ids=["old_source"])
        atomic_json(provenance,original)
    sheet=tmp_path/"samples.tsv"
    row={"sample_id":"new-1.v2","project":"newproject","condition":"one","assay":"microc","biological_replicate":"1",
         "subject_id":"subject1","reference_id":"custom","matrix_path":str(tmp_path/"not-created.cool"),"provenance_path":str(provenance),"pairs_path":str(pairs)}
    write_table(sheet,[row],list(row))
    membership=tmp_path/"technical.tsv"
    source={"sample_id":"new-1.v2","library_id":"lib1","source_id":"src1","source_kind":"independent_library","dedup_group_id":"d1","dedup_scope":"joint_library","source_checksum":file_hash(pairs)}
    write_table(membership,[source],list(source))
    config.update(PROJECT="newproject",REFERENCE_ID="custom",REFERENCE_REGISTRY=str(ref),SAMPLESHEET=str(sheet),TECHNICAL_MEMBERSHIP=str(membership),
                  RUN_DIFFERENTIAL=False,RUN_TADS=False,RUN_INSULATION=False,THREADS_BALANCE=2,TOTAL_CPUS=2,MEMORY_MB_BALANCE=4000,PIXEL_CHUNKSIZE=1000)
    return pairs,provenance,len(i)


def test_exact_pairs_construct_400_600_800_with_conservation_and_new_weights(tmp_path,config):
    pairs,provenance,retained=setup_pairs(tmp_path,config)
    before={str(p):file_hash(p) for p in (pairs,provenance,Path(config["REFERENCE_REGISTRY"]),Path(config["SAMPLESHEET"]))}
    result=prepare_from_pairs(config,tmp_path/"new",[400,600,800])
    registry=json.loads(Path(result["matrix_registry"]).read_text())
    for grid in (400,600,800):
        clr=cooler.Cooler(registry["samples"]["new-1.v2"][str(grid)])
        assert clr.binsize==grid and clr.info["sum"]==retained
        native=clr.bins()[:]
        assert np.isfinite(native.loc[(native.end-native.start)==grid,"weight"]).all()
        with clr.open("r") as handle:assert bool(handle["bins/weight"].attrs["converged"])
    # Exact one-based position 713 becomes 712: its 600 bp bin is 600,
    # which cannot be obtained by splitting the original 1000 bp bin.
    clr=cooler.Cooler(registry["samples"]["new-1.v2"]["600"])
    assert clr.matrix(balance=False)[1,1]==6
    preparation=json.loads((tmp_path/"new/preparation.json").read_text())
    assert all(r["records"]==retained+2 and r["excluded_pair_type"]==r["excluded_unmapped"]==1 for r in preparation["rows"])
    assert before=={p:file_hash(p) for p in before}
    assert result["scientific_admission"]=="not_performed"
    active=read_json(tmp_path/"new/new-1.v2.provenance.json")
    schema=read_json(Path(__file__).resolve().parents[1]/"schemas/provenance.schema.json")
    jsonschema.validate(active,schema)
    assert active["library_ids"]==["lib1"] and active["source_ids"]==["src1"]
    assert active["count_units"]=="unique_contacts" and active["filters"]["min_distance_bp"]==0
    assert set(active["balancing"])=={"400","600","800"}


@pytest.mark.parametrize("field,value",[("coordinate_system","unknown"),("deduplicated",False),("complete",False),("source_sha256","stale")])
def test_incomplete_pairs_evidence_fails_before_creating_output(tmp_path,config,field,value):
    _,provenance,_=setup_pairs(tmp_path,config)
    data=json.loads(provenance.read_text());data["pairs_provenance"][field]=value;atomic_json(provenance,data)
    with pytest.raises(ContractError):prepare_from_pairs(config,tmp_path/"new",[400])
    assert not (tmp_path/"new").exists()


def test_wrong_coordinate_convention_never_publishes_registry(tmp_path,config):
    pairs,provenance,_=setup_pairs(tmp_path,config)
    content=pairs.read_text().replace("r0\tcontigA\t113\tcontigA\t113", "r0\tcontigA\t0\tcontigA\t113")
    pairs.write_text(content)
    data=json.loads(provenance.read_text());data["pairs_provenance"]["source_sha256"]=file_hash(pairs);atomic_json(provenance,data)
    with pytest.raises(ContractError,match="coordinate_system"):prepare_from_pairs(config,tmp_path/"new",[400])
    assert not (tmp_path/"new/matrix_registry.json").exists()
    assert (tmp_path/"new/preparation.failed.json").exists()


def test_existing_destination_is_never_overwritten(tmp_path,config):
    setup_pairs(tmp_path,config);destination=tmp_path/"new";destination.mkdir();sentinel=destination/"sentinel";sentinel.write_text("keep")
    with pytest.raises(ContractError,match="NEW directory"):prepare_from_pairs(config,destination,[400])
    assert sentinel.read_text()=="keep"


def test_pairs_header_sizes_must_match_reference(tmp_path,config):
    pairs,provenance,_=setup_pairs(tmp_path,config)
    pairs.write_text(pairs.read_text().replace("#genome_assembly: custom", "#genome_assembly: custom\n#chromsize: contigA 50000"))
    data=json.loads(provenance.read_text());data["pairs_provenance"]["source_sha256"]=file_hash(pairs);atomic_json(provenance,data)
    with pytest.raises(ContractError,match="chromosome sizes contradict"):prepare_from_pairs(config,tmp_path/"new",[400])
    assert not (tmp_path/"new").exists()


@pytest.mark.parametrize("field,value",[("mapq",None),("mapq",True),("chromosome_policy",""),("min_distance_bp",1000)])
def test_pairs_active_filter_declarations_fail_closed(tmp_path,config,field,value):
    _,provenance,_=setup_pairs(tmp_path,config)
    original=read_json(provenance)
    original["pairs_provenance"]["upstream_filters"][field]=value
    atomic_json(provenance,original)
    with pytest.raises(ContractError,match="filters|conflicting aliases"):
        prepare_from_pairs(config,tmp_path/"new",[2000])
    assert not (tmp_path/"new").exists()


@pytest.mark.parametrize("old_minimum,minimum",[(8000,1000),(1000,8000)])
def test_pairs_active_cutoff_controls_real_preflight_and_preserves_history(tmp_path,config,old_minimum,minimum):
    from microc2compare.science import inspect_inputs
    pairs,provenance,_=setup_pairs(tmp_path,config,minimum,old_minimum)
    original=read_json(provenance);before={str(p):file_hash(p) for p in (pairs,provenance)}
    result=prepare_from_pairs(config,tmp_path/"new",[2000])
    active=read_json(tmp_path/"new/new-1.v2.provenance.json")
    assert active["source_provenance"]==original
    assert active["filters"]=={"min_distance_bp":minimum,"mapq":30,"chromosome_policy":"contigA","retained_pair_types":["UU"]}
    assert active["dedup_scope"]=="joint_library" and set(active["balancing"])=={"2000"}
    config.update(SAMPLESHEET=result["samplesheet"],MATRIX_REGISTRY=result["matrix_registry"],
                  ANALYSIS_PROFILE="existing",ENVIRONMENT_MODE="development",PYTHON_BOOTSTRAP=sys.executable,
                  OUTDIR=str(tmp_path/"results"),SCRATCH_ROOT=str(tmp_path/"scratch"),MAX_MEMORY_MB=2000,
                  RUN_QC=False,RUN_HICREP=False,RUN_COMPARTMENTS=False,RUN_LOOPS=False,RUN_AGGREGATES=False,
                  RUN_CHROMOSIGHT=False,RUN_CHROMOSIGHT_DISCOVERY=False,RUN_BROWSER_EXPORT=False,
                  LOOP_RESOLUTION=2000,LOOP_QUANTIFICATION_RESOLUTIONS="2000",_requested_analyses=["quantify"],
                  LOOP_MIN_DISTANCE_BP=minimum-1)
    with pytest.raises(ContractError,match="below upstream contact cutoff"):
        inspect_inputs(config)
    config["LOOP_MIN_DISTANCE_BP"]=minimum
    assert inspect_inputs(config)["status"]=="validated"
    assert before=={path:file_hash(path) for path in before}


def test_invalid_prepared_entry_never_publishes_registry(tmp_path,config,monkeypatch):
    import h5py
    from microc2compare import pairs_preparation
    setup_pairs(tmp_path,config)
    validate=pairs_preparation._validate_prepared_outputs
    def corrupt_then_validate(c,updated,registry,sizes,membership):
        target=registry["samples"]["new-1.v2"]["2000"]
        with h5py.File(target,"r+") as handle:
            handle["bins/weight"].attrs["converged"]=False
        return validate(c,updated,registry,sizes,membership)
    monkeypatch.setattr(pairs_preparation,"_validate_prepared_outputs",corrupt_then_validate)
    with pytest.raises(ContractError,match="lack convergence"):
        prepare_from_pairs(config,tmp_path/"new",[2000])
    assert not (tmp_path/"new/matrix_registry.json").exists()
    assert (tmp_path/"new/preparation.failed.json").exists()
