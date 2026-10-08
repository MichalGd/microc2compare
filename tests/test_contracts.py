import csv
import json
import pytest
from microc2compare.contracts import ContractError,concurrency,formula_terms,normalize,read_table
from microc2compare.science import design_matrix


def test_defaults_reuse(config,tmp_path):
    c=normalize({k:str(v).lower() if isinstance(v,bool) else str(v) for k,v in config.items()},str(tmp_path))
    assert c["RUN_POOL"] is False and c["RUN_DEPTH_SENSITIVITY"] is False
    assert not c["BALANCE_NEW_DERIVED_MAPS"] and not c["KEEP_POOLED_MAPS"]
    assert c["SADDLE_PARTITION_SOURCE"]=="replicate_rank_mean"


def test_csv_quotes_crlf(tmp_path):
    path=tmp_path/"sample table.csv"
    path.write_bytes(b'id,path\r\nS1,"a,b with spaces"\r\n')
    assert read_table(path,("id","path"))[0]["path"]=="a,b with spaces"


@pytest.mark.parametrize("formula",["~condition+system('x')","~batch*condition","~condition+condition","condition"])
def test_reject_executable_or_unsupported_formula(formula):
    with pytest.raises(ContractError):formula_terms(formula)


def test_resource_admission(config):
    assert concurrency(config,8,64000,matrix=True)==2
    with pytest.raises(ContractError):concurrency(config,32,1000)
    with pytest.raises(ContractError):concurrency(config,1,200000)


def test_paired_df_and_confounding(config):
    sample=[{"sample_id":f"{s}_{c}","subject_id":s,"condition":c,"batch":c} for s in ("a","b") for c in ("control","treated")]
    cmp={"comparison_id":"test","reference_condition":"control","case_condition":"treated","design":"~subject_id+condition","inference":"preliminary"}
    with pytest.raises(ContractError,match="degrees"):design_matrix(sample,cmp,config)
    cmp["design"]="~batch+condition"
    with pytest.raises(ContractError,match="rank"):design_matrix(sample,cmp,config)


@pytest.mark.parametrize("inference", ["descriptive", "preliminary", "formal"])
def test_case_minus_reference(config, inference):
    sample=[{"sample_id":str(i),"condition":"control" if i<3 else "treated"} for i in range(6)]
    cmp={"comparison_id":"test","reference_condition":"control","case_condition":"treated","design":"~condition","inference":inference}
    d=design_matrix(sample,cmp,config)
    assert d["contrast"]==[0,1] and d["residual_df"]==4
    assert json.loads(json.dumps(d, allow_nan=False)) == d
