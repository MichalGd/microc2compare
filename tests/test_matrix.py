import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
from microc2compare.contracts import ContractError,file_hash
from microc2compare.features import LOOP_COLUMNS,canonical_loops,consensus_loops,quantify
from microc2compare.matrix import MaskedCooler,expected_table,matrix_uri,raw_summary,weights


@pytest.fixture
def tiny(tmp_path):
    bins=pd.DataFrame({"chrom":["chr1"]*4,"start":[0,10000,20000,30000],"end":[10000,20000,30000,40000]})
    pixels=pd.DataFrame({"bin1_id":[0,0,1,1,2],"bin2_id":[1,2,2,3,3],"count":[2,3,5,7,11]})
    path=tmp_path/"upstream.cool";cooler.create_cooler(str(path),bins,pixels,ordered=True)
    with h5py.File(path,"r+") as h:h["bins"].create_dataset("weight",data=[1,2,3,4])
    return path,bins


def test_stored_weights_and_input_immutability(tiny):
    path,bins=tiny;before=file_hash(path)
    clr=MaskedCooler(str(path),np.array([1,2,np.nan,4]))
    assert clr.matrix(balance=False)[0:1,1:2][0,0]==2
    assert clr.matrix(balance=True)[0:1,1:2][0,0]==4
    assert np.isnan(clr.matrix(balance=True)[0:1,2:3][0,0])
    assert clr.bins()["weight"][:].isna().iloc[2]
    assert file_hash(path)==before


@pytest.mark.parametrize("replace", [False, True])
@pytest.mark.parametrize("pairs", [[(1, 3)], [(3, 3), (1, 2), (1, 3), (1, 2)], [], [(0, 0)]])
def test_cooler_annotation_with_virtual_masked_bins(tiny, replace, pairs):
    path, _ = tiny
    before = file_hash(path)
    clr = MaskedCooler(str(path), np.array([1, 2, np.nan, 4]))
    pixels = pd.DataFrame(pairs, columns=["bin1_id", "bin2_id"], dtype="int64")
    pixels["count"] = np.arange(len(pixels)) + 1
    pixels.index = pd.Index(np.arange(len(pixels)) + 10, name="pixel_id")
    # A fully materialized tiny table is an oracle only in this test.
    expected = cooler.annotate(pixels, clr.bins()[:], replace=replace)
    actual = cooler.annotate(pixels, clr.bins(), replace=replace)
    pd.testing.assert_frame_equal(actual, expected)
    assert file_hash(path) == before


def test_masked_bin_selector_column_selection_and_fetch(tiny):
    path, _ = tiny
    before = file_hash(path)
    clr = MaskedCooler(str(path), np.array([1, 2, np.nan, 4]))
    selected = clr.bins()[["start", "weight"]]
    assert list(selected.columns) == ["start", "weight"]
    assert selected.shape == (4,)
    # Native Cooler fetch/slice can use different integer index widths on Windows.
    pd.testing.assert_frame_equal(selected.fetch("chr1:10000-40000"), selected[1:4], check_index_type=False)
    pd.testing.assert_series_equal(clr.bins()["weight"].fetch("chr1:10000-40000"), clr.bins()["weight"][1:4], check_index_type=False)
    assert list(selected[1:4].index) == [1, 2, 3]
    assert np.isnan(selected[2].weight.iloc[0])
    assert file_hash(path) == before


def test_cooler_annotation_reads_only_requested_bin_ranges(tiny, monkeypatch):
    path, _ = tiny
    clr = MaskedCooler(str(path), np.array([1, 2, np.nan, 4]))
    selector = clr.bins()[["chrom", "start", "end"]]
    original = selector._slice
    observed = []

    def bounded(fields, lo, hi):
        observed.append((lo, hi))
        assert hi - lo == 1
        return original(fields, lo, hi)

    monkeypatch.setattr(selector, "_slice", bounded)
    pixels = pd.DataFrame({"bin1_id": [1], "bin2_id": [3]})
    result = cooler.annotate(pixels, selector, replace=True)
    assert observed == [(1, 2), (3, 4)]
    assert result.start1.iloc[0] == 10000 and result.start2.iloc[0] == 30000


def test_expected_denominators_include_zero(tiny,config):
    path,bins=tiny;clr=cooler.Cooler(str(path));mask=np.ones(4,dtype=bool);labels=np.zeros(4,dtype=int)
    config["INSULATION_IGNORE_DIAGS"]=0
    views=pd.DataFrame({"chrom":["chr1"],"start":[0],"end":[40000],"name":["chr1"]})
    ex=expected_table(clr,bins,mask,labels,views,np.ones(4),config).set_index("dist")
    assert ex.loc[1,"n_valid"]==3 and ex.loc[1,"count.sum"]==18
    assert ex.loc[2,"n_valid"]==2 and ex.loc[2,"count.avg"]==5
    assert ex.loc[3,"n_valid"]==1 and ex.loc[3,"count.avg"]==0


def test_exact_raw_counts_exposure_zero_vs_na(tiny,config):
    path,bins=tiny;config["LOOP_MIN_DISTANCE_BP"]=10000;config["LOOP_MAX_DISTANCE_BP"]=30000
    features=pd.DataFrame([["chr1",0,10000,"chr1",30000,40000,"zero"],["chr1",10000,20000,"chr1",30000,40000,"seven"]],columns=LOOP_COLUMNS)
    values,L,reasons=quantify(cooler.Cooler(str(path)),features,bins,np.ones(4,dtype=bool),np.zeros(4,dtype=int),config)
    assert list(values)==[0,7] and L==28 and not reasons
    mask=np.array([True,False,True,True]);values,L,reasons=quantify(cooler.Cooler(str(path)),features,bins,mask,np.zeros(4,dtype=int),config)
    assert values.iloc[0]==0 and pd.isna(values.iloc[1]) and L==14


def test_missing_resolution_is_not_rebuilt(tiny):
    with pytest.raises(ValueError,match="upstream"):matrix_uri(str(tiny[0]),100000)


def test_nonconverged_weights_are_not_repaired(tiny,config):
    prov={"balancing":{"10000":{"convention":"multiplicative","converged":False}}}
    with pytest.raises(ContractError,match="convergence"):weights(cooler.Cooler(str(tiny[0])),config,prov,10000)


def test_consensus_does_not_chain():
    frame=pd.DataFrame([["chr1",a,a+10000,"chr1",a+100000,a+110000] for a in (0,10000,20000)],columns=LOOP_COLUMNS[:6])
    result=consensus_loops({"S1":frame},"ref",10000,1)
    assert len(result)==2


def test_symmetric_square_pixels_count_once(tiny,tmp_path,config):
    path,bins=tiny;upper=cooler.Cooler(str(path)).pixels()[:]
    lower=upper.rename(columns={'bin1_id':'bin2_id','bin2_id':'bin1_id'})[upper.columns]
    pixels=pd.concat([upper,lower]).sort_values(['bin1_id','bin2_id'])
    square=tmp_path/'square.cool';cooler.create_cooler(str(square),bins,pixels,ordered=True,symmetric_upper=False)
    clr=cooler.Cooler(str(square));config['LOOP_MIN_DISTANCE_BP']=10000;config['LOOP_MAX_DISTANCE_BP']=30000
    summary=raw_summary(clr,bins,np.ones(4,dtype=bool),np.zeros(4,dtype=int),config)
    assert clr.info['sum']==56 and summary['raw_total']==summary['eligible_cis_exposure']==28
