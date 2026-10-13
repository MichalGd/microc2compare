"""Numerical regression against the former per-distance dictionary algorithm."""
import cooler
import numpy as np
import pandas as pd
import pytest
from scipy.signal import fftconvolve

from microc2compare.contracts import file_hash
from microc2compare.matrix import expected_table, pixel_chunks

EXPECTED_COLUMNS = ["region1", "region2", "dist", "distance_bp", "n_valid", "count.sum", "balanced.sum", "count.avg", "balanced.avg"]


def legacy_expected_table(clr, bins, mask, labels, views, values, c):
    # Freeze the previous algorithm here as the numerical regression oracle.
    rows, accumulators = [], {}
    for vid, view in enumerate(views.itertuples()):
        ids = np.flatnonzero(labels == vid)
        if not len(ids):
            continue
        valid = mask[ids].astype(float)
        denominator = np.rint(fftconvolve(valid, valid[::-1], mode="full")[len(valid)-1:]).astype(np.int64)
        accumulators[vid] = {"denom": denominator, "raw": np.zeros(len(ids)), "balanced": np.zeros(len(ids))}
    for pixels in pixel_chunks(clr, c["PIXEL_CHUNKSIZE"]):
        i,j,count = [pixels[k].to_numpy() for k in ("bin1_id", "bin2_id", "count")]
        valid = mask[i] & mask[j] & (labels[i] >= 0) & (labels[i] == labels[j])
        for vid in np.unique(labels[i[valid]]):
            select = valid & (labels[i] == vid)
            d = j[select]-i[select]
            a = accumulators[vid]
            a["raw"] += np.bincount(d, weights=count[select], minlength=len(a["denom"]))
            a["balanced"] += np.bincount(d, weights=count[select]*values[i[select]]*values[j[select]], minlength=len(a["denom"]))
    for vid,a in accumulators.items():
        name = views.iloc[vid]["name"]
        for d,denominator in enumerate(a["denom"]):
            usable = denominator > 0 and d >= c["INSULATION_IGNORE_DIAGS"]
            rows.append({"region1":name,"region2":name,"dist":d,"distance_bp":d*clr.binsize,
                         "n_valid":int(denominator),"count.sum":a["raw"][d],"balanced.sum":a["balanced"][d],
                         "count.avg":a["raw"][d]/denominator if usable else np.nan,
                         "balanced.avg":a["balanced"][d]/denominator if usable else np.nan})
    return pd.DataFrame(rows, columns=EXPECTED_COLUMNS)


def fixture_map(tmp_path, square):
    bins = cooler.binnify(pd.Series({"contigA":2400,"contigB":1600}),400)
    # Deliberate absent/observable-zero pixels, diagonals, between-view cis and
    # trans contacts. Square storage also contains lower representatives.
    pixels = pd.DataFrame({"bin1_id":[0,0,0,0,1,1,2,3,4,4,6,6,7,8],
                           "bin2_id":[0,1,2,7,1,2,5,5,4,5,6,7,9,9],
                           "count":[9,2,3,101,5,7,103,107,11,13,17,19,23,29]})
    if square:
        lower = pixels.loc[pixels.bin1_id != pixels.bin2_id].rename(columns={"bin1_id":"bin2_id","bin2_id":"bin1_id"})[pixels.columns]
        pixels = pd.concat([pixels,lower],ignore_index=True).sort_values(["bin1_id","bin2_id"])
    path = tmp_path/("square.cool" if square else "upper.cool")
    cooler.create_cooler(str(path),bins,pixels,ordered=True,symmetric_upper=not square)
    views = pd.DataFrame([["contigA",0,1200,"left"],["contigA",1600,2400,"right"],
                          ["contigB",0,1600,"other"],["contigB",0,0,"empty"]], columns=["chrom","start","end","name"])
    labels = np.array([0,0,0,-1,1,1,2,2,2,2])
    return cooler.Cooler(str(path)),bins,views,labels,path


@pytest.mark.parametrize("square",[False,True])
@pytest.mark.parametrize("cutoff",[0,1,2,20])
@pytest.mark.parametrize("mask_pattern",["full","holes","fully_masked_view"])
def test_vectorized_expected_is_exactly_equal_across_storage_masks_and_views(tmp_path,config,square,cutoff,mask_pattern):
    clr,bins,views,labels,path = fixture_map(tmp_path,square)
    mask = np.ones(len(bins),dtype=bool)
    if mask_pattern == "holes":
        mask[[1,5,7]]=False
    elif mask_pattern == "fully_masked_view":
        mask[6:]=False
    values = np.linspace(.7,1.9,len(bins));values[~mask]=np.nan
    config.update(PIXEL_CHUNKSIZE=3,INSULATION_IGNORE_DIAGS=cutoff)
    before = file_hash(path)
    old = legacy_expected_table(clr,bins,mask,labels,views,values,config)
    new = expected_table(clr,bins,mask,labels,views,values,config)
    # Cooler scalar resolutions can infer int32 under Windows in the old
    # construction. Keep safe int64 distances while requiring identical values.
    pd.testing.assert_frame_equal(new,old,check_exact=True,check_dtype=False)
    assert new.dist.dtype==new.distance_bp.dtype==np.dtype("int64")
    assert new.groupby("region1",sort=False).size().to_dict()=={"left":3,"right":2,"other":4}
    assert not new.region1.eq("empty").any()
    if cutoff == 0 and mask_pattern == "full":
        distant = new.loc[(new.region1=="other") & (new.dist==3)].iloc[0]
        assert distant.n_valid==1 and distant["count.sum"]==distant["count.avg"]==0
    if mask_pattern == "fully_masked_view":
        missing = new.loc[new.region1=="other"]
        assert missing.n_valid.eq(0).all() and missing["count.avg"].isna().all()
    assert file_hash(path)==before


def test_no_registered_view_bins_preserves_empty_columns(tmp_path,config):
    clr,bins,views,labels,_ = fixture_map(tmp_path,False)
    labels[:]=-1
    args = (clr,bins,np.ones(len(bins),dtype=bool),labels,views,np.ones(len(bins)),config)
    pd.testing.assert_frame_equal(expected_table(*args),legacy_expected_table(*args),check_exact=True)
    assert list(expected_table(*args))==EXPECTED_COLUMNS
