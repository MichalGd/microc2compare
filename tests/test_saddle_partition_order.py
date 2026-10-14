"""Exact saddle partitions are independent of caller observation-list order."""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata

from microc2compare import science
from microc2compare.contracts import ContractError


def partition_fixture(tmp_path, config, monkeypatch, tied=False):
    c=dict(config,run_dir=str(tmp_path),COMPARTMENT_RESOLUTION=10000,
           SADDLE_PARTITION_SOURCE='replicate_rank_mean')
    # Manifest order intentionally differs from alphabetic observation order.
    samples=[{'sample_id':sid} for sid in ('zeta','alpha','eta','beta','delta','gamma')]
    monkeypatch.setattr(science,'validate_metadata',lambda c:{'samples':samples})
    base=pd.DataFrame({'chrom':['chrSynthetic']*200+['chrShort']*45,
                       'start':list(np.arange(200)*10000)+list(np.arange(45)*10000),
                       'end':list(np.arange(1,201)*10000)+list(np.arange(1,46)*10000)})
    labels=np.r_[np.zeros(200,dtype=int),np.ones(45,dtype=int)]
    mask=np.ones(len(base),dtype=bool);mask[7]=False
    usable=np.setdiff1d(np.arange(200),np.r_[7,np.arange(20,26)])
    directory=tmp_path/'00_metadata/masks';directory.mkdir(parents=True)
    base.to_csv(directory/'10000.bins.tsv',sep='\t',index=False)
    np.savez(directory/'10000.npz',labels=labels,mask=mask)
    directory=tmp_path/'03_features/compartments';directory.mkdir(parents=True)
    rng=np.random.default_rng(0)
    ranks=[]
    for index,sample in enumerate(samples):
        table=base.copy();values=np.ones(len(base))*1000
        values[usable]=4 if tied else rng.permutation(len(usable))
        values[20+index]=np.nan
        table['oriented_score']=values
        table.to_csv(directory/(sample['sample_id']+'.tsv'),sep='\t',index=False)
        ranks.append(rankdata(values[usable],method='average')/len(usable))
    # This is the unchanged normal manifest-order rank mean, before sorting.
    expected_rank=np.stack(ranks).mean(axis=0)
    expected_partition=np.full(len(base),-1,dtype=int)
    ordered=usable[np.lexsort((usable,expected_rank))]
    expected_partition[ordered]=np.arange(len(usable))*50//len(usable)
    return c,samples,usable,expected_rank,expected_partition


@pytest.mark.parametrize('tied',[False,True])
def test_saddle_partition_manifest_order_is_exact_under_permutations(tmp_path,config,monkeypatch,tied):
    c,samples,usable,rank,partition=partition_fixture(tmp_path,config,monkeypatch,tied)
    expected,labels=science.saddle_partition(c,samples)
    np.testing.assert_array_equal(expected.partition,partition)
    np.testing.assert_array_equal(expected.reference_rank.to_numpy()[usable],rank)
    # A masked bin, each observation's unavailable bin and a short view stay out.
    assert expected.loc[np.r_[7,np.arange(20,26),np.arange(200,245)],'partition'].eq(-1).all()
    assert expected.loc[np.r_[7,np.arange(20,26),np.arange(200,245)],'reference_rank'].isna().all()
    for changed in (samples[::-1],samples[2:]+samples[:2],sorted(samples,key=lambda s:s['sample_id'])):
        actual,again=science.saddle_partition(c,changed)
        pd.testing.assert_frame_equal(actual,expected,check_exact=True)
        np.testing.assert_array_equal(again,labels)
    if tied:
        assert np.ptp(expected.reference_rank.to_numpy()[usable])==0
        assert expected.partition.to_numpy()[usable].tolist()==partition[usable].tolist()


def test_saddle_partition_keeps_subset_in_manifest_order(tmp_path,config,monkeypatch):
    c,samples,_,_,_=partition_fixture(tmp_path,config,monkeypatch)
    first,_=science.saddle_partition(c,samples[::2])
    again,_=science.saddle_partition(c,samples[::2][::-1])
    pd.testing.assert_frame_equal(first,again,check_exact=True)


@pytest.mark.parametrize('invalid',['duplicate','unknown'])
def test_saddle_partition_rejects_invalid_observation_identity(tmp_path,config,monkeypatch,invalid):
    c,samples,_,_,_=partition_fixture(tmp_path,config,monkeypatch)
    requested=samples+[samples[0]] if invalid=='duplicate' else [{'sample_id':'outside_manifest'}]
    with pytest.raises(ContractError,match='Saddle observation'):
        science.saddle_partition(c,requested)
