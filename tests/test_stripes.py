"""Adapter tests and an optional actual pinned-StripePy CPU acceptance test."""
from pathlib import Path
import importlib.metadata
import json
import os
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import cooler
import numpy as np
import pandas as pd
import pytest

from microc2compare import stripes as module
from microc2compare.contracts import atomic_json, file_hash


def config(tmp_path):
    return dict(run_dir=str(tmp_path/"run"), SCRATCH_ROOT=str(tmp_path/"scratch"), REFERENCE_ID="newAssembly",
                STRIPE_RESOLUTIONS="4000", STRIPE_NORMALIZATION="none", STRIPE_MASK_POLICY="common_segments",
                STRIPE_GENOMIC_BELT_BP=5000000, STRIPE_MAX_WIDTH_BP=100000, STRIPE_MIN_RELATIVE_CHANGE=0,
                THREADS_STRIPES=1, PIXEL_CHUNKSIZE=1000)


def native(**changes):
    row = dict(seed=3, left_bound=3, right_bound=3, top_bound=3, bottom_bound=7,
               rel_change=100, inner_mean=2, outer_mean=1, top_persistence=.2, cfx_of_variation=.4)
    row.update(changes)
    return row


def segment():
    return dict(segment_id="segment000000", chrom="customChr", start=8000, end=48000, view="chosen", bin_count=10)


def observation(sid="fresh", condition="new condition"):
    return dict(sample_id=sid, condition=condition)


def test_inclusive_geometry_orientation_and_enrichment(tmp_path):
    c = config(tmp_path)
    row, reason = module.canonicalize(native(), segment(), "LT", observation(), 4000, c)
    assert reason is None
    assert (row["anchor_start"], row["anchor_end"], row["body_start"], row["body_end"]) == (20000,24000,20000,40000)
    assert row["width_bp"] == 4000 and row["direction"] == "downstream"
    upstream, _ = module.canonicalize(native(top_bound=0,bottom_bound=3), segment(), "UT", observation(), 4000, c)
    assert upstream["direction"] == "upstream" and upstream["body_end"] == row["anchor_end"]
    assert module.canonicalize(native(inner_mean=.5), segment(), "LT", observation(),4000,c)[1] == "nonpositive_contrast"
    assert module.canonicalize(native(outer_mean=0), segment(), "LT", observation(),4000,c)[1] == "undefined_contrast"
    c["STRIPE_MIN_RELATIVE_CHANGE"] = 101
    assert module.canonicalize(native(), segment(), "LT", observation(),4000,c)[1] == "below_relative_change_threshold"
    with pytest.raises(ValueError, match="direction/seed"):
        module.canonicalize(native(), segment(), "UT", observation(),4000,c)


def fixture(tmp_path, monkeypatch, n=100):
    c = config(tmp_path); bins = cooler.binnify(pd.Series({"customChr": n*4000}), 4000)
    i,j = np.triu_indices(n)
    count = np.maximum(1, 200//(j-i+1)).astype(np.int64)
    count[(i==25)&(j>25)] += 100
    path = tmp_path/"fresh.cool"
    cooler.create_cooler(str(path),bins,pd.DataFrame(dict(bin1_id=i,bin2_id=j,count=count)),ordered=True,assembly="newAssembly")
    clr = cooler.Cooler(str(path)); mask=np.ones(n,bool);mask[n//2]=False
    labels = np.zeros(n,dtype=int); views=pd.DataFrame([dict(chrom="customChr",start=0,end=n*4000,name="chosen")])
    from microc2compare import science
    monkeypatch.setattr(science,"context",lambda *args:(clr,bins,mask,labels,views,np.ones(n)))
    monkeypatch.setattr(module,"validate_metadata",lambda c:{"samples":[observation()]})
    meta=Path(c["run_dir"])/"00_metadata/masks";meta.mkdir(parents=True)
    np.savez(meta/"4000.npz",mask=mask)
    views.to_csv(meta/"4000.views.tsv",sep="\t",index=False)
    atomic_json(meta.parent/"inputs.json",{str(path):{"sha256":file_hash(path)}})
    return c,clr,bins,mask,labels,views,path


def test_segment_input_counts_mask_and_source_identity(tmp_path, monkeypatch):
    c,clr,bins,mask,labels,views,source=fixture(tmp_path,monkeypatch)
    segments=module.segments_for_context(bins,mask,labels,views,4000)
    assert segments.bin_count.tolist()==[50,49]
    before=file_hash(source);temporary=tmp_path/"tool.cool"
    stats=module.write_segment_input(clr,bins,segments,temporary,c)
    result=cooler.Cooler(str(temporary)); pixels=result.pixels()[:];target_bins=result.bins()[:]
    assert (target_bins.chrom.iloc[pixels.bin1_id].to_numpy()==target_bins.chrom.iloc[pixels.bin2_id].to_numpy()).all()
    original=clr.pixels()[:]
    same=(original.bin1_id<50)&(original.bin2_id<50) | (original.bin1_id>50)&(original.bin2_id>50)
    assert stats["raw_contact_sum"]==int(original.loc[same,"count"].sum())==int(result.info["sum"])
    assert len(target_bins)==99 and before==file_hash(source)


def test_missing_version_actionable_and_wrong_version_rejected(monkeypatch):
    monkeypatch.setattr(importlib.metadata,"version",lambda name:"9.0")
    with pytest.raises(ValueError,match="requires stripepy-hic==1.3.0"):module.probe()


def test_fake_caller_products_are_explicit_and_cleanup(tmp_path, monkeypatch):
    c,_,_,_,_,_,source=fixture(tmp_path,monkeypatch)
    monkeypatch.setattr(module,"probe",lambda:dict(software="stripepy-hic",version="1.3.0",normalization="raw_counts",cpu_only=True))
    commands=[]
    def run(command,stdout_path=None):
        commands.append(command)
        if stdout_path:
            Path(stdout_path).write_text("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\n")
        else:
            Path(command[command.index("--output-file")+1]).write_bytes(b"fake API fixture; not a real caller")
    monkeypatch.setattr(module,"run_caller",run)
    class FakeResult:
        chromosomes={"segment000000":200000,"segment000001":196000}
        resolution=4000
        def __init__(self,path):pass
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def get(self,chrom,field,location):
            return pd.DataFrame([native()]) if chrom=="segment000000" and location=="LT" else pd.DataFrame()
    monkeypatch.setitem(sys.modules,"stripepy.data_structures",SimpleNamespace(ResultFile=FakeResult))
    before=file_hash(source);outputs=module.stripes(c,"fresh")
    assert all(p.exists() and p.stat().st_size for p in outputs)
    table=pd.read_csv(Path(c["run_dir"])/"03_features/stripes/fresh.4000.tsv",sep="\t")
    assert table.direction.tolist()==["downstream"] and not {"PValue","FDR"}&set(table.columns)
    assert file_hash(source)==before and not list(Path(c["SCRATCH_ROOT"]).rglob("*.cool"))
    assert commands[0][commands[0].index("--min-chrom-size")+1]=="1"
    assert commands[0][commands[0].index("--nproc")+1]=="1"


def test_summary_arbitrary_conditions_and_many_to_many_links(tmp_path,monkeypatch):
    c=config(tmp_path);directory=Path(c["run_dir"])/"03_features/stripes";directory.mkdir(parents=True)
    samples=[observation("A","alpha"),observation("B","alpha"),observation("C","other"),observation("D","single")]
    monkeypatch.setattr(module,"validate_metadata",lambda c:{"samples":samples})
    for sample in samples:
        row,_=module.canonicalize(native(),segment(),"LT",sample,4000,c)
        pd.DataFrame([] if sample["sample_id"]=="D" else [row],columns=module.STRIPE_COLUMNS).to_csv(directory/f"{sample['sample_id']}.4000.tsv",sep="\t",index=False)
    outputs=module.stripes_summary(c)
    assert len(outputs)==3
    summary=pd.read_csv(directory/"summary.tsv",sep="\t")
    assert summary.biological_replicates_in_condition.tolist()==[2,2,1,1]
    assert summary.stripe_count.tolist()==[1,1,1,0]
    links=pd.read_csv(directory/"matches.tsv",sep="\t")
    assert len(links)==3 and links.body_reciprocal_overlap.eq(1).all()


def has_real_caller():
    try:return importlib.metadata.version("stripepy-hic")==module.SUPPORTED_VERSION
    except importlib.metadata.PackageNotFoundError:return False


@pytest.mark.skipif(os.name != "posix", reason="AF_UNIX multiprocessing paths are POSIX-specific")
@pytest.mark.parametrize("inherited", ["short", "deep_ascii", "deep_unicode"])
def test_real_manager_socket_under_inherited_temp_and_private_cleanup(tmp_path, monkeypatch, inherited):
    # A short explicitly inherited root must be honored. Long and multibyte
    # roots exercise the pipeline failure while leaving parent scratch intact.
    with tempfile.TemporaryDirectory(prefix="m2c-test-", dir="/tmp") as owned:
        root = Path(owned)
        if inherited != "short":
            root /= "stage-scratch-" * 9 if inherited == "deep_ascii" else "\u03bb" * 45
            root.mkdir()
        for key in ("TMPDIR", "TMP", "TEMP"):
            monkeypatch.setenv(key, str(root))
        monkeypatch.setattr(tempfile, "tempdir", str(root))
        output = tmp_path / "ipc.json"
        script = (
            "import json,multiprocessing as mp,os,stat,tempfile\n"
            "if __name__ == '__main__':\n"
            " with mp.Manager() as manager:\n"
            "  value=manager.dict(answer=42)\n"
            "  print(json.dumps(dict(tmp=tempfile.gettempdir(),"
            "variables=[os.environ[key] for key in ('TMPDIR','TMP','TEMP')],"
            "mode=stat.S_IMODE(os.stat(tempfile.gettempdir()).st_mode),"
            "address=manager.address,answer=value['answer'])))\n"
        )
        module.run_caller([sys.executable, "-c", script], output)
        value = json.loads(output.read_text())
        temporary = Path(value["tmp"])
        assert value["answer"] == 42 and value["mode"] == 0o700
        assert value["variables"] == [str(temporary)] * 3
        assert len(os.fsencode(value["address"])) <= 103
        assert temporary.parent == (root.resolve() if inherited == "short" else Path("/tmp").resolve())
        assert not temporary.exists() and root.is_dir()
        assert tempfile.tempdir == str(root)
        assert all(os.environ[key] == str(root) for key in ("TMPDIR", "TMP", "TEMP"))


@pytest.mark.skipif(os.name != "posix", reason="Private POSIX caller temporary directory")
def test_caller_private_temp_cleanup_on_failure(tmp_path, monkeypatch):
    inherited = tmp_path / ("stage-scratch-" * 9)
    inherited.mkdir()
    monkeypatch.setenv("TMPDIR", str(inherited))
    output = tmp_path / "failed-ipc.json"
    with pytest.raises(subprocess.CalledProcessError):
        module.run_caller([sys.executable, "-c",
            "import json,os; print(json.dumps(os.environ['TMPDIR'])); raise SystemExit(17)"], output)
    temporary = Path(json.loads(output.read_text()))
    assert not temporary.exists() and inherited.is_dir()
    assert os.environ["TMPDIR"] == str(inherited)


@pytest.mark.skipif(not has_real_caller(),reason="Optional pinned StripePy environment is not installed")
@pytest.mark.parametrize("threads", [1, 2])
def test_actual_stripepy_cli_on_fresh_masked_map(tmp_path,monkeypatch,threads):
    c,_,_,_,_,_,source=fixture(tmp_path,monkeypatch,n=160)
    c["THREADS_STRIPES"] = threads
    if os.name == "posix":
        inherited = tmp_path / ("stage-scratch-" * 9)
        inherited.mkdir()
        for key in ("TMPDIR", "TMP", "TEMP"):
            monkeypatch.setenv(key, str(inherited))
    before=file_hash(source)
    outputs=module.stripes(c,"fresh")
    assert all(p.exists() and p.stat().st_size for p in outputs)
    assert file_hash(source)==before
    assert not list(Path(c["SCRATCH_ROOT"]).rglob("*.cool"))
    calls=pd.read_csv(Path(c["run_dir"])/"03_features/stripes/fresh.4000.tsv",sep="\t")
    assert len(calls) >= 1
    assert (calls.width_bp>0).all() and (calls.inner_mean>calls.outer_mean).all()
    assert not ((calls.body_start<320000)&(calls.body_end>320000)).any()
