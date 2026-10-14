import importlib.util
import sys
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import atomic_json, file_hash, normalize, read_json, validate_metadata
from microc2compare.science import aggregates, register
from microc2compare.state import initialize

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("module",["apa","ada","saddle"])
def test_separate_aggregate_outputs_and_read_only_maps(tmp_path,module,monkeypatch):
    spec=importlib.util.spec_from_file_location("generator",ROOT/"scripts/generate_demo.py")
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    data=tmp_path/"data";generator.generate(data,development=True)
    schema=read_json(ROOT/"schemas/config.schema.json")
    values={k:str(v["default"]).lower() if isinstance(v["default"],bool) else str(v["default"]) for k,v in schema["properties"].items()}
    values.update(PROJECT="demo",REFERENCE_ID="demo",SAMPLESHEET="samples.csv",REFERENCE_REGISTRY="references.json",
                  OUTDIR="results",SCRATCH_ROOT="scratch",RUN_DIFFERENTIAL="false",RUN_HICREP="false",RUN_QC="false",
                  RUN_INSULATION="true" if module=="ada" else "false",RUN_TADS="false",RUN_LOOPS="false",
                  RUN_COMPARTMENTS="true" if module=="saddle" else "false",RUN_AGGREGATES="true",AGGREGATE_MODULES=module)
    c=normalize(values,data);root=Path(c["run_dir"])
    if module=="saddle":
        for path in data.glob("*.mcool"):
            with h5py.File(path,"r+") as handle:del handle["resolutions/10000"]
    before={p.name:file_hash(p) for p in data.glob("*.mcool")}
    initialize(c,1,"aggregate-test")
    register(c,"run")
    if module=="ada":
        path=root/"03_features/boundary_derived_domains/consensus.tsv";path.parent.mkdir(parents=True)
        pd.DataFrame([dict(chrom="chrDemo",start=1000000,end=1500000,window_bp=100000,feature_id="domain")]).to_csv(path,sep="\t",index=False)
    if module=="saddle":
        path=root/"03_features/compartments";path.mkdir(parents=True)
        bins=pd.read_csv(root/"00_metadata/masks/100000.bins.tsv",sep="\t")
        for sample in validate_metadata(c)["samples"]:
            frame=bins.copy();frame["oriented_score"]=np.linspace(-1,1,len(frame))
            frame.to_csv(path/f"{sample['sample_id']}.tsv",sep="\t",index=False)
    products=aggregates(c,"run")
    assert all(p.is_file() and p.stat().st_size for p in products)
    assert {p.name for p in (root/"05_aggregates").iterdir()}=={module}
    files=list((root/f"05_aggregates/{module}").glob("*.npz"))
    assert len([p for p in files if not p.name.startswith('condition__')])==6
    if module!='apa':
        assert len([p for p in files if p.name.startswith('condition__')])==2
        for condition in ('control','treated'):
            members=['C1','C2','C3'] if condition=='control' else ['T1','T2','T3']
            stack=np.stack([np.load(root/f'05_aggregates/{module}/{sid}.npz')['mean'] for sid in members])
            with np.load(root/f'05_aggregates/{module}/condition__{condition}.npz') as aggregate:
                np.testing.assert_allclose(aggregate['mean'],np.nanmean(stack,axis=0),equal_nan=True)
    if module=="ada":
        scores=pd.read_csv(root/"04_quantification/domain_scores/strength.tsv",sep="\t")
        assert len(scores)==6 and scores.mean_intra_domain_oe.gt(0).all()
    if module=="saddle":
        # Real aggregate arrays from known oriented-score fixtures exercise the
        # HTML consumer. This does not substitute for actual cooltools eigensolving.
        from microc2compare import report, fine_report
        monkeypatch.setenv('MPLCONFIGDIR',str(tmp_path/'plot-cache'))
        monkeypatch.setattr(report,'complete',lambda c,stage,scope:stage=='aggregates')
        monkeypatch.setattr(fine_report,'complete',lambda c,stage,scope:stage=='aggregates')
        atomic_json(root/'completion/run/aggregates.done.json',{'stage':'aggregates','scope':'run','outputs':[{'path':str(p.relative_to(root))} for p in products]})
        before_arrays={str(p):file_hash(p) for p in files}
        c['ANALYSIS_PROFILE']='practical_2_4kb'
        report.render(c,'run')
        page=(root/'08_reports/index.html').read_text(encoding='utf-8')
        assert 'Compartment saddle plots' in page and 'segregation_ratio' in page
        for unit in ('C1','C2','C3','T1','T2','T3','condition__control','condition__treated'):
            assert f'src="figures/saddle.100000.{unit}.png"' in page
            assert f'saddle.100000.{unit}.svg' in page
        assert before_arrays=={str(p):file_hash(p) for p in files}
        (root/'05_aggregates/saddle/condition__control.npz').unlink()
        with pytest.raises(ValueError,match='full report requires current saddle arrays'):
            report.render(c,'run')
    assert before=={p.name:file_hash(p) for p in data.glob("*.mcool")}
    assert not list(root.rglob("*.cool")) and not list(root.rglob("*.mcool"))
