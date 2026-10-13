from pathlib import Path
from importlib.resources import files
import numpy as np
import pandas as pd
import pytest
from microc2compare import science
from microc2compare.contracts import ContractError,atomic_json,file_hash
from microc2compare.chromosomes import chromosome_selection
from microc2compare.fine_loops import cross_resolution_catalogue,map_candidates,chromosight_scores
from microc2compare.features import LOOP_COLUMNS
from microc2compare.annotation import annotate_anchors,read_gtf,read_ccre
from microc2compare.state import scope_identity,initialize


@pytest.mark.parametrize('assembly,count',[('mm39',19),('hg38',22)])
def test_autosomes_default_and_sex_opt_in(assembly,count,config):
    names=[f'chr{i}' for i in range(1,count+1)]+['chrX','chrY','chrM']
    reference={'reference_id':assembly}
    assert chromosome_selection(config,reference,names,names)['chromosomes']==names[:count]
    config.update(CHROMOSOME_SCOPE='autosomes_plus_sex',SEX_CHROMOSOMES='chrX')
    assert chromosome_selection(config,reference,names,names)['chromosomes']==names[:count]+['chrX']
    with pytest.raises(ContractError):chromosome_selection(config,reference,names,names[:count])


def calls(res,starts):
    return pd.DataFrame([['chr1',a,a+res,'chr1',b,b+res,'source'] for a,b in starts],columns=LOOP_COLUMNS)


def test_complete_linkage_retains_source_and_avoids_transitive_clusters(config):
    table,members=cross_resolution_catalogue(config,[('A',2000,calls(2000,[(0,20000),(4000,24000)])),('B',2000,calls(2000,[(2000,22000)]))])
    assert len(table)==2 and len(members)==3
    assert table.source_count.sum()==3
    assert table.inference_status.eq('descriptive_same_data').all()
    shuffled,_=cross_resolution_catalogue(config,[('B',2000,calls(2000,[(2000,22000)])),('A',2000,calls(2000,[(4000,24000),(0,20000)]))])
    pd.testing.assert_frame_equal(table,shuffled)


def test_coarse_alias_collisions_are_one_measurement_and_finer_projections_labelled():
    candidates=calls(2000,[(0,20000),(2000,22000)]);candidates.feature_id=['a','b'];candidates['representative_resolution']=2000
    measured,aliases=map_candidates(candidates,'mm39',4000)
    assert len(measured)==1 and aliases.collapsed_candidate_count.eq(2).all()
    coarse=calls(4000,[(0,20000)]);coarse['representative_resolution']=4000
    fine,aliases=map_candidates(coarse,'mm39',2000)
    assert aliases.coarse_to_fine_projection.all() and fine.start1.iloc[0]==0


def test_real_chromosight_kernel_score_and_ineligible_candidate(tmp_path,config,monkeypatch):
    pytest.importorskip('chromosight.utils.detection')
    kernel=np.loadtxt(str(files('chromosight.kernels').joinpath('artificial_template_loops_type1.txt')))
    window=np.pad(kernel,1,constant_values=1.)
    class Map:
        chromsizes={'chr1':1000000}
    monkeypatch.setattr(science,'context',lambda *args:(Map(),None,None,None,None,None))
    monkeypatch.setattr(science,'expected_table',lambda *args:pd.DataFrame(columns=['region1','region2','dist','balanced.avg']))
    monkeypatch.setattr(science,'oe_window',lambda *args:window.copy())
    config.update(run_dir=str(tmp_path),LOOP_RESOLUTION=2000,LOOP_QUANTIFICATION_RESOLUTIONS='2000',RUN_LOOPS=True)
    table=calls(2000,[(100000,300000),(120000,320000)]);table.feature_id=['yes','excluded'];table['A']=[10,np.nan]
    science.save(table,tmp_path/'04_quantification/discovered_loops/2000/counts.tsv')
    products=chromosight_scores(config,'A');frame=pd.read_csv(products[0],sep='\t')
    assert frame.score.iloc[0]==pytest.approx(1.,abs=1e-10)
    assert frame.status.tolist()==['measured','ineligible_candidate']
    assert frame.PValue.isna().all() and frame.FDR.isna().all()


def test_both_strands_all_overlaps_and_ccre_evidence(tmp_path,config):
    gtf=tmp_path/'genes.gtf';gtf.write_text('chr1\tx\tgene\t1001\t5000\t.\t+\t.\tgene_id "A"; gene_name "Alpha";\nchr1\tx\texon\t1001\t1200\t.\t+\t.\tgene_id "A";\nchr1\tx\tgene\t8001\t10000\t.\t-\t.\tgene_id "B";\n')
    bed=tmp_path/'ccre.bed';bed.write_text('chr1\t0\t1500\te1\tdELS\nchr1\t500\t2000\te2\tpELS\n')
    genes,exons=read_gtf(gtf,['chr1']);ccre=read_ccre(bed,['chr1'],{})
    anchors=pd.DataFrame([['one','chr1',900,1100],['two','chr1',9900,10100],['three','chr1',20000,21000]],columns=['anchor_id','chrom','start','end'])
    result=annotate_anchors(anchors,genes,exons,ccre,config,True)
    assert result.gene_context.tolist()==['promoter','promoter','distal_intergenic']
    assert result.ccre_all_ids.iloc[0]=='e1;e2' and result.enhancer_like.iloc[0]
    assert result.promoter_gene_ids.iloc[1]=='B'
    assert result.nearest_tss_signed_distance_bp.iloc[2]<0


def test_effective_view_change_rejected_before_writing_saved_metadata(tmp_path,config):
    sizes=tmp_path/'sizes';sizes.write_text('chr1\t100000\n');views=tmp_path/'views';views.write_text('chr1\t0\t100000\tchr1\n')
    registry=tmp_path/'ref.json';atomic_json(registry,{'reference_id':'mm39','chromsizes':str(sizes),'views':str(views),'chromosome_classes':{'autosomes':['chr1']}})
    config.update(REFERENCE_ID='mm39',REFERENCE_REGISTRY=str(registry),run_dir=str(tmp_path/'run'))
    atomic_json(tmp_path/'run/00_metadata/config.json',config);atomic_json(tmp_path/'run/00_metadata/chromosome_policy.json',scope_identity(config))
    before=file_hash(tmp_path/'run/00_metadata/config.json');views.write_text('chr1\t0\t80000\tchr1\n')
    with pytest.raises(ContractError,match='Effective chromosome'):initialize(config,1,'test')
    assert before==file_hash(tmp_path/'run/00_metadata/config.json')


def test_complete_fixed_candidate_path_counts_scores_annotations_apa_and_resume(tmp_path,config,monkeypatch):
    pytest.importorskip('chromosight.utils.detection')
    import importlib.util,sys
    from microc2compare.contracts import normalize,enabled,STAGES
    from microc2compare.state import scopes,publish_completion,complete
    from microc2compare.report import render
    monkeypatch.setenv('MPLCONFIGDIR',str(tmp_path/'plot-cache'))
    source=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('fine_generator',source/'scripts/generate_fine_loop_demo.py')
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    data=tmp_path/'data';generator.generate(data,development=True)
    values={key:str(value).lower() if isinstance(value,bool) else str(value) for key,value in config.items()}
    values.update(PROJECT='demo',REFERENCE_ID='demo',SAMPLESHEET='samples.csv',TECHNICAL_MEMBERSHIP='technical_membership.tsv',REFERENCE_REGISTRY='references.json',INDEPENDENT_FEATURES_REGISTRY='independent.json',ANNOTATION_REGISTRY='annotation.json',OUTDIR='results',SCRATCH_ROOT='scratch',ENVIRONMENT_MODE='development',PYTHON_BOOTSTRAP=sys.executable,
                  LOOP_CALLER='mustache',LOOP_RESOLUTION='2000',LOOP_DISCOVERY_RESOLUTIONS='2000,4000',LOOP_QUANTIFICATION_RESOLUTIONS='2000,4000',LOOP_MIN_DISTANCE_BP='5000',LOOP_MAX_DISTANCE_BP='500000',ANALYSIS_PROFILE='practical_2_4kb',TOTAL_CPUS='4',MAX_MEMORY_MB='16000',M2C_ANALYSES='chromosight,annotation,apa')
    for kind in ('QUANTIFY','AGGREGATES','LOOPS','REPORT','PREFLIGHT'):values.update({f'THREADS_{kind}':'1',f'MEMORY_MB_{kind}':'4000'})
    c=normalize(values,data);initialize(c,1,'local-test')
    before={p.name:file_hash(p) for p in data.glob('*.mcool')}
    for stage in STAGES:
        if not enabled(c,stage):continue
        for unit in scopes(c,stage):
            products=render(c,unit) if stage=='report' else science.IMPLEMENTATIONS[stage](c,unit)
            publish_completion(c,stage,unit,products)
    root=Path(c['run_dir'])
    for res in (2000,4000):
        counts=pd.read_csv(root/f'04_quantification/independent_loops/{res}/counts.tsv',sep='\t')
        assert len(counts)==1 and counts[['C1','C2','C3','T1','T2','T3']].notna().all().all()
        for sid in ('C1','C2','C3','T1','T2','T3'):
            scores=pd.read_csv(root/f'04_quantification/chromosight/{res}/{sid}.tsv',sep='\t')
            assert scores.status.eq('measured').all()
            assert (root/f'08_reports/figures/apa.{res}.{sid}.svg').is_file()
    assert 'promoter_enhancer_like' in (root/'03_features/annotation/loops.tsv').read_text()
    assert 'APA aggregates' in (root/'08_reports/index.html').read_text()
    initialize(c,2,'local-test')
    assert all(complete(c,stage,unit) for stage in STAGES if enabled(c,stage) for unit in scopes(c,stage))
    assert before=={p.name:file_hash(p) for p in data.glob('*.mcool')}
    assert not list(root.rglob('*.cool')) and not list(root.rglob('*.mcool'))


def test_volcano_export_for_each_grid_without_changing_statistics(tmp_path,config,monkeypatch):
    from microc2compare.volcano import render_loop_volcanoes
    monkeypatch.setenv('MPLCONFIGDIR',str(tmp_path/'plot-cache'))
    config.update(run_dir=str(tmp_path),LOOP_RESOLUTION=2000,LOOP_QUANTIFICATION_RESOLUTIONS='2000,4000')
    comparison={'comparison_id':'case_vs_ref','modules':'loops','case_condition':'case','reference_condition':'ref','inference':'formal'}
    paths=[];hashes={}
    for res in (2000,4000):
        directory=tmp_path/f'06_comparisons/case_vs_ref/{res}';directory.mkdir(parents=True)
        table=directory/'loops.tsv';pd.DataFrame([{'feature_id':'a','logFC':1.,'PValue':.001,'FDR':.01,'retained':True,'inference_status':'formal','validation_status':'uncalibrated'}]).to_csv(table,sep='\t',index=False)
        params=directory/'parameters.json';atomic_json(params,{'inference':'formal','validation_status':'uncalibrated','fdr':.05})
        paths.extend([table,params]);hashes[str(table)]=file_hash(table)
    products,plots,status=render_loop_volcanoes(config,[comparison],paths)
    assert len(plots)==2 and {r['resolution'] for r in status}=={2000,4000}
    assert all(r['status']=='plotted' for r in status)
    assert {p.name for p in products if p.suffix=='.svg'}=={'volcano.case_vs_ref.2000.svg','volcano.case_vs_ref.4000.svg'}
    assert hashes=={str(p):file_hash(p) for p in paths if p.name=='loops.tsv'}


def test_explicit_preparation_conserves_raw_counts_and_keeps_source_unchanged(tmp_path):
    import cooler,h5py,subprocess,sys,json
    bins=cooler.binnify(pd.Series({'chr1':80000}),1000)
    i,j=np.triu_indices(len(bins));pixels=pd.DataFrame({'bin1_id':i,'bin2_id':j,'count':np.full(len(i),10,dtype=np.int64)})
    source=tmp_path/'source.cool';cooler.create_cooler(str(source),bins,pixels,ordered=True,assembly='mm39')
    recipe={'cis_only':False,'ignore_diags':0,'mad_max':5,'min_nnz':1,'min_count':0,'tol':1e-5}
    with h5py.File(source,'r+') as handle:
        ds=handle['bins'].create_dataset('weight',data=np.ones(len(bins)));ds.attrs.update(recipe);ds.attrs['converged']=True
    provenance=tmp_path/'source.json';atomic_json(provenance,{'reference_id':'mm39','balancing':{'1000':{'convention':'multiplicative','converged':True,'recipe':recipe}}})
    sheet=tmp_path/'samples.tsv';pd.DataFrame([{'sample_id':'one','reference_id':'mm39','matrix_path':str(source),'provenance_path':str(provenance)}]).to_csv(sheet,sep='\t',index=False)
    before=file_hash(source);destination=tmp_path/'prepared'
    script=Path(__file__).resolve().parents[1]/'scripts/prepare_resolutions.py'
    run=subprocess.run([sys.executable,'-I','-B',str(script),'--samples',str(sheet),'--out',str(destination),'--resolutions','2000','--threads','1'],text=True,capture_output=True,timeout=120)
    assert run.returncode==0,run.stdout+run.stderr
    derived=cooler.Cooler(str(destination/'one.2000.cool'))
    assert derived.info['sum']==cooler.Cooler(str(source)).info['sum']
    assert file_hash(source)==before
    with derived.open('r') as handle:assert bool(handle['bins/weight'].attrs['converged'])
    registry=json.loads((destination/'matrix_registry.json').read_text())
    assert registry['samples']['one']['2000']==str(destination/'one.2000.cool')


def test_one_kilobase_requires_reviewed_admission(config):
    from microc2compare.fine_validation import validate
    config.update(LOOP_RESOLUTION=1000,LOOP_DISCOVERY_RESOLUTIONS='1000',LOOP_QUANTIFICATION_RESOLUTIONS='1000',LOOP_CALLER='mustache')
    with pytest.raises(ContractError,match='reviewed RESOLUTION_ELIGIBILITY'):validate(config,[],{'effective_scope':{'digest':'test'}})


def test_effect_agreement_uses_the_recorded_descriptive_fallback_universe(tmp_path,config,monkeypatch):
    from microc2compare import fine_report
    monkeypatch.setattr(fine_report,'complete',lambda *args:True)
    config.update(run_dir=str(tmp_path),LOOP_RESOLUTION=2000,LOOP_QUANTIFICATION_RESOLUTIONS='2000,4000',
                  RUN_QC=False,RUN_LOOPS=False,RUN_AGGREGATES=False,RUN_CHROMOSIGHT=False,RUN_ANNOTATION=False,RUN_DIFFERENTIAL=True)
    comparison={'comparison_id':'case_vs_ref','modules':'loops','inference':'descriptive'}
    products=[]
    for resolution in (2000,4000):
        directory=tmp_path/f'06_comparisons/case_vs_ref/{resolution}'
        products.append(science.save(pd.DataFrame([{'feature_id':'measurement','log2_allocation_ratio':1.}]),directory/'loops.tsv'))
        params=directory/'parameters.json';atomic_json(params,{'loop_universe':'independent_loops'});products.append(params)
        products.append(science.save(pd.DataFrame([{'candidate_id':'original','measurement_id':'measurement','collapsed_candidate_count':1,'coarse_to_fine_projection':False}]),
                                    tmp_path/f'04_quantification/independent_loops/{resolution}/candidate_mapping.tsv'))
        # A stale discovery alias must never substitute for the recorded measurement universe.
        science.save(pd.DataFrame([{'candidate_id':'wrong','measurement_id':'measurement'}]),
                     tmp_path/f'04_quantification/discovered_loops/{resolution}/candidate_mapping.tsv')
    atomic_json(tmp_path/'completion/differential/run.done.json',{'stage':'differential','scope':'run','outputs':[{'path':str(p.relative_to(tmp_path))} for p in products]})
    report=tmp_path/'08_reports';report.mkdir()
    outputs,_,_=fine_report.render(config,[],{'comparisons':[comparison]})
    effects=pd.read_csv(report/'case_vs_ref.resolution_effects.tsv',sep='\t')
    assert effects.candidate_id.tolist()==['original','original']
    assert effects.resolution.tolist()==[2000,4000]
    assert (report/'case_vs_ref.resolution_effects.tsv') in outputs
