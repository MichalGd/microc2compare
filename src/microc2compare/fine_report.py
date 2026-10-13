"""Numerical-product-backed plots and explanations for fine-loop analysis."""
import html
from pathlib import Path
import numpy as np
import pandas as pd
from .contracts import IDENTIFIER,atomic_json,digest,read_json,require
from .state import complete
from .analyses import aggregate_modules,requires_saddle,has_loop_discovery
from .resolutions import quantification_resolutions


def render(c, samples, metadata):
    import matplotlib.pyplot as plt
    root=Path(c['run_dir']);report=root/'08_reports';figures=report/'figures'
    outputs,plots,sections=[],[],[]
    current=set()
    for record_path in (root/'completion').glob('*/*.done.json'):
        record=read_json(record_path)
        if record['stage']=='report':continue
        if complete(c,record['stage'],record['scope']):current.update((root/r['path']).resolve() for r in record['outputs'])
    contacts=[]
    for path in sorted((root/'01_qc/contacts').glob('*.json')):
        if path.resolve() not in current:continue
        import json
        record=json.loads(path.read_text())
        if 'raw_total' in record:contacts.append({'sample_id':path.stem,**record})
    if contacts and c['RUN_QC']:
        frame=pd.DataFrame(contacts)
        fields=[k for k in ('sample_id','resolution','raw_total','raw_cis','selected_view_cis','distance_eligible_cis_before_mask','eligible_cis_exposure','mask_excluded_distance_contacts') if k in frame]
        table=report/'contact_retention.tsv';frame[fields].to_csv(table,sep='\t',index=False);outputs.append(table)
        sections.append('<h2>Contact retention</h2><p>Native totals include every stored contact. Selected-view cis contacts have both ends in the same registered chromosome view. The distance restriction is applied next, then the common usable-bin mask. These nested denominators explain why eligible loop contacts are fewer than upstream valid pairs. They do not discard or change the imported files.</p>'+frame[fields].to_html(index=False))
    if has_loop_discovery(c):
        sections.append('<h2>Multi-resolution loop catalogue</h2><p>Source calls at each resolution remain available. Corresponding calls form one catalogue entry using complete-linkage matching at both anchors. Each resolution has one measurement per unique mapped pixel; candidate aliases identify collapsed entries. Mapping uses anchor midpoints. A coarse call projected onto a finer grid is a labelled measurement projection, not new evidence of finer localisation. Caller enrichment significance is separate from biological differential significance.</p>')
        if c['RUN_LOOPS'] and c['LOOP_CALLER']=='mustache':
            from .loop_profiles import mustache_geometry
            cutoffs=', '.join(f"{r} bp grid: best-case occupancy floor {mustache_geometry(c,r)['effective_best_case_min_distance_bp']} bp" for r in __import__('microc2compare.resolutions',fromlist=['loop_resolutions']).loop_resolutions(c))
            sections.append('<p>Mustache 1.3.3 suppresses its first diagonals: '+html.escape(cutoffs)+'. The requested measurement minimum is '+str(c['LOOP_MIN_DISTANCE_BP'])+' bp. Suppressed short-range calls are not proof of absent interactions. This is an adapted Hsieh-inspired analysis, not an exact publication reproduction.</p>')
    if c.get('RUN_STRIPES'):
        from .stripe_report import render as render_stripes
        extra_outputs, extra_plots, extra_sections = render_stripes(c, samples, current)
        outputs.extend(extra_outputs); plots.extend(extra_plots); sections.extend(extra_sections)
    if c.get('RUN_BROWSER_EXPORT'):
        browser_audits = [p for p in current if p.name == 'export_audit.json']
        if browser_audits:
            import os
            directory = sorted(browser_audits)[-1].parent
            relative = os.path.relpath(directory, report).replace('\\','/')
            sections.append('<h2>Browser exports</h2><p>UCSC custom tracks, per-observation tracks, condition summaries and coordinate profiles are available in <a href="'+html.escape(relative+'/README.md')+'">the browser README</a> and <a href="'+html.escape(relative+'/inventory.tsv')+'">track inventory</a>. HiGlass conversion and actual ingestion/rendering have separate statuses.</p>')
    for module in sorted(aggregate_modules(c)):
        grids=quantification_resolutions(c) if module=='apa' else [c['ADA_RESOLUTION'] if module=='ada' else c['COMPARTMENT_RESOLUTION']]
        for res in grids:
            directory=root/'05_aggregates'/module/(str(res) if module=='apa' else '')
            files=sorted(p for p in directory.glob('*.npz') if p.resolve() in current)
            condition_labels={}
            if module=='saddle':
                for condition in {s['condition'] for s in samples}:
                    identifier=condition if IDENTIFIER.fullmatch(condition) else 'condition_'+digest(condition)[:12]
                    condition_labels['condition__'+identifier]=condition+' (equal-replicate condition mean)'
                if requires_saddle(c):
                    expected={s['sample_id'] for s in samples}|set(condition_labels)
                    require(expected<={p.stem for p in files},
                            'The full report requires current saddle arrays for every biological replicate and condition; resume the aggregates stage first')
            arrays=[np.load(p)['mean'] for p in files]
            positive=np.concatenate([a[np.isfinite(a)&(a>0)] for a in arrays]) if arrays else np.array([])
            limit=max(.5,float(np.percentile(np.abs(np.log2(positive)),95))) if len(positive) else 1.
            details=[]
            for path,data in zip(files,arrays):
                with np.load(path) as archive:
                    coverage=next((archive[key] for key in ('valid_replicate_count','valid_pixel_count','valid_feature_count','valid_domain_count') if key in archive),np.isfinite(data).astype(int))
                good=np.isfinite(data)&(data>0)
                transformed=np.full(data.shape,np.nan);transformed[good]=np.log2(data[good])
                fig,ax=plt.subplots(figsize=(5,4))
                if good.any():
                    extent=[-c['APA_FLANK_BP']/1000,c['APA_FLANK_BP']/1000,c['APA_FLANK_BP']/1000,-c['APA_FLANK_BP']/1000] if module=='apa' else None
                    im=ax.imshow(transformed,cmap='coolwarm',vmin=-limit,vmax=limit,extent=extent)
                    fig.colorbar(im,ax=ax,label='log2 observed / expected')
                else:ax.text(.5,.5,'No finite positive aggregate\nReview coverage and candidate eligibility',ha='center',va='center',transform=ax.transAxes)
                axis='Anchor-centred position (kb)' if module=='apa' else 'Scaled domain position' if module=='ada' else 'Shared compartment quantile (B to A)'
                unit=condition_labels.get(path.stem,path.stem)
                ax.set(title=f'{module.upper()} {unit} ({res} bp)',xlabel=axis,ylabel=axis)
                fig.tight_layout()
                for suffix in ('png','svg'):
                    picture=figures/f'{module}.{res}.{path.stem}.{suffix}';fig.savefig(picture,dpi=150,bbox_inches='tight');outputs.append(picture)
                    if suffix=='png':plots.append((ax.get_title(),picture))
                plt.close(fig)
                details.append({'unit':path.stem,'resolution':res,'finite_pixels':int(np.isfinite(data).sum()),'positive_pixels':int(good.sum()),'total_pixels':data.size,'coverage_sum':int(coverage.sum()),'colour_limit_log2':limit,'zero_OE_not_visible_on_log_scale':int((data==0).sum())})
            if details:
                info=report/f'{module}.{res}.plot_metadata.json';atomic_json(info,details);outputs.append(info)
        explanation={'apa':'APA averages contact patterns centred on fixed loop anchors. The centre measures focal enrichment; neighbouring pixels show its context. All candidates are measured in every biological observation, including samples with no individual discovery call.',
                     'ada':'ADA averages domains after scaling them to a common shape. It complements the separately configured insulation/TAD calls; it does not change their grid or establish event significance.',
                     'saddle':'Saddle plots group contacts by a shared oriented compartment ordering. B–B is at the upper left, A–A at the lower right, and the other corners show mixed A–B contacts. Red means higher than expected and blue lower than expected. Blank pixels lack a finite positive value on the log scale; use the numerical arrays and coverage counts to distinguish observed zeros from unavailable measurements. The per-replicate segregation ratio is (AA + BB) / (2 × AB). These summaries are descriptive and carry no differential P value; orientation and phasing quality are required.'}[module]
        heading='Compartment saddle plots' if module=='saddle' else module.upper()+' aggregates'
        sections.append('<h2>'+heading+'</h2><p>'+explanation+' Plots show log2 observed/expected with a shared colour scale within each module and resolution. Condition summaries weight biological replicates equally. NPZ files preserve linear values, observed zeros and contribution counts; missing pixels are unavailable. Condition means are descriptive, not additional replicates.</p>')
        if module=='saddle':
            summary=root/'05_aggregates/saddle/summary.tsv'
            if summary.resolve() in current:sections.append(pd.read_csv(summary,sep='\t').to_html(index=False))
    if c['RUN_CHROMOSIGHT']:
        sections.append('<h2>Fixed-candidate Chromosight scores</h2><p>These scores measure similarity to the pinned Chromosight loop kernel in bounded observed/expected windows. Masked, flat and edge windows are flagged. This adapted scorer produces descriptive correlations, not biological P values or FDR. Every eligible fixed candidate is scored per replicate, regardless of whether that replicate called it.</p>')
        for res in quantification_resolutions(c):
            paths=sorted(p for p in (root/f'04_quantification/chromosight/{res}').glob('*.tsv') if p.resolve() in current)
            frames=[pd.read_csv(p,sep='\t') for p in paths]
            if frames:
                all_scores=pd.concat(frames);fig,ax=plt.subplots(figsize=(7,4))
                for sid,frame in all_scores.groupby('sample_id'):
                    values=frame.loc[frame.status.eq('measured'),'score'].dropna()
                    if len(values):ax.hist(values,bins=np.linspace(-1,1,31),histtype='step',label=sid)
                ax.set(xlabel='Kernel correlation',ylabel='Measured candidate pixels',title=f'Chromosight ({res} bp)')
                if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=7)
                fig.tight_layout();picture=figures/f'chromosight.{res}.png';fig.savefig(picture,dpi=150);plt.close(fig);outputs.append(picture);plots.append((f'Chromosight scores at {res} bp',picture))
                counts=all_scores.groupby(['sample_id','status']).size().reset_index(name='candidates');sections.append(counts.to_html(index=False))
    if c['RUN_ANNOTATION']:
        path=root/'03_features/annotation/loops.tsv'
        sections.append('<h2>Loop-anchor annotation</h2><p>Both anchors retain promoter, exon, gene-body and cCRE overlaps, ambiguous matches and nearest-TSS distances. Promoter–promoter or enhancer-like–promoter labels describe annotation evidence; they do not establish regulation, direction or cell-type-specific enhancer activity. cCRE provenance and any assembly lift are recorded.</p>')
        if path.resolve() in current:
            frame=pd.read_csv(path,sep='\t');column='candidate_category'
            if column in frame:sections.append(frame[column].value_counts().rename_axis('category').reset_index(name='candidates').to_html(index=False))
    # Match the same catalogue IDs between resolution-specific effect tables.
    for comparison in metadata.get('comparisons',[]):
        if not c['RUN_DIFFERENTIAL'] or 'loops' not in comparison['modules'].split(';'):continue
        tables=[]
        for res in quantification_resolutions(c):
            result=root/f"06_comparisons/{comparison['comparison_id']}/{res}/loops.tsv"
            kind='discovered_loops' if comparison['inference']=='descriptive' else 'independent_loops'
            parameters=result.parent/'parameters.json'
            if parameters.resolve() in current:
                kind=read_json(parameters).get('loop_universe',kind)
            if kind not in ('discovered_loops','independent_loops'):continue
            mapping=root/f'04_quantification/{kind}/{res}/candidate_mapping.tsv'
            if result.resolve() not in current or mapping.resolve() not in current:continue
            effects=pd.read_csv(result,sep='\t');aliases=pd.read_csv(mapping,sep='\t')
            column='logFC' if 'logFC' in effects else 'log2_allocation_ratio'
            if column not in effects:continue
            table=aliases.merge(effects[['feature_id',column]],left_on='measurement_id',right_on='feature_id',validate='many_to_one')
            table['resolution']=res;table=table.rename(columns={column:'effect'});tables.append(table)
        if len(tables)>1:
            long=pd.concat(tables,ignore_index=True);path=report/f"{comparison['comparison_id']}.resolution_effects.tsv";long.to_csv(path,sep='\t',index=False);outputs.append(path)
            sections.append('<h2>Resolution effect agreement</h2><p>The downloadable table joins effects through original candidate IDs. Alias collisions and coarse-to-fine projections remain labelled. Resolutions are correlated sensitivity measurements; their P values are not independent confirmations and are not pooled.</p>')
    return outputs,plots,sections
