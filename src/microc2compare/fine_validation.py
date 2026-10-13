"""Explicit admission for optional ultra-fine grids and sourced method claims."""
from pathlib import Path
from .contracts import file_hash, load_structured, require, resolve
from .resolutions import loop_resolutions, quantification_resolutions, chromosight_discovery_resolutions


def validate(c, samples, ref):
    from .analyses import needs_loop_counts, aggregate_modules
    active_counts = needs_loop_counts(c) or 'apa' in aggregate_modules(c)
    grids=sorted(set((loop_resolutions(c) if c['RUN_LOOPS'] else [])+(quantification_resolutions(c) if active_counts else [])+(chromosight_discovery_resolutions(c) if c.get('RUN_CHROMOSIGHT_DISCOVERY') else [])))
    active=c['RUN_LOOPS'] or c.get('RUN_CHROMOSIGHT_DISCOVERY') or active_counts
    if c['ANALYSIS_PROFILE']=='practical_2_4kb' and active:
        require(c['LOOP_CALLER']=='mustache' and {2000,4000}<=set(loop_resolutions(c)) and {2000,4000}<=set(quantification_resolutions(c)),
                'practical_2_4kb requires Mustache discovery and quantification at both 2000 and 4000 bp')
    if active and any(r<=1000 for r in grids):
        require(c['RESOLUTION_ELIGIBILITY_REGISTRY'], '1-kb/sub-kb analysis needs a reviewed RESOLUTION_ELIGIBILITY_REGISTRY; start with 2/4 kb')
        path=Path(c['RESOLUTION_ELIGIBILITY_REGISTRY']); registry=load_structured(path)
        require(registry.get('reference_id')==c['REFERENCE_ID'] and registry.get('scope_digest')==ref['effective_scope']['digest'], 'Fine-grid eligibility reference/scope mismatch')
        audit=resolve(registry.get('audit_path',''),path.parent)
        require(audit and registry.get('audit_sha256')==file_hash(audit), 'Fine-grid eligibility must bind a resolution audit by SHA256')
        evidence=load_structured(audit)
        require(evidence.get('reference_registry_sha256')==file_hash(c['REFERENCE_REGISTRY']), 'Fine-grid reference assets changed since audit')
        require(set(evidence['samples'])=={s['sample_id'] for s in samples}, 'Fine-grid audit sample membership mismatch')
        for res in (r for r in grids if r<=1000):
            decision=registry.get('resolutions',{}).get(str(res),{})
            require(decision.get('approved') is True and decision.get('reviewer') and decision.get('retention_review') and decision.get('reproducibility_review'), f'{res} bp needs explicit coverage/contact-retention and replicate-reproducibility review')
            rows=[r for r in evidence['rows'] if r['resolution']==res]
            require(len(rows)==len(samples) and {r['sample_id'] for r in rows}=={s['sample_id'] for s in samples} and all(r['status']=='metadata_usable' for r in rows), f'{res} bp audit contains missing/duplicate samples or unconverged inputs')
            for row in rows:
                from .matrix import matrix_uri
                from .resolutions import sample_matrix
                sample=next(s for s in samples if s['sample_id']==row['sample_id'])
                require(row['uri']==matrix_uri(sample_matrix(c,sample,res),res),'Fine-grid current input differs from reviewed audit')
                p=Path(row['uri'].split('::')[0]); stat=p.stat()
                require(row.get('file_size')==stat.st_size and row.get('file_mtime_ns')==stat.st_mtime_ns, 'Fine-grid input changed since eligibility audit')
    if c['ANALYSIS_PROFILE']=='hsieh_inspired':
        require(c['HSIEH_METHODS_REGISTRY'], 'hsieh_inspired requires HSIEH_METHODS_REGISTRY; practical_2_4kb remains an explicitly adapted profile')
        registry=load_structured(c['HSIEH_METHODS_REGISTRY'])
        require(registry.get('sources') and registry.get('differences_from_publication') and registry.get('verified_parameters'), 'Record primary sources, verified settings and explicit differences from the publication')
        for key, record in registry['verified_parameters'].items():
            require(key in c and record.get('source') and record.get('location') and record.get('value')==c[key], f'Hsieh parameter/source mismatch: {key}')
    return {'profile':c['ANALYSIS_PROFILE'],'resolutions':grids,'exact_publication_reproduction':False}
