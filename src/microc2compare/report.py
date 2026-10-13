"""Standalone HTML and bounded plots; scientific results remain downloadable."""
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .contracts import atomic_json, enabled, read_json, validate_metadata
from .state import STAGES, complete, scopes
from .analyses import aggregate_modules


def render(c,scope):
    root=Path(c["run_dir"]);directory=root/"08_reports";directory.mkdir(parents=True,exist_ok=True)
    figures=directory/"figures";figures.mkdir(exist_ok=True)
    outputs=[];plots=[]
    tad_tables=[]
    metadata=validate_metadata(c)
    samples=metadata["samples"]
    volcano_statuses=[]
    if c["RUN_DIFFERENTIAL"]:
        from .volcano import render_loop_volcanoes
        record=root/"completion/run/differential.done.json"
        current=[root/entry["path"] for entry in read_json(record)["outputs"]] if record.exists() and complete(c,"differential","run") else []
        products,pictures,volcano_statuses=render_loop_volcanoes(c,metadata.get("comparisons",[]),current)
        outputs.extend(products);plots.extend(pictures)
    statuses=[]
    for stage in STAGES[:-1]:
        if not enabled(c,stage):status="skipped_config"
        else:status="validated" if all(complete(c,stage,s) for s in scopes(c,stage)) else "missing_or_stale"
        statuses.append({"stage":stage,"status":status})
    statuses.append({"stage":"hicrep","status":"validated" if c["RUN_HICREP"] and complete(c,"hicrep","run") else "missing_or_stale" if c["RUN_HICREP"] else "skipped_config"})
    if c["RUN_QC"]:
        fig,ax=plt.subplots(figsize=(7,4))
        for sample in samples:
            path=root/f"01_qc/ps/{sample['sample_id']}.{c['PS_RESOLUTION']}.tsv"
            if not path.exists():continue
            frame=pd.read_csv(path,sep="\t")
            grouped=frame.groupby("distance_bp")[["balanced.sum","n_valid"]].sum()
            signal=grouped["balanced.sum"]/grouped.n_valid
            good=(grouped.index>0)&np.isfinite(signal)&(signal>0)
            ax.loglog(grouped.index[good],signal[good],label=sample["sample_id"])
        ax.set(xlabel="Distance (bp)",ylabel="Balanced contacts per eligible pixel")
        ax.legend(fontsize=7);fig.tight_layout()
        path=figures/"ps.png";fig.savefig(path,dpi=150);plt.close(fig);outputs.append(path);plots.append(("Distance decay",path))
    summary=root/"01_qc/hicrep/summary.tsv"
    if c["RUN_HICREP"] and summary.exists():
        frame=pd.read_csv(summary,sep="\t")
        if len(frame):
            ids=[s["sample_id"] for s in samples];matrix=np.full((len(ids),len(ids)),np.nan);np.fill_diagonal(matrix,1)
            for row in frame.itertuples():matrix[ids.index(row.sample1),ids.index(row.sample2)]=matrix[ids.index(row.sample2),ids.index(row.sample1)]=row.weighted_scc
            fig,ax=plt.subplots(figsize=(6,5));im=ax.imshow(matrix,vmin=-1,vmax=1,cmap="coolwarm")
            ax.set_xticks(range(len(ids)),ids,rotation=90);ax.set_yticks(range(len(ids)),ids);fig.colorbar(im,ax=ax,label="HiCRep SCC")
            fig.tight_layout();path=figures/"hicrep.png";fig.savefig(path,dpi=150);plt.close(fig);outputs.append(path);plots.append(("Replicate concordance",path))
    from .fine_report import render as render_fine
    products,pictures,fine_sections=render_fine(c,samples,metadata)
    outputs.extend(products);plots.extend(pictures)
    # One bounded representative locus per replicate; read original counts/weights.
    if enabled(c,"tads"):
        tad_record=root/"completion/run/tads.done.json"
        tad_paths=[root/entry["path"] for entry in read_json(tad_record)["outputs"]] if tad_record.exists() and complete(c,"tads","run") else []
        for path in sorted(tad_paths):
            if not path.name.startswith("tad_events.") or path.suffix!=".tsv":continue
            if path.name.endswith(".summary.tsv"): continue
            frame=pd.read_csv(path,sep="\t")
            if len(frame): tad_tables.append(frame)
        if tad_tables:
            events=pd.concat(tad_tables,ignore_index=True)
            primary=events.loc[events.window_bp==c["INSULATION_PRIMARY_WINDOW_BP"]]
            if len(primary):
                counts=primary.event.value_counts()
                fig,ax=plt.subplots(figsize=(8,4));counts.plot.bar(ax=ax)
                ax.set(xlabel="Descriptive TAD event",ylabel="Event groups (all declared comparisons)");fig.tight_layout()
                path=figures/"tad_events.png";fig.savefig(path,dpi=150);plt.close(fig);outputs.append(path);plots.append(("TAD changes: primary insulation window",path))
                changes=primary.loc[primary.event.isin(["growing","shrinking","shifted","stable","strength_change"]),"length_change_bp"].dropna()
                if len(changes):
                    fig,ax=plt.subplots(figsize=(6,4));ax.hist(changes/1000,bins=min(20,len(changes)))
                    ax.set(xlabel="Case minus reference domain length (kb)",ylabel="Matched domains");fig.tight_layout()
                    path=figures/"tad_lengths.png";fig.savefig(path,dpi=150);plt.close(fig);outputs.append(path);plots.append(("TAD size changes",path))
    if samples and (c["RUN_QC"] or enabled(c,"tads")):
        from .science import context
        from .matrix import bounded_raw
        for sample in samples[:6]:
            resolution=c["INSULATION_RESOLUTION"] if enabled(c,"tads") else c["PS_RESOLUTION"]
            clr,bins,mask,labels,views,weights=context(c,sample["sample_id"],resolution)
            view=views.iloc[0];region=(view.chrom,int(view.start),min(int(view.end),int(view.start)+500000))
            if tad_tables:
                combined=pd.concat(tad_tables,ignore_index=True)
                primary=combined.loc[combined.window_bp==c["INSULATION_PRIMARY_WINDOW_BP"]]
                interesting=primary.loc[~primary.event.isin(["stable","uncertain"])]
                if len(interesting) or len(primary):
                    event=(interesting if len(interesting) else primary).iloc[0]
                    view=views.loc[(views.chrom==event.chrom)&(views.name==event.region)].iloc[0]
                    start=max(int(view.start),(int(event.start)-50000)//resolution*resolution)
                    region=(event.chrom,start,min(int(view.end),start+500000))
            lo,hi=clr.extent(region);data=bounded_raw(clr,region,region,c)*weights[lo:hi,None]*weights[None,lo:hi]
            fig,ax=plt.subplots(figsize=(4,4));im=ax.imshow(np.log1p(data),cmap="magma",origin="lower",extent=[region[1]/1e6,region[2]/1e6,region[1]/1e6,region[2]/1e6])
            tad_path=root/f"03_features/tads/samples/{sample['sample_id']}.{c['INSULATION_PRIMARY_WINDOW_BP']}.tsv"
            if enabled(c,"tads") and tad_path.exists():
                domains=pd.read_csv(tad_path,sep="\t")
                boundaries=sorted(set(domains.loc[domains.chrom==region[0],"start"])|set(domains.loc[domains.chrom==region[0],"end"]))
                for boundary in boundaries:
                    if region[1]<=boundary<=region[2]:ax.axvline(boundary/1e6,color="cyan",ls="--",lw=.7);ax.axhline(boundary/1e6,color="cyan",ls="--",lw=.7)
                ax.set(xlabel="Position (Mb); cyan = called domain boundaries",ylabel="Position (Mb)")
            ax.set_title(f"{sample['sample_id']} {region[0]}:{region[1]}-{region[2]}");fig.colorbar(im,ax=ax,label="log(1 + balanced contacts)")
            fig.tight_layout();path=figures/("locus."+sample["sample_id"]+".png");fig.savefig(path,dpi=150);plt.close(fig);outputs.append(path);plots.append(("Representative locus",path))
    validation="uncalibrated" if c["ALLOW_UNCALIBRATED_DEVELOPMENT"] else "software_run; consult per-comparison calibration status"
    normalization="inherited_upstream_weights" if any(c[k] for k in ("RUN_QC","RUN_COMPARTMENTS","RUN_INSULATION","RUN_LOOPS","RUN_AGGREGATES")) else "raw_counts_only"
    coordinate_scopes=[read_json(path) for path in sorted((root/"00_metadata/masks").glob("*.scope.json"))]
    summary_path=directory/"run_summary.json"
    atomic_json(summary_path,{"project":c["PROJECT"],"run_id":c["RUN_ID"],"validation_status":validation,"stages":statuses,"normalization":normalization,"coordinate_scopes":coordinate_scopes,"production_readiness":"requires calibration, real-data and server pilot validation"})
    outputs.append(summary_path)
    methods=directory/"methods.md"
    methods.write_text("# Methods\n\nImported raw counts and existing upstream weights are read without copying/rebalancing. Expected tables use the common eligible-pixel denominator, including observable zeros. Pools and thinning are optional descriptive products.\n\nFormal loops use independent frozen single-pixel candidates, global eligible-cis exposures and edgeR robust quasi-likelihood with normalization factors one. Insulation uses the fixed primary grid with direct limma empirical-Bayes models; boundary/domain/compartment summaries remain descriptive. Biological FDR families are separate from caller enrichment significance.\n\n[Cooler](https://cooler.readthedocs.io/en/latest/api.html), [cooltools](https://cooltools.readthedocs.io/), [HiCRep](https://github.com/dejunlin/hicrep), [edgeR](https://bioconductor.org/packages/edgeR/), [limma](https://bioconductor.org/packages/limma/).\n",encoding="utf-8")
    outputs.append(methods)
    if coordinate_scopes:
        with methods.open("a",encoding="utf-8") as handle:
            handle.write("\nAnalysis chromosomes: "+", ".join(coordinate_scopes[0]["analysis_chromosomes"])+". Exact registered bin coordinates and the intersection of usable bins are shared across biological observations. Stored weights are not recalculated. Common coordinates do not equalize upstream balancing scopes. QC native contact totals can include excluded chromosomes; eligible exposures and distance decay use only the common registered views. X/Y are retained when requested and available, with sparse or masked measurements unavailable rather than zero. Interpret sex-chromosome comparisons with sample sex/copy number and coverage.\n")
    esc=html.escape
    sections=[f"<h1>microc2compare — {esc(c['PROJECT'])}</h1><p>Run {esc(c['RUN_ID'])}; {esc(validation)}.</p>",
              "<p>Existing raw counts were read directly. Enabled balanced analyses reuse validated stored weights. Descriptive results carry no biological P values. Caller significance and differential FDR have different meanings.</p>","<h2>Stage status</h2>",pd.DataFrame(statuses).to_html(index=False,escape=True),"<h2>Biological samples</h2>",pd.DataFrame(samples)[["sample_id","condition","assay","subject_id"]].to_html(index=False,escape=True),"<h2>Figures</h2>"]
    sections.insert(-1,"<h2>Analysis settings</h2><p>Requested analyses: "+esc(", ".join(c.get("_requested_analyses",["configured full pipeline"])))+".</p><details><summary>Configuration snapshot, paths and parameter values</summary>"+pd.DataFrame([{"setting":k,"value":str(v)} for k,v in sorted(c.items())]).to_html(index=False,escape=True)+"</details>")
    identities=read_json(root/"00_metadata/inputs.json")
    if coordinate_scopes:
        sections.insert(-1,"<h2>Chromosome scope and inherited weights</h2><p>Analysis chromosomes: "+esc(", ".join(coordinate_scopes[0]["analysis_chromosomes"]))+". Matching coordinates and common usable bins are verified. Unselected sequences are excluded from comparisons without modifying the input maps.</p><p>Stored weights retain their upstream balancing scope. Selecting common chromosomes does not make those scopes identical. Native QC totals describe the entire input map; eligible exposures and distance decay describe the common views. For X/Y, review sample sex, copy number and coverage. Missing or masked measurements are not zeros.</p><details><summary>Per-resolution and per-sample chromosome scope</summary><pre>"+esc(json.dumps(coordinate_scopes,indent=2))+"</pre></details>")
    sections.insert(-1,"<details><summary>Registered input files and SHA256 identities</summary>"+pd.DataFrame(list(identities.values())).to_html(index=False,escape=True)+"</details>")
    if enabled(c,"tads"):
        sections.insert(-1,"<h2>TAD changes</h2><p>Domains are intervals between strong insulation boundaries. Growing/shrinking refers to length in base pairs. Event groups are descriptive; coverage failures and replicate disagreement are uncertain, not biological disappearance. Merge/split members are counted once per event group. No event-level biological P values or FDR are available. Mean log2 insulation and boundary prominence are distinct from intra-domain contact enrichment.</p>")
        if tad_tables:
            events=pd.concat(tad_tables,ignore_index=True)
            sections.insert(-1,events.groupby(["comparison_id","window_bp","event"]).size().reset_index(name="event_groups").to_html(index=False,escape=True))
            sections.insert(-1,"<details><summary>First 50 events; download complete tables below</summary>"+events.head(50).to_html(index=False,escape=True)+"</details>")
        else: sections.insert(-1,"<p>No supported event groups were available. Empty tables remain valid products; consult coverage, unsupported domains and sensitivity results.</p>")
    if volcano_statuses:
        sections.insert(-1,"<h2>Differential loop volcano plots</h2><p>Each point is a retained tested loop. Positive log2 effects indicate greater relative contact allocation in the case condition; negative effects indicate less. Height shows biological model P values. Colours use BH FDR across the complete retained testing family, not loop-caller significance. The optional effect threshold changes highlighting only. No absolute genome-wide contact change is inferred. Preliminary results remain preliminary; UNCALIBRATED DEVELOPMENT results are not approved for formal conclusions.</p>")
        for entry in volcano_statuses:
            sections.insert(-1,"<p><strong>"+esc(entry["comparison_id"])+"</strong>: "+esc(entry["reason"])+"</p>")
        sections.insert(-1,pd.DataFrame(volcano_statuses).to_html(index=False,escape=True))
    sections.extend(fine_sections)
    for title,path in plots:sections.append(f'<figure><figcaption>{esc(title)}</figcaption><img loading="lazy" src="{esc(str(path.relative_to(directory)).replace(chr(92),chr(47)))}" width="650"></figure>')
    sections.append("<h2>Download results</h2><ul>")
    current={p.resolve() for p in outputs}
    current.update(p.resolve() for p in (root/"00_metadata").rglob("*") if p.is_file())
    for record_path in (root/"completion").glob("*/*.done.json"):
        record=read_json(record_path)
        if record["stage"]=="report" or not enabled(c,record["stage"]):continue
        if record["stage"]=="hicrep" and not c["RUN_HICREP"]:continue
        if complete(c,record["stage"],record["scope"]):
            current.update((root/entry["path"]).resolve() for entry in record["outputs"] if not entry.get("discarded_intermediate"))
    for path in sorted(root.rglob("*")):
        if path.resolve() not in current:continue
        if not path.is_file() or path.suffix not in (".tsv",".json",".bed",".bedpe",".bedgraph",".npz",".md",".png",".svg"):continue
        if any(x in path.parts for x in ("completion","logs",".run-lock")):continue
        relative=Path("..")/path.relative_to(root)
        purpose=output_purpose(path.relative_to(root))
        sections.append(f'<li><a href="{esc(str(relative).replace(chr(92),chr(47)))}">{esc(str(path.relative_to(root)))}</a> — {esc(purpose)}</li>')
    sections.append("</ul>")
    page=directory/"index.html"
    page.write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>microc2compare report</title><style>body{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:1rem}table{border-collapse:collapse}td,th{padding:.4rem;border:1px solid #ddd}img{max-width:100%}figure{margin:2rem 0}</style>'+"\n".join(sections)+"</html>",encoding="utf-8")
    outputs.append(page)
    return outputs


def output_purpose(path):
    name=str(path).replace("\\","/")
    if "08_reports/volcano/" in name: return "Differential loop plot points or availability/count/exclusion metadata; original model statistics are preserved"
    if "figures/volcano." in name: return "Differential loop volcano export; biological P values on height, model BH FDR for colours"
    if name.endswith(".bedgraph"): return "Continuous browser track; excluded measurements are omitted"
    if name.endswith(".bed"): return "Zero-based half-open browser intervals; fourth column is the feature/event ID"
    if "tad_events" in name: return "Descriptive TAD changes, coordinates, support and coverage; no biological FDR"
    if "tads/sensitivity" in name: return "TAD event sensitivity to windows and boundary thresholds"
    if "tads/conditions" in name: return "Condition domains; unsupported tables retain replicate disagreement"
    if "tads/samples" in name: return "Per-biological-replicate insulation-defined domains"
    if "01_qc/hicrep" in name: return "Raw-count contact-map concordance; NA means undefined"
    if "01_qc/ps" in name: return "Distance-specific raw/balanced contacts and eligible-pixel denominators"
    if "01_qc/contacts" in name: return "Contact totals, excluded bins and eligible exposure"
    if "03_features/compartments" in name: return "Compartment vectors, eigenvalues or orientation QC"
    if "03_features/insulation" in name: return "Fixed-grid insulation, boundary strength and coverage"
    if "03_features/loops" in name: return "Loop discovery/consensus; caller significance is not differential FDR"
    if "03_features/stripes" in name: return "Stripe geometry, native caller descriptors and common-segment coverage; descriptive support only"
    if "browser_tracks" in name: return "UCSC tracks and optional HiGlass exports; see inventory for format and conversion status"
    if "04_quantification" in name: return "Raw feature counts, exposures, exclusions or domain O/E measurements"
    if "05_aggregates" in name: return "Per-replicate observed/expected aggregates and contribution counts"
    if "06_comparisons" in name: return "Comparison inputs, parameters or effects; inspect inference/validation status"
    if "00_metadata" in name: return "Configuration, processing evidence, input identities or analysis masks"
    return "Analysis product or methods record; see the module/output guide"
