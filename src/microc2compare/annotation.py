"""Both-anchor gene/cCRE annotations following the ATACseq2tracks conventions."""
from pathlib import Path
import bisect
import gzip
import re
import numpy as np
import pandas as pd
from .contracts import atomic_json, digest, file_hash, load_structured, require, resolve


def resources(c):
    registry=load_structured(c["ANNOTATION_REGISTRY"])
    require(registry.get("reference_id")==c["REFERENCE_ID"], "Annotation reference mismatch")
    for key in ("gtf","ccre"):
        if not registry.get(key):
            require(key!="gtf", "Annotation requires a GTF")
            continue
        path=resolve(registry[key],Path(c["ANNOTATION_REGISTRY"]).parent)
        require(registry.get(key+"_sha256")==file_hash(path), f"Annotation fingerprint mismatch: {key}")
        registry[key]=path
    require(registry.get("gtf_release"), "Record annotation GTF release/source")
    require(not registry.get("ccre") or registry.get("ccre_source"), "Record cCRE source and assembly-conversion provenance")
    return registry


def read_gtf(path, selected):
    genes,exons=[],[]
    opener=gzip.open if str(path).endswith(".gz") else open
    with opener(path,"rt",encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip(): continue
            f=line.rstrip().split("\t")
            require(len(f)==9, "GTF must have nine columns")
            if f[0] not in selected or f[2] not in ("gene","exon"): continue
            attrs=dict(re.findall(r'(\w+)\s+"([^\"]*)"',f[8]))
            gid=attrs.get("gene_id")
            require(gid and f[6] in ("+","-"), "GTF gene IDs and strands are required")
            row={"chrom":f[0],"start":int(f[3])-1,"end":int(f[4]),"gene_id":gid,"gene_name":attrs.get("gene_name",gid),"strand":f[6]}
            require(0<=row["start"]<row["end"], "Invalid GTF interval")
            (genes if f[2]=="gene" else exons).append(row)
    require(genes, "No GTF gene records overlap the selected chromosome scope")
    fields=["chrom","start","end","gene_id","gene_name","strand"]
    return pd.DataFrame(genes,columns=fields).drop_duplicates(),pd.DataFrame(exons,columns=fields).drop_duplicates()


def read_ccre(path, selected, registry):
    records=[];opener=gzip.open if str(path).endswith(".gz") else open
    with opener(path,"rt",encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#","track","browser")): continue
            f=line.rstrip().split("\t")
            require(len(f)>=4,"cCRE BED requires coordinates and an ID")
            if f[0] not in selected: continue
            classes=[]
            if "ccre_class_column" in registry:
                index=int(registry["ccre_class_column"])
                require(4<=index<len(f),"Configured cCRE class column must follow BED ID")
                classes=re.split(r"[;,]",f[index])
            else:
                for value in f[4:]:
                    classes.extend(x for x in re.split(r"[;,]",value) if x in ("PLS","pELS","dELS","CTCF-only","DNase-H3K4me3"))
            require(classes, "cCRE classes unavailable; declare the zero-based ccre_class_column")
            require(0<=int(f[1])<int(f[2]),'Invalid cCRE interval')
            records.append({"chrom":f[0],"start":int(f[1]),"end":int(f[2]),"ccre_id":f[3],"ccre_class":";".join(sorted(set(classes)))})
    return pd.DataFrame(records,columns=["chrom","start","end","ccre_id","ccre_class"])


def annotate_anchors(anchors, genes, exons, ccre, c, ccre_enabled):
    import bioframe
    promoters=genes.copy()
    tss=np.where(genes.strand.eq("+"),genes.start,genes.end-1)
    promoters["start"]=np.maximum(0,np.where(genes.strand.eq("+"),tss-c["ANNOTATION_PROMOTER_UPSTREAM"],tss-c["ANNOTATION_PROMOTER_DOWNSTREAM"]))
    promoters["end"]=np.where(genes.strand.eq("+"),tss+c["ANNOTATION_PROMOTER_DOWNSTREAM"]+1,tss+c["ANNOTATION_PROMOTER_UPSTREAM"]+1)
    def overlaps(frame):
        if frame.empty:return {}
        joined=bioframe.overlap(anchors,frame,how="inner",suffixes=("","_hit"))
        joined=joined.rename(columns={key+'_hit':key for key in frame if key not in ('chrom','start','end')})
        return {aid: group for aid,group in joined.groupby("anchor_id",sort=False)}
    prom,exon,body,reg=map(overlaps,(promoters,exons,genes,ccre))
    nearest={}
    for chrom,group in genes.assign(tss=tss).groupby("chrom",sort=False):
        rows=group.sort_values(["tss","gene_id"]).to_dict("records")
        nearest[chrom]=([r["tss"] for r in rows],rows)
    def values(hits,key):return ";".join(sorted(set(hits[key].dropna().astype(str)))) if hits is not None else ""
    rows=[]
    for a in anchors.itertuples():
        p,e,g,r=prom.get(a.anchor_id),exon.get(a.anchor_id),body.get(a.anchor_id),reg.get(a.anchor_id)
        context="promoter" if p is not None else "exon" if e is not None else "intron" if g is not None else "distal_intergenic"
        row={"anchor_id":a.anchor_id,"chrom":a.chrom,"start":a.start,"end":a.end,"gene_context":context,
             "promoter_overlap":p is not None,"exon_overlap":e is not None,"gene_body_overlap":g is not None,
             "promoter_gene_ids":values(p,"gene_id"),"promoter_gene_names":values(p,"gene_name"),
             "overlapping_gene_ids":values(g,"gene_id"),"nearest_gene_id":"","nearest_gene_name":"","nearest_tss_signed_distance_bp":np.nan,
             "nearest_gene_ids_all_ties":"","ccre_status":"no_overlap" if ccre_enabled else "disabled",
             "ccre_all_ids":values(r,"ccre_id"),"ccre_all_classes":values(r,"ccre_class"),"ccre_primary_id":"","ccre_primary_class":"",
             "enhancer_like":False if ccre_enabled else pd.NA}
        if a.chrom in nearest:
            positions,items=nearest[a.chrom]; midpoint=(a.start+a.end)/2
            at=bisect.bisect_left(positions,midpoint)
            probes=[i for i in (at-1,at) if 0<=i<len(items)]
            minimum=min(abs(items[i]["tss"]-midpoint) for i in probes)
            ties=[]
            for point in {items[i]["tss"] for i in probes if abs(items[i]["tss"]-midpoint)==minimum}:
                ties.extend(items[bisect.bisect_left(positions,point):bisect.bisect_right(positions,point)])
            gene=min(ties,key=lambda x:x["gene_id"])
            row.update(nearest_gene_id=gene["gene_id"],nearest_gene_name=gene["gene_name"],
                       nearest_tss_signed_distance_bp=(midpoint-gene["tss"])*(1 if gene["strand"]=="+" else -1),
                       nearest_gene_ids_all_ties=";".join(sorted({x["gene_id"] for x in ties})))
        if r is not None:
            hits=r.copy();hits["overlap"]=np.minimum(a.end,hits.end_hit)-np.maximum(a.start,hits.start_hit)
            priority={"PLS":0,"pELS":1,"dELS":2,"DNase-H3K4me3":3,"CTCF-only":4}
            hits["priority"]=hits.ccre_class.map(lambda x:min(priority.get(y,99) for y in x.split(";")))
            primary=hits.sort_values(["overlap","priority","ccre_id"],ascending=[False,True,True]).iloc[0]
            row.update(ccre_status="overlap",ccre_primary_id=primary.ccre_id,ccre_primary_class=primary.ccre_class,
                       enhancer_like=bool({"pELS","dELS"}&set(row["ccre_all_classes"].split(";"))))
        rows.append(row)
    return pd.DataFrame(rows)


def run(c,scope):
    from .science import save, external_features, validate_metadata
    from .matrix import reference
    root=Path(c["run_dir"]);directory=root/"03_features/annotation"
    registry=resources(c);selected=reference(c)[0]["effective_scope"]["chromosomes"]
    catalogue=root/"03_features/loops/catalogue.tsv"
    features=pd.read_csv(catalogue,sep="\t") if c["RUN_LOOPS"] else external_features(c,validate_metadata(c)["samples"])[0]
    anchors={};links=[]
    for f in features.itertuples():
        link={"feature_id":f.feature_id}
        for k in (1,2):
            chrom,start,end=getattr(f,f"chrom{k}"),getattr(f,f"start{k}"),getattr(f,f"end{k}")
            if chrom not in selected:continue
            aid="anchor_"+digest([c["REFERENCE_ID"],chrom,start,end])[:20]
            anchors[aid]={"anchor_id":aid,"chrom":chrom,"start":start,"end":end};link[f"anchor{k}_id"]=aid
        if "anchor1_id" in link and "anchor2_id" in link:links.append(link)
    anchor_table=pd.DataFrame(list(anchors.values()),columns=["anchor_id","chrom","start","end"])
    genes,exons=read_gtf(registry["gtf"],selected)
    ccre=read_ccre(registry["ccre"],selected,registry) if registry.get("ccre") else pd.DataFrame(columns=["chrom","start","end","ccre_id","ccre_class"])
    if len(anchor_table):
        annotated=annotate_anchors(anchor_table,genes,exons,ccre,c,bool(registry.get("ccre")))
    else:
        # Header-only annotation remains a valid result for an empty call set.
        annotated=pd.DataFrame(columns=["anchor_id","chrom","start","end","gene_context","promoter_overlap","enhancer_like"])
    pairs=pd.DataFrame(links,columns=["feature_id","anchor1_id","anchor2_id"])
    for k in (1,2):
        pairs=pairs.merge(annotated.rename(columns={x:f"anchor{k}_{x}" for x in annotated if x!="anchor_id"}),left_on=f"anchor{k}_id",right_on="anchor_id",how="left",validate="many_to_one").drop(columns="anchor_id")
    def category(row):
        a,b=bool(row.anchor1_promoter_overlap),bool(row.anchor2_promoter_overlap)
        ea,eb=row.anchor1_enhancer_like,row.anchor2_enhancer_like
        if a and b:return "promoter_promoter"
        if (a and pd.notna(eb) and eb) or (b and pd.notna(ea) and ea):return "promoter_enhancer_like"
        if pd.notna(ea) and pd.notna(eb) and ea and eb:return "enhancer_like_enhancer_like"
        return "other_or_unresolved"
    pairs["candidate_category"]=pairs.apply(category,axis=1) if len(pairs) else pd.Series(dtype=str)
    outputs=[save(annotated,directory/"anchors.tsv"),save(pairs,directory/"loops.tsv")]
    directory.mkdir(parents=True,exist_ok=True)
    bed=directory/'anchors.bed';anchor_table[['chrom','start','end','anchor_id']].to_csv(bed,sep='\t',index=False,header=False)
    bedpe=directory/'loops.bedpe';features[['chrom1','start1','end1','chrom2','start2','end2','feature_id']].to_csv(bedpe,sep='\t',index=False,header=False)
    if bed.stat().st_size==0:bed.write_text('# No eligible anchors\n')
    if bedpe.stat().st_size==0:bedpe.write_text('# No eligible candidates\n')
    outputs.extend([bed,bedpe])
    info=directory/"resources.json";atomic_json(info,{"registry":registry,"chromosomes":selected,"promoter_upstream":c["ANNOTATION_PROMOTER_UPSTREAM"],"promoter_downstream":c["ANNOTATION_PROMOTER_DOWNSTREAM"],"nearest_gene_definition":"nearest anchor midpoint to strand-aware gene TSS; all ties retained","coordinate_convention":"zero-based half-open; GTF converted from one-based inclusive","interpretation":"cCRE/nearest-gene annotations do not establish cell-type activity or causal regulation"});outputs.append(info)
    return outputs
