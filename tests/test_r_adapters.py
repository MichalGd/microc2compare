import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO=Path(__file__).resolve().parents[1]


def command(script,*args):
    prefix=os.environ.get('M2C_TEST_R_PREFIX')
    if prefix:
        conda=shutil.which('conda')
        if not conda:pytest.skip('Conda required for explicit R test prefix')
        return [conda,'run','--no-capture-output','-p',prefix,'Rscript','--vanilla',str(REPO/'scripts'/script),*map(str,args)]
    r=os.environ.get('M2C_TEST_RSCRIPT') or shutil.which('Rscript')
    if not r:pytest.skip('Rscript/edgeR/limma runtime unavailable')
    return [r,'--vanilla',str(REPO/'scripts'/script),*map(str,args)]


@pytest.mark.r
def test_real_r_effect_direction_offsets_and_na(tmp_path):
    rng=np.random.default_rng(97);ids=['C1','C2','C3','T1','T2','T3'];n=1200
    design={'sample_ids':ids,'matrix':[[1,0]]*3+[[1,1]]*3,'contrast':[0,1]}
    params={'design':design,'validation_status':'uncalibrated','inference':'formal','fdr':.05,'min_total':20,'min_mean_at_min_exposure':5,'min_features':1000}
    p=tmp_path/'params.json';p.write_text(json.dumps(params))
    means=np.full((n,6),80.);means[:60,3:]*=3
    counts=pd.DataFrame(rng.negative_binomial(20,20/(20+means)),columns=ids);counts.insert(0,'feature_id',['f'+str(i) for i in range(n)])
    path=tmp_path/'counts.tsv';counts.to_csv(path,sep='\t',index=False)
    ex=pd.DataFrame({'sample_id':ids,'eligible_cis_exposure':[1000000]*6});e=tmp_path/'ex.tsv';ex.to_csv(e,sep='\t',index=False)
    output=tmp_path/'loops.tsv';subprocess.run(command('differential_loops.R',path,e,p,output),check=True)
    result=pd.read_csv(output,sep='\t').set_index('feature_id')
    assert result.loc[['f'+str(i) for i in range(60)],'logFC'].median()>1
    assert result.loc[['f'+str(i) for i in range(60,n)],'logFC'].abs().median()<.2
    assert result.SE_log2_approx_QL.gt(0).all() and result.FDR.between(0,1).all()
    assert np.isfinite(result.SE_log2_approx_QL).all()
    assert (result.CI_upper_pointwise95_approx>=result.CI_lower_pointwise95_approx).all()
    # Independent reference: edgeR's unchanged fit/test and the analytical
    # two-group information variance. This catches field/scalar/scale errors.
    reference=tmp_path/'edgeR-reference.tsv'
    reference_script=tmp_path/'reference.R'
    reference_script.write_text('''
suppressPackageStartupMessages(library(edgeR))
a <- commandArgs(trailingOnly=TRUE)
d <- read.delim(a[1],check.names=FALSE)
ids <- c("C1","C2","C3","T1","T2","T3")
x <- cbind(1,c(0,0,0,1,1,1))
y <- DGEList(as.matrix(d[,ids]),lib.size=rep(1e6,6),norm.factors=rep(1,6))
y <- estimateDisp(y,x,robust=TRUE)
f <- glmQLFit(y,x,robust=TRUE)
t <- glmQLFTest(f,contrast=c(0,1))
q <- if(!is.null(f$s2.post)) f$s2.post else f$var.post
phi <- rep_len(f$dispersion,nrow(y))
if(!is.null(f$average.ql.dispersion)) phi <- phi/f$average.ql.dispersion
w <- f$fitted.values/(1+phi*f$fitted.values)
se <- sqrt(q*(1/rowSums(w[,1:3])+1/rowSums(w[,4:6])))/log(2)
out <- data.frame(feature_id=d$feature_id,t$table,
                  FDR=p.adjust(t$table$PValue,"BH"),reference_SE=se)
write.table(out,a[2],sep="\\t",row.names=FALSE,quote=FALSE)
''')
    # command() resolves the same R runtime; use the absolute fixture script.
    subprocess.run(command(str(reference_script),path,reference),check=True)
    ref=pd.read_csv(reference,sep='\t').set_index('feature_id').loc[result.index]
    for col in ('logFC','PValue','F','FDR'):
        np.testing.assert_allclose(result[col],ref[col],rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(result.SE_log2_approx_QL,ref.reference_SE,rtol=1e-10,atol=1e-12)
    ex['eligible_cis_exposure']*=7;ex.to_csv(e,sep='\t',index=False)
    scaled=tmp_path/'scaled.tsv';subprocess.run(command('differential_loops.R',path,e,p,scaled),check=True)
    compare=pd.read_csv(scaled,sep='\t').set_index('feature_id')
    np.testing.assert_allclose(result.logFC,compare.loc[result.index,'logFC'],atol=.005)
    scores=pd.DataFrame(rng.normal(0,.1,(n,6)),columns=ids);scores.loc[:59,ids[3:]]+=1
    scores.loc[60,'C1']=np.nan;scores.insert(0,'feature_id',['f'+str(i) for i in range(n)])
    s=tmp_path/'scores.tsv';scores.to_csv(s,sep='\t',index=False)
    out=tmp_path/'ins.tsv';subprocess.run(command('differential_insulation.R',s,p,out),check=True)
    ins=pd.read_csv(out,sep='\t').set_index('feature_id')
    assert ins.loc[['f'+str(i) for i in range(60)],'case_minus_reference_log2_insulation'].median()>.8
    assert pd.isna(ins.loc['f60','PValue']) and pd.isna(ins.loc['f60','FDR'])
    assert (ins.CI_upper_pointwise95.dropna()>=ins.CI_lower_pointwise95.dropna()).all()


@pytest.mark.r
def test_r_loop_field_shapes_and_information_variance():
    subprocess.run(command('..'+os.sep+'tests'+os.sep+'test_loop_uncertainty.R',
                           REPO/'scripts/differential_loop_uncertainty.R'),check=True)
