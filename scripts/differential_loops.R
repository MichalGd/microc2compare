#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(edgeR); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4) stop("Usage: differential_loops.R counts.tsv exposures.tsv parameters.json output.tsv")
script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if(length(script_arg)!=1L) stop("Cannot locate the loop uncertainty helper")
source(file.path(dirname(normalizePath(sub("^--file=","",script_arg))),
                 "differential_loop_uncertainty.R"))
d <- read.delim(args[1], check.names=FALSE, na.strings="NA")
ex <- read.delim(args[2], check.names=FALSE)
p <- fromJSON(args[3], simplifyVector=TRUE)
ids <- p$design$sample_ids
x <- as.matrix(p$design$matrix)
contrast <- as.numeric(p$design$contrast)
if(qr(x)$rank!=ncol(x) || nrow(x)-ncol(x)<2) stop("Nonestimable design or residual df <2")
y <- as.matrix(d[,ids,drop=FALSE])
L <- ex$eligible_cis_exposure[match(ids,ex$sample_id)]
if(any(!is.finite(L)|L<=0)) stop("Invalid global eligible-cis exposures")
valid <- rowSums(!is.finite(y))==0
if(any(y[is.finite(y)]<0) || any(y[is.finite(y)]!=floor(y[is.finite(y)]))) stop("Response must contain nonnegative raw integer counts")
abundance <- rowMeans(sweep(y,2,L,"/"), na.rm=FALSE)*min(L)
keep <- valid & rowSums(y,na.rm=TRUE)>=p$min_total & abundance>=p$min_mean_at_min_exposure
if(sum(keep)<p$min_features) stop("Too few retained independent features for prespecified dispersion-learning gate")
dge <- DGEList(counts=y[keep,,drop=FALSE], lib.size=L, norm.factors=rep(1,length(L)))
dge <- estimateDisp(dge,x,robust=TRUE)
fit <- glmQLFit(dge,x,robust=TRUE)
test <- glmQLFTest(fit,contrast=contrast)
results <- as.data.frame(test$table)
results$FDR <- p.adjust(results$PValue,"BH")
results$feature_id <- d$feature_id[keep]
fields <- loop_model_fields(fit,nrow(dge))
results$dispersion <- fields$nb
results$working_NB_dispersion <- fields$working_nb
results$QL_dispersion <- fields$ql
results$QL_dispersion_field <- fields$ql_field
results$QL_df_total <- test$df.total
# Approximate moderated Wald uncertainty on log2 allocation effects. The QL
# F-test remains the inferential test; these intervals are explicitly approximate.
se <- loop_approximate_se(fit,x,contrast,fields)
results$SE_log2_approx_QL <- se
critical <- qt(0.975,test$df.total)
results$CI_lower_pointwise95_approx <- results$logFC-critical*se
results$CI_upper_pointwise95_approx <- results$logFC+critical*se
results$uncertainty_method <- "approximate_QL_information_Wald; not simultaneous"
for(j in seq_along(ids)) {
  results[[paste0("fitted_count_",ids[j])]] <- fit$fitted.values[,j]
  results[[paste0("Pearson_residual_",ids[j])]] <- (y[keep,j]-fit$fitted.values[,j])/sqrt(fit$fitted.values[,j]+fields$working_nb*fit$fitted.values[,j]^2)
}
out <- merge(d,results,by="feature_id",all.x=TRUE,sort=FALSE)
out$retained <- keep[match(out$feature_id,d$feature_id)]
out$filter_reason <- ifelse(out$retained,NA,ifelse(!valid[match(out$feature_id,d$feature_id)],"missing_measurement","prespecified_count_filter"))
out$inference_status <- p$inference
out$validation_status <- p$validation_status
write.table(out,args[4],sep="\t",quote=FALSE,row.names=FALSE,na="NA")
