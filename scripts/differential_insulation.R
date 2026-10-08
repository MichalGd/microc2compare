#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(limma);library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("Usage: differential_insulation.R scores.tsv parameters.json output.tsv")
d <- read.delim(args[1],check.names=FALSE,na.strings="NA")
p <- fromJSON(args[2],simplifyVector=TRUE)
y <- as.matrix(d[,p$design$sample_ids,drop=FALSE])
x <- as.matrix(p$design$matrix)
contrast <- matrix(as.numeric(p$design$contrast),ncol=1)
if(qr(x)$rank!=ncol(x) || nrow(x)-ncol(x)<2) stop("Nonestimable design or residual df <2")
keep <- rowSums(!is.finite(y))==0
if(sum(keep)<2) stop("Insufficient complete fixed-grid insulation observations")
fit <- lmFit(y[keep,,drop=FALSE],x)
fit <- contrasts.fit(fit,contrast)
fit <- eBayes(fit,robust=TRUE,trend=TRUE)
se <- fit$stdev.unscaled[,1]*sqrt(fit$s2.post)
effect <- fit$coefficients[,1]
critical <- qt(0.975,fit$df.total)
out <- d
out$case_minus_reference_log2_insulation <- NA_real_
out$SE_moderated <- NA_real_
out$CI_lower_pointwise95 <- NA_real_;out$CI_upper_pointwise95 <- NA_real_
out$t_moderated <- NA_real_;out$df_total <- NA_real_
out$PValue <- NA_real_;out$FDR <- NA_real_
out$case_minus_reference_log2_insulation[keep] <- effect
out$SE_moderated[keep] <- se
out$CI_lower_pointwise95[keep] <- effect-critical*se
out$CI_upper_pointwise95[keep] <- effect+critical*se
out$t_moderated[keep] <- fit$t[,1];out$df_total[keep] <- fit$df.total
out$PValue[keep] <- fit$p.value[,1]
out$FDR[keep] <- p.adjust(fit$p.value[,1],"BH")
out$filter_reason <- ifelse(keep,NA,"missing_or_ineligible_fixed_grid_score")
baseline <- lmFit(y[keep,,drop=FALSE],x)
fitted <- baseline$coefficients%*%t(x)
for(j in seq_along(p$design$sample_ids)) {
  name <- p$design$sample_ids[j]
  out[[paste0("fitted_score_",name)]] <- NA_real_
  out[[paste0("residual_",name)]] <- NA_real_
  out[[paste0("fitted_score_",name)]][keep] <- fitted[,j]
  out[[paste0("residual_",name)]][keep] <- y[keep,j]-fitted[,j]
}
out$inference_status <- p$inference;out$validation_status <- p$validation_status
write.table(out,args[3],sep="\t",quote=FALSE,row.names=FALSE,na="NA")
