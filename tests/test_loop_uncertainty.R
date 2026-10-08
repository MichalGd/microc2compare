# Base-R unit regressions; the separate adapter test exercises real edgeR.
source(commandArgs(trailingOnly=TRUE)[1])
expect_error <- function(expression,pattern) {
  error <- tryCatch({force(expression); NULL},error=function(e) e)
  stopifnot(inherits(error,"error"),grepl(pattern,conditionMessage(error),fixed=TRUE))
}
n <- 2L
mu <- rbind(rep(40,6),rep(80,6))
design <- cbind(1,c(0,0,0,1,1,1))
fit <- list(s2.post=c(2,3),dispersion=.2,average.ql.dispersion=2,fitted.values=mu)
fields <- loop_model_fields(fit,n)
stopifnot(identical(fields$nb,c(.2,.2)),identical(fields$working_nb,c(.1,.1)),
          fields$ql_field=="s2.post")
se <- loop_approximate_se(fit,design,c(0,1),fields)
# For equal fitted counts within two independent groups of three,
# Var(case-reference) = QL dispersion * 2 / (3 * NB information weight).
expected <- sqrt(c(2,3)*2/(3*c(40,80)/(1+.1*c(40,80))))/log(2)
stopifnot(length(se)==n,isTRUE(all.equal(se,expected,tolerance=1e-12)))
legacy <- fit
legacy$var.post <- legacy$s2.post; legacy$s2.post <- NULL
legacy$average.ql.dispersion <- NULL; legacy$dispersion <- c(.1,.1)
old <- loop_model_fields(legacy,n)
stopifnot(old$ql_field=="var.post",
          isTRUE(all.equal(loop_approximate_se(legacy,design,c(0,1),old),expected)))
both <- fit; both$var.post <- c(999,999)
stopifnot(identical(loop_model_fields(both,n)$ql,c(2,3)))
scalar <- fit; scalar$s2.post <- 2
stopifnot(identical(loop_model_fields(scalar,n)$ql,c(2,2)))
zero_nb <- fit; zero_nb$dispersion <- 0
stopifnot(all(loop_model_fields(zero_nb,n)$working_nb==0))
for(value in list(NULL,numeric(0),c(1,2,3),c(NA,1),c(Inf,1),c(0,1),c(-1,1),matrix(1,2,1))) {
  bad <- fit; bad$s2.post <- value
  expect_error(loop_model_fields(bad,n),"Invalid edgeR fit field")
}
for(value in list(NULL,numeric(0),c(1,2,3),c(NA,1),c(Inf,1),c(-1,1),matrix(1,2,1))) {
  bad <- fit; bad$dispersion <- value
  expect_error(loop_model_fields(bad,n),"Invalid edgeR fit field")
}
for(value in list(c(1,2),NA_real_,Inf,0,-1,matrix(1,1,1))) {
  bad <- fit; bad$average.ql.dispersion <- value
  expect_error(loop_model_fields(bad,n),"Invalid edgeR average.ql.dispersion")
}
zero_mu <- fit; zero_mu$fitted.values[1,] <- 0
overflow <- fit; overflow$average.ql.dispersion <- 1e-320
expect_error(loop_model_fields(overflow,n),"Invalid edgeR working NB dispersion")
zero_se <- loop_approximate_se(zero_mu,design,c(0,1),fields)
stopifnot(is.na(zero_se[1]),is.finite(zero_se[2]))
bad_mu <- fit; bad_mu$fitted.values[1,1] <- NA
expect_error(loop_approximate_se(bad_mu,design,c(0,1),fields),"Invalid edgeR fitted values")
expect_error(loop_approximate_se(fit,design,c(1,2,3),fields),"Invalid information-matrix")
cat("Loop uncertainty field/shape/scaling and analytical-variance regressions passed\n")
