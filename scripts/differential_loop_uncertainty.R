# Read edgeR fit fields without changing its fit or inferential test.
# edgeR 4.4 uses s2.post; older fits use var.post. See the edgeR manual.
loop_model_fields <- function(fit, n) {
  feature_vector <- function(value, name, positive=FALSE) {
    if(!is.numeric(value) || !is.null(dim(value)) ||
       !(length(value) %in% c(1L,n)) || any(!is.finite(value)) ||
       any(if(positive) value<=0 else value<0)) {
      stop(paste("Invalid edgeR fit field:",name))
    }
    if(length(value)==1L) value <- rep(value,n)
    as.numeric(value)
  }
  field <- if(!is.null(fit$s2.post)) "s2.post" else "var.post"
  ql <- feature_vector(fit[[field]],field,positive=TRUE)
  nb <- feature_vector(fit$dispersion,"dispersion")
  scale <- fit$average.ql.dispersion
  if(is.null(scale)) scale <- 1
  if(!is.numeric(scale) || length(scale)!=1L || !is.null(dim(scale)) ||
     !is.finite(scale) || scale<=0) stop("Invalid edgeR average.ql.dispersion")
  # The edgeR v4 adjusted-deviance fit/test uses NB dispersion divided by
  # average.ql.dispersion. This is the dispersion for the information matrix.
  working_nb <- nb/scale
  if(any(!is.finite(working_nb))) stop("Invalid edgeR working NB dispersion")
  list(nb=nb, working_nb=working_nb, ql=ql, ql_field=field)
}

loop_approximate_se <- function(fit, design, contrast, fields) {
  mu <- fit$fitted.values
  if(!is.matrix(mu) || !is.numeric(mu) ||
     nrow(mu)!=length(fields$ql) || ncol(mu)!=nrow(design) ||
     any(!is.finite(mu)) || any(mu<0)) stop("Invalid edgeR fitted values")
  if(!is.matrix(design) || any(!is.finite(design)) ||
     length(contrast)!=ncol(design) || any(!is.finite(contrast)))
    stop("Invalid information-matrix design or contrast")
  vapply(seq_len(nrow(mu)),function(i) {
    information <- crossprod(design,design*(mu[i,]/(1+fields$working_nb[i]*mu[i,])))
    covariance <- tryCatch(solve(information),error=function(e)
      matrix(NA_real_,ncol(design),ncol(design)))
    variance <- as.numeric(t(contrast)%*%covariance%*%contrast)*fields$ql[i]
    sqrt(variance)/log(2)
  },numeric(1))
}
