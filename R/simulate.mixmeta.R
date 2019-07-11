###
### R routines for the R package mixmeta (c)
#
simulate.mixmeta <-
function(object, nsim=1, seed=NULL, ...) {
#
################################################################################
#
  # RECOVER OBJECTS
  y <- as.matrix(object$fitted.values)
  S <- object$S
  Psi <- object$Psi
  random <- object$random
  data <- object$model
  addarg <- c(object$control,list(contrasts=object$contrasts))
#
  # FILL MISSING
  nay <- is.na(y)
  k <- ncol(y)
  if(any(nay)) {
    yS <- inputna(y,S)
    y <- yS[,seq(k),drop=FALSE]
    S <- yS[,-seq(k),drop=FALSE]
  }
#
################################################################################
#
  # SIMULATE
  sim <- mixmetaSim(y,S,Psi,random,data,nsim,seed,addarg)
#
  # INCLUDE PARTIAL MISSING
  if(nsim==1L) sim <- list(sim)
  na.action <- object$na.action
  if(any(nay)||!is.null(na.action)) sim <- lapply(sim, function(x) {
    x[nay] <- NA
    napredict(na.action,x)
  })
  if(nsim==1L) sim <- sim[[1]]
#
  sim
}
