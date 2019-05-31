###
### R routines for the R package mixmeta (c)
#
qtest.mixmeta <-
function(object, ...) {
#
################################################################################
# RUN THE RELATED FIXED-EFFECTS MODEL
#
  # EXTRACT MODEL FRAME
  mf <- model.frame(object)
#
  # DEFINE ORDER AND MAIN GROUPING FACTOR (FIRST GROUP, ORDERED)
  groups <- getGroups(object$random,mf)
  ord <- do.call(order,lapply(seq(ncol(groups)),function(i) groups[,i]))
  #
  # RE-ORDER
  groups <- groups[ord,,drop=FALSE]
  mf <- mf[ord,,drop=FALSE]
#
  # EXTRACT MATRICES
  int <- attr(object$terms,"intercept")==1L
  y <- as.matrix(model.response(mf,"numeric"))
  X <- model.matrix(object)[ord,,drop=FALSE]
  S <- if(!is.null(object$S)) as.matrix(object$S)[ord,,drop=FALSE] else NULL
  nay <- is.na(y)
#
  # DIMENSIONS
  m <- object$dim$m
  k <- object$dim$k
  p <- object$dim$p
  gp <- groups[,1]
#
  # TRANSFORM IN LISTS
  nalist <- lapply(seq(m),function(i) c(t(nay[gp%in%i,])))
  ylist <- lapply(seq(m),function(i) c(t(y[gp%in%i,]))[!nalist[[i]]])
  Xlist <- lapply(seq(m),function(i)
    (X[gp%in%i,,drop=FALSE]%x%diag(k))[!nalist[[i]],,drop=FALSE])
  Slist <- getSlist(S,nay,groups,m,k,object$control$addSlist,
    object$control$checkPD)
#
  # FIT GLS (WITHOUT RANDOM PART)
  gls <- glsfit(Xlist,ylist,Slist,onlycoef=FALSE)
#
################################################################################
# COMPUTE THE STATS
#
  # GLOBAL
  Q <- drop(crossprod(gls$invtUy-gls$invtUX%*%gls$coef))
  df <- with(object$df,nall-fixed)
#
  # IF MULTIVARIATE, ADD OUTCOME-SPECIFIC
  if(k>1L) {
    coef <- matrix(gls$coef,ncol=k,byrow=TRUE)
    indS <- diag(xpndMat(seq(k*(k+1)/2)))
    Q <- c(Q,sapply(seq(k), function(i) sum((y[,i]-X%*%coef[,i])^2/S[,indS[i]])))
    df <- c(df,colSums(!nay,na.rm=TRUE)-p)
  }
#
  pvalue <- sapply(seq(length(Q)),function(i) 1-pchisq(Q[i],df[i]))
  names(Q) <- names(df) <- names(pvalue) <-
    if(k>1L) c(".all",object$lab$k) else object$lab$k
#
  qstat <- list(Q=Q,df=df,pvalue=pvalue,residual=p-int>0L,k=k)
  class(qstat) <- "qtest.mixmeta"
#
  qstat
}
