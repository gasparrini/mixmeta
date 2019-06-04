###
### R routines for the R package mixmeta (c)
#
predict.mixmeta <-
function(object, newdata, se=FALSE, ci=FALSE, vcov=FALSE, ci.level=0.95, format,
  aggregate="stat", na.action=na.pass, ...) {
#
################################################################################
# CHECK ARGUMENTS AND SET DEFAULTS
#
  if(missing(format)) format <- ifelse(vcov&&object$dim$k>1,"list","matrix")
  format <- match.arg(format,c("matrix","list"))
  aggregate <- match.arg(aggregate,c("stat","outcome"))
#
################################################################################
# CREATE THE OBJECTS
#
  # DETERMINE TYPE OF PREDICTION
  new <- !missing(newdata) && !is.null(newdata)
#
  # MODEL MATRIX AND OFFSET
  X <- if(!new) model.matrix(object) else {
    tt <- getFixTerms(object$formula, object$terms)
    contr <- getContrXlev(tt, object$contrasts)
    xlev <- getContrXlev(tt, object$xlev)
    mf <- model.frame(tt, newdata, na.action=na.action, xlev=xlev)
    if(!is.null(cl <- attr(tt, "dataClasses"))) .checkMFClasses(cl, mf)
    X <- model.matrix(tt, mf, contr)
  }
  offset <- if(new) model.offset(mf) else object$offset
#
  # DIMENSIONS AND NAMES
  n <- nrow(X)
  k <- object$dim$k
  q <- object$dim$q
  nm <- rownames(X)
# 
  # INFO ON MISSING
  nay <- if(new) matrix(FALSE, nrow(X), k) else
    is.na(as.matrix(object$residuals))
  na.action <- if(new) attr(mf,"na.action") else object$na.action
#
  # TRANSFORM X IN LIST
  Xlist <- lapply(seq(n),function(i) X[i,,drop=FALSE]%x%diag(k))
  nalist <- lapply(Xlist,function(x) rep(FALSE,k))
#
################################################################################
#
  # COMPUTE PREDICTION, IGNORING OFFSET
  fitlist <- lapply(seq(n),function(i) drop(Xlist[[i]]%*%object$coefficients))
#
  # COMPUTE (CO)VARIANCE AND STANDARD ERRORS
  Vlist <- lapply(Xlist,function(X) X%*%tcrossprod(object$vcov,X))
  stderrlist <- lapply(Vlist,function(x) sqrt(diag(x)))
#
  # COMPACT INTO A MATRIX, AND DEAL WITH OFFSET
  fit <- rbindList(fitlist, k)
  if(!is.null(offset)) fit <- fit+offset
  V <- rbindList(lapply(Vlist, vechMat), k*(k+1)/2)
  stderr <- rbindList(stderrlist, k)
  colnames(fit) <- colnames(stderr) <- object$lab$k
#
  # PAD WITH MISSING
  if(!is.null(na.action)) {
    fit <- napredict(na.action,fit)
    V <- napredict(na.action,V)
    stderr <- napredict(na.action,stderr)
    if(class(na.action)%in%c("exclude","pass")) {
      nm <- napredict(na.action,nm)
      nm[na.action] <- names(na.action)
    }
  }
#
  # ADD COMPONENTS
  zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
  cilb <- fit-zvalci*stderr
  ciub <- fit+zvalci*stderr
#
################################################################################
# AGGREGATE AND RETURN
#
  # IF ONLY 1 OUTCOME, FORCE THE AGGREGATION ON IT
  if(format=="matrix" && k==1) {
    rownames(fit) <- nm
    if(se) fit <- cbind(fit,stderr)
    if(ci) fit <- cbind(fit,cilb,ciub)
    if(vcov) fit <- cbind(fit,V)
    colnames(fit) <- c("fit","se","ci.lb","ci.ub","vcov")[c(TRUE,se,ci,ci,vcov)]
    if(ncol(fit)==1L) fit <- drop(fit)
  } else if(format=="matrix" && any(se,ci,vcov) && aggregate=="stat") {
    fit <- list(fit=fit)
    if(se) fit$se <- stderr
    if(ci) fit[c("ci.lb","ci.ub")] <- list(cilb,ciub)
    if(vcov) fit$vcov <- V
    for(i in seq(fit)) rownames(fit[[i]]) <- nm
  # WHEN AGGREGATE ON OUTCOME (BUT NO vcov)
  } else if(format=="matrix" && any(se,ci) && !vcov && aggregate=="outcome") {
    rownames(fit) <- nm
    fit <- lapply(colnames(fit), function(j) cbind(fit=fit[,j],se=stderr[,j],
      ci.lb=cilb[,j],ci.ub=ciub[,j])[,c(TRUE,se,ci,ci)])
    names(fit) <- object$lab$k
  # WHEN LIST, OR vcov AND AGGREGATE ON OUTCOME
  } else if(format=="list" || (vcov && aggregate=="outcome")) {
    fit <- lapply(seq(nrow(fit)),function(i) {
      temp <- list(fit=fit[i,])
      if(se) temp$se <- stderr[i,]
      if(ci) temp[c("ci.lb","ci.ub")] <- list(cilb[i,],ciub[i,])
      if(vcov) {
        temp$vcov <- xpndMat(V[i,])
        dimnames(temp$vcov) <- list(object$lab$k,object$lab$k)
      }
      temp <- lapply(temp,function(x) if(any(is.na(x))) NA else x)
      return(dropList(temp))
    })
    names(fit) <- nm
  } else rownames(fit) <- nm
#
  # SIMPLIFY IF ONLY ONE PREDICTION
  if(is.matrix(fit)) fit <- drop(fit)
  dropList(fit)
}
