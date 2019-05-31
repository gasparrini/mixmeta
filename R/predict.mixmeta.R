###
### R routines for the R package mixmeta (c)
#
predict.mixmeta <-
function(object, newdata, se=FALSE, ci=FALSE, vcov=FALSE,
  interval=c("confidence","prediction"), ci.level=0.95, format,
  aggregate=c("stat","y"), na.action=na.pass, ...) {
#
################################################################################
# CHECK ARGUMENTS AND SET DEFAULTS
#
  interval <- match.arg(interval,c("confidence","prediction"))
  if(missing(format)) format <- ifelse(vcov&&object$dim$k>1,"list","matrix")
  format <- match.arg(format,c("matrix","list"))
  aggregate <- match.arg(aggregate,c("stat","y"))
#
################################################################################
# CREATE THE OBJECTS
#
  # DETERMINE TYPE OF PREDICTION
  new <- !missing(newdata) && !is.null(newdata)
  predint <- interval=="prediction"
#
  # MODEL FRAME (GROUPING VARIABLES ARE NOT CONSIDERED)
  formula <- object$formula[c(1L,3L)]
  mf <- if(new) {
    fixterms <- attr(terms(formula),"term.labels")
    xlev <- object$xlevels[names(object$xlevels%in%fixterms)]
    model.frame(formula,newdata,na.action=na.action,xlev=xlev)
  } else object$model
  if(!is.null(class <- attr(terms(object),"dataClasses")))
    .checkMFClasses(class,mf)
  X <- if(new) model.matrix(formula,mf,object$contrasts) else
    model.matrix(object)
  Z <- if(predint) getZ(object$random,mf,object$contrasts) else NULL
  offset <- model.offset(mf)
  if(!is.null(offset) && length(offset)!=NROW(mf)) stop("wrong offset")
  nay <- if(new) matrix(FALSE,nrow(X),object$dim$k) else
    is.na(as.matrix(object$residuals))
  na.action <- attr(mf,"na.action")
#
  # DIMENSIONS AND NAMES
  n <- nrow(X)
  k <- object$dim$k
  q <- object$dim$q
  nm <- rownames(X)
#
  # TRANSFORM X IN LIST
  Xlist <- lapply(seq(n),function(i) X[i,,drop=FALSE]%x%diag(object$dim$k))
#
  # TRANSFORM Z IN LIST, PLUS OTHER OBJECTS FOR RANDOM PART (ONLY IF NEEDED)
  # NB: Zlist DIFFERENT THAN USUAL, AS NO ACCOUNT FOR GROUPING
  # NB: OTHER OBJECTS ONLY DEFINED FOR CALLING getSigmalist LATER
  Zlist <- if(predint&&!is.null(Z))
    getZlist(Z,matrix(FALSE,n,k),sapply(rep(n,length(q)),seq),n,k,q) else NULL
  Psi <- if(predint) object$Psi else NULL
  rep <- if(predint) matrix(1,n,length(Z)) else NULL
  nalist <- lapply(Xlist,function(x) rep(FALSE,object$dim$k))
#
################################################################################
#
  # COMPUTE PREDICTION, ACCOUNTING FOR OFFSET
  fitlist <- lapply(seq(n),function(i) {
    fit <- Xlist[[i]]%*%object$coefficients
    if(!is.null(offset)) fit <- fit+offset[i]
    return(drop(fit))})
#
  # COMPUTE (CO)VARIANCE AND STANDARD ERRORS
  # ADD RANDOM PART IF PREDICTION INTERVALS (NB: TRICK USING Sigmalist)
  Vlist <- lapply(Xlist,function(X) X%*%tcrossprod(object$vcov,X))
  Vlist <- getSigmalist(Zlist,nalist,Psi,Vlist)
  stderrlist <- lapply(Vlist,function(x) sqrt(diag(x)))
#
  # COMPACT
  fit <- do.call(rbind,fitlist)
  V <- do.call(rbind, lapply(Vlist,vechMat))
  stderr <- do.call(rbind,stderrlist)
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
  if(format=="matrix" && object$dim$k==1) {
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
  } else if(format=="matrix" && any(se,ci) && !vcov && aggregate=="y") {
    rownames(fit) <- nm
    fit <- lapply(colnames(fit), function(j) cbind(fit=fit[,j],se=stderr[,j],
      ci.lb=cilb[,j],ci.ub=ciub[,j])[,c(TRUE,se,ci,ci)])
    names(fit) <- object$lab$k
  # WHEN LIST, OR vcov AND AGGREGATE ON OUTCOME
  } else if(format=="list" || (vcov && aggregate=="y")) {
    fit <- lapply(seq(nrow(fit)),function(i) {
      temp <- list(fit=fit[i,])
      if(se) temp$se <- stderr[i,]
      if(ci) temp[c("ci.lb","ci.ub")] <- list(cilb[i,],ciub[i,])
      if(vcov) {
        temp$vcov <- xpndMat(V[i,])
        dimnames(temp$vcov) <- list(object$lab$k,object$lab$k)
      }
      temp <- lapply(temp,function(x) if(any(is.na(x))) NA else x)
      return(if(length(temp)>1L) temp else temp[[1]])
    })
    names(fit) <- nm
  } else rownames(fit) <- nm
#
  # SIMPLIFY IF ONLY ONE PREDICTION
  if(is.matrix(fit)) fit <- drop(fit)
  if(is.list(fit) && length(fit)==1L) fit[[1]] else fit
}
