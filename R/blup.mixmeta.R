###
### R routines for the R package mixmeta (c)
#
blup.mixmeta <-
function(object, se=FALSE, pi=FALSE, vcov=FALSE, pi.level=0.95, type="outcome",
  level, format, aggregate="stat", ...) {
#
################################################################################
# CHECK ARGUMENTS AND SET DEFAULTS
#
  if(missing(format)) format <- ifelse(vcov&&object$dim$k>1,"list","matrix")
  type <- match.arg(type,c("outcome","deviation"))
  format <- match.arg(format,c("matrix","list"))
  aggregate <- match.arg(aggregate,c("stat","outcome"))
  if(missing(level)) level <- length(object$dim$q)
  if(object$method=="fixed") level <- 0L
  if(level<0L||level>length(object$dim$q))
    stop("'level' non compatible with random levels")
#
################################################################################
# CREATE DESIGN MATRIX X, y AND S
#
  # EXTRACT MODEL FRAME
  mf <- model.frame(object)
  na.action <- object$na.action
  nm <- rownames(mf)
#
  # DEFINE ORDER AND MAIN GROUPING FACTOR (FIRST GROUP, ORDERED)
  groups <- getGroups(object$random,mf)
  ord <- do.call(order,lapply(seq(ncol(groups)),function(i) groups[,i]))
#
  # RE-ORDER, KEEP ONLY USED LEVELS
  groups <- groups[ord,seq(max(1,level)),drop=FALSE]
  mf <- mf[ord,,drop=FALSE]
#
  # EXTRACT MATRICES (EXCLUDE OFFSET IF ANY)
  y <- as.matrix(model.response(mf,"numeric"))
  offset <- model.offset(mf)
  if(!is.null(offset)) y <- y - offset
  X <- model.matrix(object)[ord,,drop=FALSE]
  Z <- if(level>0L) getZ(object$random,mf,object$contrasts) else NULL
  if(!is.null(Z)&&is.list(Z)) Z <- if(level==1L) Z[[1L]] else Z[seq(level)]
#
  # S AND OFFSET
  S <- if(!is.null(object$S)) as.matrix(object$S)[ord,,drop=FALSE] else NULL
#
################################################################################
# RE-CREATE LISTS
#
  #  DEFINE DIMENSIONS AND REPEATED MEASURES
  m <- object$dim$m
  k <- object$dim$k
  q <- if(level==0L) 0 else object$dim$q[seq(level)]
  gp <- groups[,1]
  rep <- do.call(cbind,lapply(seq(ncol(groups)), function(i)
    tapply(groups[,i],gp, function(xi) length(unique(xi)))))
#
  # FILL MISSING
  nay <- is.na(y)
  if(any(nay)) {
    yS <- inputna(y,S)
    y <- yS[,seq(k),drop=FALSE]
    S <- yS[,-seq(k),drop=FALSE]
    nay[nay] <- FALSE
  }
#
  # TRANSFORM IN LISTS (NB: NO MISSING, nalist CREATED ONLY FOR Zlist)
  nalist <- lapply(seq(m),function(i) rep(FALSE,k))
  ylist <- lapply(seq(m),function(i) c(t(y[gp%in%i,])))
  Xlist <- lapply(seq(m),function(i) X[gp%in%i,,drop=FALSE]%x%diag(k))
  Zlist <- getZlist(Z,nay,groups,m,k,q)
  Slist <- getSlist(S,nay,groups,m,k,object$control$addSlist,
    object$control$checkPD)
#
  # PREDICTED VALUES (IGNORING OFFSET) AND RESIDUALS
  predlist <- lapply(seq(m),function(i) Xlist[[i]]%*%object$coefficients)
  reslist <- mapply(function(y,pred) y-pred,ylist,predlist,SIMPLIFY=FALSE)
#
  # COMPUTE RANDOM PARTS
  if(level>0L && !is.null(Psi <- object$Psi) && is.list(Psi))
    Psi <- Psi[seq(level)]
  ZPZlist <- if(level==0L) NULL else getZPZlist(Zlist,nalist,Psi)
  if(!is.null(ZPZlist)) {
    Ulist <- mapply(function(ZPZ,S) chol(ZPZ+S),ZPZlist,Slist,SIMPLIFY=FALSE)
    invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
  }
#
################################################################################
#
  # COMPUTE THE COMPONENTS (POINT ESTIMATES, (CO)VARIANCE AND STANDARD ERRORS)
  complist <- lapply(seq(m),function(i) {
    # FIXED PART (SET TO 0 IF DEVIATIONS)
    blup <- if(type=="deviation") array(0,dim(predlist[[i]])) else predlist[[i]]
    V <- if(type=="deviation") matrix(0,nrow(Xlist[[i]]),nrow(Xlist[[i]])) else
      Xlist[[i]]%*%tcrossprod(object$vcov,Xlist[[i]])
    # RANDOM PART
    if(!is.null(ZPZlist)) {
      ZPZinvSigma <- ZPZlist[[i]] %*% tcrossprod(invUlist[[i]])
      blup <- blup + ZPZinvSigma%*%reslist[[i]]
      V <- V + ZPZlist[[i]] - ZPZinvSigma%*%ZPZlist[[i]]
    }
    # BIND (ADD OFFSET IF ANY)
    stderr <- sqrt(diag(V))
    seqlist <- lapply(seq(length(blup)/k),function(i) c(i*k-k+1,i*k))
    blup <- rbindList(lapply(seqlist, function(x) blup[x[1]:x[2],]), k)
    V <- rbindList(lapply(seqlist, function(x) vechMat(V[x[1]:x[2],x[1]:x[2]])),
      k*(k+1)/2)
    stderr <- rbindList(lapply(seqlist, function(x) stderr[x[1]:x[2]]), k)
    # RETURN
    return(list(blup,V,stderr))
  })
#
  # EXTRACT, COMPACT IN MATRICES (ADD OFFSET IF ANY), RE-ORDER, NAMES
  blup <- rbindList(lapply(complist,"[[",1), k)
  if(!is.null(offset) && type=="outcome") blup <- blup+offset
  blup <- blup[order(ord),,drop=FALSE]
  V <- rbindList(lapply(complist,"[[",2), k*(k+1)/2)[order(ord),,drop=FALSE]
  stderr <- rbindList(lapply(complist,"[[",3), k)[order(ord),,drop=FALSE]
  colnames(blup) <- colnames(stderr) <- object$lab$k
#
  # PAD WITH MISSING OBSERVATIONS
  if(!is.null(na.action)) {
    blup <- napredict(na.action,blup)
    V <- napredict(na.action,V)
    stderr <- napredict(na.action,stderr)
    if(class(na.action)%in%c("exclude","pass")) {
      nm <- napredict(na.action,nm)
      nm[na.action] <- names(na.action)
    }
  }
#
  # ADD COMPONENTS
  zvalpi <- qnorm((1-pi.level)/2,lower.tail=FALSE)
  pilb <- blup-zvalpi*stderr
  piub <- blup+zvalpi*stderr
#
################################################################################
# AGGREGATE AND RETURN
#
  # IF ONLY 1 OUTCOME, FORCE THE AGGREGATION ON IT
  if(format=="matrix" && k==1) {
    rownames(blup) <- nm
    if(se) blup <- cbind(blup,stderr)
    if(pi) blup <- cbind(blup,pilb,piub)
    if(vcov) blup <- cbind(blup,V)
    colnames(blup) <- c("blup","se","pi.lb","pi.ub","vcov")[c(TRUE,se,pi,pi,vcov)]
    if(ncol(blup)==1L) blup <- drop(blup)
  } else if(format=="matrix" && any(se,pi,vcov) && aggregate=="stat") {
    blup <- list(blup=blup)
    if(se) blup$se <- stderr
    if(pi) blup[c("pi.lb","pi.ub")] <- list(pilb,piub)
    if(vcov) blup$vcov <- V
    for(i in seq(blup)) rownames(blup[[i]]) <- nm
  # WHEN AGGREGATE ON OUTCOME (BUT NO vcov)
  } else if(format=="matrix" && any(se,pi) && !vcov && aggregate=="outcome") {
    rownames(blup) <- nm
    blup <- lapply(colnames(blup), function(j) cbind(blup=blup[,j],se=stderr[,j],
      pi.lb=pilb[,j],pi.ub=piub[,j])[,c(TRUE,se,pi,pi)])
    names(blup) <- object$lab$k
  # WHEN LIST, OR vcov AND AGGREGATE ON OUTCOME
  } else if(format=="list" || (vcov && aggregate=="outcome")) {
    blup <- lapply(seq(nrow(blup)),function(i) {
      temp <- list(blup=blup[i,])
      if(se) temp$se <- stderr[i,]
      if(pi) temp[c("pi.lb","pi.ub")] <- list(pilb[i,],piub[i,])
      if(vcov) {
        temp$vcov <- xpndMat(V[i,])
        dimnames(temp$vcov) <- list(object$lab$k,object$lab$k)
      }
      temp <- lapply(temp,function(x) if(any(is.na(x))) NA else x)
      return(dropList(temp))
    })
    names(blup) <- nm
  } else rownames(blup) <- nm
#
  blup
}
