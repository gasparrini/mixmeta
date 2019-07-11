###
### R routines for the R package mixmeta (c) Antonio Gasparrini 2012-2014
#
mixmeta.fit <-
function(X, Z, y, S, groups, method, bscov, control) {
#
################################################################################
# DIMENSIONS, NAMES, AND OPTIONAL AUGMENTATION
#
  # SET DIMENSIONS AND MISSING
  y <- as.matrix(y)
  nay <- is.na(y)
  k <- ncol(y)
  n <- nrow(y)
  gp <- groups[,1]
  m <- length(unique(gp))
  p <- ncol(X)
  q <- if(method=="fixed") 0 else if(is.null(Z)) 1 else
    if(!is.list(Z)) ncol(Z) else sapply(Z,ncol)
  nall <- sum(!nay)
  # STORE NAMES
  nk <- colnames(y)
  if(k>1L && is.null(nk)) nk <- paste("y",seq(k),sep="")
  np <- colnames(X)
  nq <- if(method=="fixed" || is.null(Z)) NULL else
    if(length(q)==1L) colnames(Z) else lapply(Z,colnames)
#
  # MISSING REPLACEMENT THROUGH DATA AUGMENTATION
  if(control$inputna) {
    augdata <- inputna(y,S,inputvar=control$inputvar)
    y <- augdata[,seq(k),drop=FALSE]
    S <- augdata[,-seq(k),drop=FALSE]
    nay[nay] <- FALSE
  }
#
################################################################################
#  CREATE LISTS
#
  #  REPEATED MEASURES
  rep <- do.call(cbind,lapply(seq(ncol(groups)), function(i)
    tapply(groups[,i],gp, function(xi) length(unique(xi)))))
#
  # TRANSFORM nay AND y, VECTORIZING THEM BY ROW (MULTIPLE OUTCOMES)
  nalist <- lapply(seq(m),function(i) c(t(nay[gp%in%i,])))
  ylist <- lapply(seq(m),function(i) c(t(y[gp%in%i,]))[!nalist[[i]]])
#
  # TRANSFORM X (EXPAND BY MULTIPLE OUTCOMES)
  Xlist <- lapply(seq(m),function(i)
    (X[gp%in%i,,drop=FALSE]%x%diag(k))[!nalist[[i]],,drop=FALSE])
#
  # TRANSFORM Z (EXPAND BY k AND LIST BY q)
  Zlist <- getZlist(Z,nay,groups,m,k,q)
#
  # TRANSFORM S (OPTIONALLY FROM INFO PROVIDED IN control)
  Slist <- getSlist(S,nay,groups,m,k,control$addSlist,control$checkPD)
#
################################################################################
# FIT THE MODEL
#
  # SELECT THE ESTIMATION METHOD
  fun <- paste("mixmeta",method,sep=".")
  fit <- do.call(fun,list(Xlist=Xlist,Zlist=Zlist,ylist=ylist,Slist=Slist,
    nalist=nalist,rep=rep,k=k,m=m,p=p,q=q,nall=nall,bscov=bscov,control=control))
#
################################################################################
# SET OBJECTS
#
  # MESSAGE OF NON-CONVERGENCE
  if(!is.null(fit$converged)&&!fit$converged) {
    warning("convergence not reached after maximum number of iterations")
  }
#
  # INCLUDE COMPONENTS
  fit$dim <- list(k=k,n=n,m=m,p=p,q=q)
  fit$df <- list(nall=nall,nobs=nall-(method=="reml")*fit$rank,
    df=nall-fit$df.residual,fixed=fit$rank,random=ifelse(method=="fixed",0,
    nall-fit$rank-fit$df.residual))
  fit$lab <- list(k=nk,p=np,q=nq)
#
  # CREATE FITTED VALUES AND RESIDUALS
  # WITH LABELS AND DROPPED DIMENSIONS
  temp <- rep(NA,prod(dim(y)))
  temp[!t(nay)] <- unlist(fit$fitted.values)
  fit$fitted.values <-
    matrix(temp,nrow(y),byrow=TRUE,dimnames=list(rownames(y),nk))
  fit$residuals <- y - fit$fitted.values
#
  # NAMES OF FIXED TERMS
  nfix <- if(k==1L) np else if(p==1L) nk else
    paste(rep(nk,p),rep(np,each=k),sep=".")
  names(fit$coefficients) <- nfix
  dimnames(fit$vcov) <- list(nfix,nfix)
#
  # NAMES OF RANDOM TERMS
  if(method!="fixed") {
    if(length(q)==1L) {
      nran <- if(k==1L) nq else if(q==1L) nk else
        paste(rep(nk,q),rep(nq,each=k),sep=".")
      dimnames(fit$Psi) <- list(nran,nran)
    } else {
      for(j in seq(q)) {
        nran <- if(k==1L) nq[[j]] else if(q[[j]]==1L) nk else
          paste(rep(nk,q[[j]]),rep(nq[[j]],each=k),sep=".")
        dimnames(fit$Psi[[j]]) <- list(nran,nran)
      }
    }
  }
#
  fit
}
