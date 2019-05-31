###
### R routines for the R package mixmeta (c)
#
mixmetaSim <-
function(y, S, Psi, random, data, nsim=1, seed=NULL, ...) {
#
################################################################################
#  SET OBJECT
#
  # ADDITIONAL ARGUMENTS
  addarg <- list(...)
#
  # PREPARE AND CHECK y
  if (missing(data)) data <- NULL
  y <- eval(substitute(y),data,parent.frame())
  if(!is.matrix(y)) y <- as.matrix(y)
  nm <- rownames(y)
  k <- ncol(y) ; n <- nrow(y)
  if(any(nay <- is.na(y))) stop("missing values not allowed in 'y'")
#
  # RESET random AND EXTRACT RANDOM TERMS
  random <- getRandom(random)
  ranform <- getFullFormula(as.formula("~1",parent.frame()),random)
  randata <- if(is.null(random)) data.frame(seq(n)) else 
    model.frame(ranform, data, na.action=na.pass)
  if(any(is.na(randata))) stop("missing values not allowed in random terms")
  if(nrow(randata)!=n) stop("'y' and random terms of different length")
#
  # DEFINE GROUPS AND RE-ORDER
  groups <- getGroups(random,randata)
  ord <- do.call(order,lapply(seq(ncol(groups)),function(i) groups[,i]))
  groups <- groups[ord,,drop=FALSE]
  y <- y[ord,,drop=FALSE]
  randata <- randata[ord,]
#
  # LIST OF DESIGN MATRICES FOR RANDOM PART (ONLY IF NEEDED)
  Z <- getZ(random,randata)
#
  # Psi
  if(!is.null(Psi)) Psi <- checkPD(Psi,force=FALSE,error=TRUE,label="Psi")
  Psi <- getList(Psi)
#
  # OTHER DIMENSIONS
  gp <- groups[,1]
  m <- length(unique(gp))
  q <- if(is.null(Psi)) 0 else if(is.null(Z)) 1 else
    if(!is.list(Z)) ncol(Z) else sapply(Z,ncol)
#
  # PRODUCE S AS A MATRIX OF VECTORIZED ENTRIES (IF NEEDED INPUT COVARIANCES)
  # CREATE ARTIFACTS FOR y AND control
  S <- eval(substitute(S),data,parent.frame())
  S <- getS(S,y,NULL,NULL,ord,Scor=addarg$Scor,checkPD=addarg$checkPD)
#
  # CHECKS
  if(n<2L) stop("at least two units must be generated")
  if(length(Psi)!=length(q) || any(sapply(seq_along(Psi), function(i)
    any(dim(Psi[[i]])!=k*q[[i]])))) stop("'random' and 'Psi' not consistent")
#
################################################################################
#  CREATE LISTS
#
  #  REPEATED MEASURES
  rep <- do.call(cbind,lapply(seq(ncol(groups)), function(i)
    tapply(groups[,i],gp, function(xi) length(unique(xi)))))
  #
  # TRANSFORM Z AND S
  predlist <- lapply(seq(m),function(i) c(t(y[gp%in%i,])))
  nalist <- lapply(predlist,is.na)
  Zlist <- getZlist(Z,nay,groups,m,k,q)
  Slist <- getSlist(S,nay,groups,m,k,addSlist=addarg$addSlist,
    checkPD=addarg$checkPD)
#
  # EIGEN TRANSFORMATION OF MARGINAL (CO)VARIANCE ERROR
  eigenlist <- lapply(getSigmalist(Zlist,nalist,Psi,Slist),eigen)
#
################################################################################
# SIMULATE
#
  # DEFINE THE SEED (FROM simulate.lm)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
#
  # SAMPLE
  sim <- lapply(seq(nsim),function(i) matrix(unlist(mapply(function(pred,eig)
    pred+eig$vec%*%diag(sqrt(eig$val),length(pred))%*%rnorm(length(pred)),
    predlist,eigenlist)),n,k,byrow=TRUE,dimnames=list(nm,NULL))[order(ord),])
#
  # RETURN
  dropList(sim)
}
