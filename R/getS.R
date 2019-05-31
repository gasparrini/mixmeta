###
### R routines for the R package mixmeta (c)
#
getS <-
function(S, y, narm=NULL, subset=NULL, ord=NULL, Scor=NULL, checkPD=NULL) {
#
################################################################################
# TRANSFORM S IN A MATRIX OF VECTORIZED VCOV MATRICES
#
  # IF NULL, IT IS ASSUMED PASSED THROUGH control, ARRANGED LATER BY getSlist
  if(is.null(S)) return(NULL)
#
  # DIMENSIONS
  k <- dim(y)[2]
  n <- dim(y)[1]
#
  # MESSAGE
  mes <- "incorrect dimensions for 'S'"
#
  # IF DATAFRAME
  if(is.data.frame(S)) S <- as.matrix(S)
#
  # IF NUMERIC
  if(is.numeric(S)) {
    # IF JUST A VECTOR, TRANSFORM IN A MATRIX
    if(is.null(dim(S))) S <- as.matrix(S)
    # IF AN ARRAY, TRANSFORM IN A MATRIX
    if(length(dim(S))==3L) S <- t(apply(S,3,vechMat))
    # FINALLY, IF A MATRIX, CHECK DIMENSIONALITY
    if(!is.null(subset)) S <- S[subset,,drop=FALSE]
    if(!is.null(narm)) S <- S[-narm,,drop=FALSE]
    # INPUT CORRELATIONS (IF NEEDED) AND TRANSFORM S
    if(dim(S)[2]==k) S <- inputcov(sqrt(S),Scor)
    if(dim(S)[1]!=n || !dim(S)[2] %in% c(k,k*(k+1)/2)) stop(mes)
  }
#
  # IF A LIST
  if(is.list(S)) {
    S <- lapply(S,as.matrix)
    if(!is.null(subset)) S <- S[subset]
    if(!is.null(narm)) S <- S[-narm]
    if(length(S)!=n) stop(mes)
    if(any(sapply(S,dim)!=k)) stop(mes)
    S <- if(k==1L) as.matrix(sapply(S,vechMat)) else t(sapply(S,vechMat))
  }
#
  # ORDER
  if(!is.null(ord)) S <- if(is.null(dim(S))) S[ord] else S[ord,,drop=FALSE]
#
  # CHECK POSITIVE-DEFINITENESS (NOT BY DEFAULT)
  if(!is.null(checkPD) && checkPD) Slist <- checkPD(lapply(seq(n), function(i)
    xpndMat(S[i,])),error=TRUE,label="S")
  #
  # NAMES
  rownames(S) <- rownames(y)
  nk <- colnames(y)
  colnames(S) <- if(!is.null(nk)) if(dim(S)[2]==k) nk else
    vechMat(outer(nk,nk,paste,sep=".")) else NULL
#
  S
}
