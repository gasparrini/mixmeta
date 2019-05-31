###
### R routines for the R package mixmeta (c)
#
getSlist <-
function(S, nay, groups, m, k, addSlist=NULL, checkPD=NULL) {
#
################################################################################
# TRANSFORM S IN A LIST OF MATRICES
#
  # IF S IS PROVIDED
  if(!is.null(S)) {
#
    # CHECK THAT S AND addSlist ARE NOT BOTH PROVIDED
    if(!is.null(addSlist)) stop("'addSlist' only allowed without 'S'")
#
    # CREATE THE LIST
    Slist <- lapply(seq(m),function(i)
      bdiagMat(lapply(which(groups[,1]%in%i),function(j)
        xpndMat(S[j,])[!nay[j,],!nay[j,],drop=FALSE])))
    if(any(is.na(unlist(Slist))))
      stop("missing pattern in 'y' and S' is not consistent")
#
    # RETURN
    return(Slist)
  }
#
  # IF NOT INSTEAD, MUST BE PROVIDED
  # CHECKS (ASSUMED CONSISTENT WITH ORDER OF GROUPING AND DIMENSIONS)
  # NOT NULL, A LIST, NO MISSING, RIGHT LENGTH, RIGHT DIMENSIONS, POS-DEF
  if(is.null(addSlist))
    stop("within-unit errors must be provided either through 'S' or 'control'")
  if(!is.list(addSlist)||length(addSlist)!=m)
    stop(paste("'addSlist' not consistent with required format and",
      "grouping length. See help(mixmeta.control)"))
  if(any(is.na(unlist(addSlist)))) stop("no missing allowed in 'addSlist'")
  ind <- sapply(seq(m), function(i)
    dim(as.matrix(addSlist[[i]]))!=sum(!nay[groups[,1]%in%i,]))
  if(any(ind)) stop("wrong dimensions in 'addSlist'. See help(mixmetaControl)")
  # CHECK POSITIVE-DEFINITENESS (BY DEFAULT)
  if(is.null(checkPD) || checkPD)
    addSlist <- checkPD(addSlist,error=TRUE,label="Slist")
#
  # RETURN
  return(addSlist)
}
