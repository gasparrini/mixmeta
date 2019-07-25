###
### R routines for the R package mixmeta (c)
#
getInitPsi <-
function(init, bscov, k, q, fix, checkPD=NULL)  {
#
################################################################################
# FUNCTION TO SET INITIAL VALUES
#
  # GENERATE DEFAULT VALUES
  initPsi <- lapply(q, function(qi) diag(0.001,k*qi))
#
  # IF PROVIDED, REPLACE
  if(!is.null(init)) {
    init <- getList(init)
    # CHECK IF NAMES MATCH (REQUIRED FOR MULTIPLE LEVELS)
    ind <- if(length(q)>1L) match(names(init),names(bscov)) else 1L
    if(length(ind)==0L || length(ind)>length(q) || length(init)>length(q) ||
        any(is.na(ind))) stop("'initPsi' does not match random components")
    # EXPAND IF VECTORIZED
    initPsi[ind] <- lapply(init,function(x) if(is.vector(x)) xpndMat(x) else x)
    # CHECK POSITIVE-DEFINITENESS (BY DEFAULT)
    if(is.null(checkPD) || checkPD)
      initPsi <- getList(checkPD(initPsi,force=FALSE,error=TRUE,label="initPsi"))
    # CHECK DIMENSIONS
    if(any(sapply(seq_along(initPsi),function(i) any(dim(initPsi[[ind[i]]])!=k*q[i]))))
      stop("wrong dimennsions in initPsi")
  }
#
  # IF FIXED, REPLACE (CONSISTENCY OF fix ALREADY CHECKED IN getPsifix)
  ind <- which(bscov%in%c("fixed"))
  if(length(ind)>0L) {
    fix <- getList(fix)
    ind2 <- if(length(fix)==1L) 1 else match(names(bscov)[ind],names(fix))
    initPsi[ind] <- fix[ind2]
  }
#
  # DROP THE LIST STRUCTURE IF ONLY ONE COMPONENT
  dropList(initPsi)
}
