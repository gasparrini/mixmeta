###
### R routines for the R package mixmeta (c)
#
getPsifix <-
function(fix, bscov, k, q, checkPD=NULL)  {
#
################################################################################
# FUNCTION TO RESET FIXED PARTS OF (CO)VARIANCE STRUCTURE
#
  # INDEX OF LEVELS WITH FIXED PARTS, RETURN NULL IF NOT NEEDED
  ind <- which(bscov%in%c("prop","cor","fixed"))
  if(length(ind)==0L) return(NULL)
#
  # CHECK WITH INFO PROVIDED
  mess <- "'Psifix' does not match random components"
  if(is.null(fix)) stop(mess)
  fix <- getList(fix)
  if(length(fix)>length(ind)) stop(mess)
  if(length(ind)>1L) {
    nm <- match(names(bscov)[ind],names(fix))
    if(any(is.na(nm))) stop(mess)
    fix <- fix[nm]
  } else if(length(fix)>1L) stop(mess)
#
  # EXPAND
  fix <- lapply(fix,function(x) if(is.vector(x)) xpndMat(x) else x)
  # CHECK POSITIVE-DEFINITENESS (BY DEFAULT)
  if(is.null(checkPD) || checkPD)
    fix <- checkPD(fix,force=FALSE,error=TRUE,label="Psifix")
  # CHECK DIMENSIONS
  if(any(sapply(seq_along(fix),function(i) any(dim(fix[[i]])!=k*q[i]))))
    stop("wrong dimennsions in Psifix")
#
  # DROP THE LIST STRUCTURE IF ONLY ONE COMPONENT
  dropList(fix)
}
