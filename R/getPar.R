###
### R routines for the R package mixmeta (c)
#
getPar <-
function(par, bscov, k, q)  {
#
################################################################################
# FUNCTION TO EXTRACT PARAMETERS FOR EACH RANDOM LEVEL
#
  # DETERMINE THE NUMBER AND CHECK
  npar <- getNpar(bscov,k,q)
  if(sum(npar)!=length(unlist(par))) stop("'par' of wrong length")
#
  # IF PROVIDED AS A SINGLE VECTOR FOR MULTIPLE LEVELS
  if(!is.list(par) && length(bscov)>1L) {
    # DEFINE STARTING/END POINTS
    end <- cumsum(npar)
    start <- c(1,end[-length(end)]+1)
    # EXTRACT
    par <- lapply(seq_along(bscov), function(i)
      if(npar[i]>0) par[start[i]:end[i]] else NULL)
    names(par) <- names(bscov)
  }
#
  par
}
