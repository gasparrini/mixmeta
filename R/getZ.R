###
### R routines for the R package mixmeta (c)
#
getZ <-
function(random, data, contrasts)  {
#
################################################################################
# FUNCTION TO DEFINE THE DESIGN MATRICES FOR THE RANDOM PART
#
  # IF random IS NULL, JUST RETURN NULL
  if(is.null(random)) return(NULL)
#
  # OTHERWISE, GENERATE THE LIST
  if(!is.list(random)) random <- list(random)
  Z <- lapply(random, function(form) {
    # REMOVE THE GROUPING FACTOR FROM THE FORMULA
    form[[2]] <- form[[2]][[2]]
    # DERIVE THE MODEL MATRIX
    model.matrix(form,data,contrasts)
  })
#
  # RETURN A SINGLE MATRIX IF SINGLE LEVEL
  if(length(Z)==1L) Z[[1L]] else Z
}
