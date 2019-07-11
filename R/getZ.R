###
### R routines for the R package mixmeta (c)
#
getZ <-
function(random, data, contrasts=NULL)  {
#
################################################################################
# FUNCTION TO DEFINE THE DESIGN MATRICES FOR THE RANDOM PART
#
  # IF random IS NULL, JUST RETURN NULL
  if(is.null(random)) return(NULL)
#
  # OTHERWISE, GENERATE THE LIST
  random <- getList(random)
  Z <- lapply(random, function(form) {
    # REMOVE THE GROUPING FACTOR FROM THE FORMULA
    form[[2]] <- form[[2]][[2]]
    # EXTRACT THE CONTRASTS
    contr <- getContrXlev(form, contrasts)
    # DERIVE THE MODEL MATRIX
    model.matrix(form, data, contr)
  })
#
  # RETURN A SINGLE MATRIX IF SINGLE LEVEL
  dropList(Z)
}
