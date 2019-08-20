###
### R routines for the R package mixmeta (c)
#
getContrXlev <-
function(formula, list)  {
#
################################################################################
# FUNCTION TO EXTRACT THE CONTRASTS/LEVELS RELATED TO A GIVEN FORMULA
#
  # IF NO CONSTRASTS, RETURN NULL
  if(is.null(list)) return(NULL)
#
  # RETURN CONSTRASTS RELATED TO TERMS IN THE FORMULA
  vars <- vapply(attr(terms(formula), "variables"), deparse, "")[-1L]
  ind <- names(list) %in% vars
  if(any(ind)) list[ind] else NULL
}
