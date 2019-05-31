###
### R routines for the R package mixmeta (c)
#
getFullFormula <-
function(formula, random)  {
#
################################################################################
# FUNCTION TO ADD RANDOM TERMS (PREDICTORS AND GROUPING VARS) IN FORMULA
#
  # IF random IS NULL
  if(is.null(random)) return(formula)
#
  # EXTRACT THE TERMS IN FORMULAE FOR FIXED AND RANDOM TERMS
  random <- getList(random)
  fixterms <- attr(terms(formula),"term.labels")
  modrandom <- lapply(random, function(x)
    formula(paste("~",gsub("|","+",deparse(x[[2]]),fixed=TRUE))))
  ranterms <- unlist(lapply(modrandom, function(x) attr(terms(x),"term.labels")))
#
  # ADD (MISSING) RANDOM TERMS TO FORMULA
  add <- unique(ranterms[!ranterms%in%fixterms])
  formula <- formula(paste(c(deparse(formula,width.cutoff=499L),add),
    collapse=" + "), env=environment(formula))
#
  formula
}
