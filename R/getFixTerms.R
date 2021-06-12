###
### R routines for the R package mixmeta (c)
#
getFixTerms <-
function(formula, terms, data=NULL)  {
#
################################################################################
# FUNCTION TO REMOVE RANDOM PART FROM TERMS OBJECT
#
  # RE-DEFINE TERMS FOR FIXED PART ONLY
  # NB: NOT USING drop.terms AS HAS ISSUES WITH RE-ORDERING AND ':' INTERACTIONS
  fixterms <- terms(formula, data=data)
#
  # IDENTIFY DIFFERENT SETS OF VARIABLES
  allvar <- vapply(attr(terms, "variables"), deparse, "")
  fixvar <- vapply(attr(fixterms, "variables"), deparse, "")
  ind <- allvar %in% fixvar
#
  # ADD ADDITIONAL ATTRIBUTES DEFINED BY model.frame, IF NEEDED
  if(!is.null(predvars <- attr(terms, "predvars")) && sum(ind)>1L) 
    attr(fixterms, "predvars") <- predvars[ind]
  if(!is.null(dataClasses <- attr(terms, "dataClasses")) && sum(ind)>1L) 
    attr(fixterms, "dataClasses") <- dataClasses[ind[-1]]
#
  fixterms
}
