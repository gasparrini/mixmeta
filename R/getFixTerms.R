###
### R routines for the R package mixmeta (c)
#
getFixTerms <-
function(formula, terms)  {
#
################################################################################
# FUNCTION TO REMOVE RANDOM PART FROM TERMS OBJECT
#
  # IDENTIFY RANDOM TERMS AND, IF ANY, EXCLUDE THEM
  ind <- which(!attr(terms,"term.labels")%in%attr(terms(formula),"term.labels"))
  if(length(ind)) terms <- terms[-ind]
#
  terms
}
