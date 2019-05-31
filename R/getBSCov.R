###
### R routines for the R package mixmeta (c)
#
getBSCov <-
function(bscov, random, method)  {
#
################################################################################
# FUNCTION TO RESET THE INFO ON BETWEEN-STUDY (CO)VARIANCE STRUCTURES
#
  # CHECK
  bscov <- match.arg(bscov,c("unstr","diag","id","cs","hcs","ar1","har1","prop",
    "cor","fixed"), several.ok=TRUE)
  if(bscov!="unstr" && !method%in%c("ml","reml","model.frame"))
    stop("structured 'bscov' only available for methods 'ml' or 'reml'")
#
  # VALUES:
  # - IF random NULL OR NOT A LIST, JUST THE FIRST VALUE
  # - IF bscov WITHOUT NAMES, REPLICATE IT
  # - IF bscov WITH NAMES, ONLY REPLACE THESE TERMS
  nm <- names(bscov)
  if(is.null(random) || !is.list(random)) {
    val <- bscov[[1]]
  } else if(is.null(nm)) {
    val <- rep(bscov,length=length(random))
    names(val) <- names(random)
  } else {
    if(!all(nm%in%names(random)))
      stop("elements in 'bscov' not corresponding to 'random'")
    val <- sapply(random,function(x) "unstr")
    val[nm] <- bscov
  }
#
  unlist(val)
}
