###
### R routines for the R package mixmeta (c)
#
getRandom <-
function(random, method=NULL, env=NULL)  {
#
################################################################################
# FUNCTION TO RESET THE FORMULA(E) DEFINING THE RANDOM PART OF THE MODEL
#
  # RETURN NULL IF MISSING/NULL, AND ONLY FOR ml/reml/ METHODS
  if(missing(random) || is.null(random)) return(NULL)
  if(!is.null(method)&&!method%in%c("ml","reml","model.frame"))
    stop("'random' only meaningful for methods 'ml' or 'reml'")
#
  # KEEP THE ENVIRONMENT, TRANSFORM IN LIST, AND STORE NAMES (ALSO EMPTY)
  random <- getList(random)
  names <- if(is.null(nm <- names(random))) rep("",length(random)) else nm
#
  # CHECK EACH FORMULA
  random <- mapply(function(form,nm) {
    # CHECK IT IS A ONE-SIDED FORMULA
    if(!inherits(form,"formula")||length(form)!=2)
      stop("'random' must be a (list of) one-sided formula(e). See help(mixmetaFormula)")
    # SPLIT THE FORMULA IF THERE IS A GROUPING SYMBOL '|'
    split <- strsplit(deparse(form,width.cutoff=499L),"|",fixed=TRUE)[[1]]
    # DEFINE THE GROUPING FACTORS (SPLIT IF '/' SYMBOL, USE NAME IF NONE)
    groups <- if(length(split)==1) nm else strsplit(split[2],"/",fixed=TRUE)[[1]]
    if(any(groups==""))
      stop("Undefined grouping factors in 'random'. See help(mixmetaFormula)")
    # REPEAT THE FORMULA FOR EACH GROUPING FACTOR, SETTING ENVIRONMENT
    form <- lapply(groups, function(x)
      as.formula(paste(split[1],x,sep="|"), env=env))
    # RENAME (REMOVING SPACES)
    names(form) <- gsub(" ","",groups,fixed=T)
    return(form)
  },random,names,SIMPLIFY=FALSE)
#
  # UNLIST
  names(random) <- NULL
  random <- unlist(random,recursive=FALSE)
  #
  # DROP THE LIST STRUCTURE IF ONLY ONE COMPONENT
  dropList(random)
}
