###
### R routines for the R package mixmeta (c)
#
getGroups <-
function(random, data)  {
#
################################################################################
# FUNCTION TO DEFINE THE GROUPING STRUCTURE FOR THE RANDOM PART
#
  # IF random IS NULL, RETURN JUST SEQUENCE OF ROWS (ONE OBS PER STUDY)
  if(is.null(random)) return(matrix(seq(nrow(data))))
#
  # EXTRACT A LIST WITH GROUPING VARIABLES
  random <- getList(random)
  groups <- lapply(random, function(form) {
    form[[2]] <- form[[2]][[3]]
    model.frame(form,data)[[1]]
  })
#
  # DEFINE GROUPING THROUGH FACTORS, ACCOUNTING FOR INTERNAL NESTING
  # NB: LEVELS FROM 2 ON ALWAYS NESTED IN 1, AS THE LATTER DEFINES THE LISTS
  groups[[1]] <- as.factor(groups[[1]])
  if((len <- length(groups))>1L) for(i in 2:len)
    groups[[i]] <- factor(paste(groups[[i-1]],groups[[i]],sep="-"))
#
  # TRANFORM IN MATRIX
  groups <- do.call(cbind,lapply(groups,unclass))
#
  groups
}
