###
### R routines for the R package mixmeta (c)
#
bdiagMat <-
function(x) {
#
################################################################################
#
  # CHECKS AND RETURNS
  if(is.matrix(x)) return(x)
  if(!all(sapply(x,is.matrix))) {
   warning("non-matrix components trasformed in matrices")
   x <- lapply(x,as.matrix)
  }
  if(length(x)==1L) return(x[[1]])
#
  # COMPUTE DIMENSIONS AND START/END VALUES
  dim <- t(sapply(x,dim))
  end <- apply(dim,2,cumsum)
  start <- apply(end,2,function(x) c(1,x[-length(x)]+1))
#
  # GENERATE INDICATOR MATRIX
  matind <- array(seq(prod(colSums(dim))),colSums(dim))
  ind <- unlist(lapply(seq(nrow(dim)),function(i)
    matind[start[i,1]:end[i,1],start[i,2]:end[i,2]]))
#
  # CREATE A 0'S MATRIX AND FILL
  mat <- matrix(0,sum(dim[,1]),sum(dim[,2]))
  mat[ind] <- unlist(x)
#
  mat
}
