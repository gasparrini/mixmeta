###
### R routines for the R package mixmeta (c)
#
glsfit <-
function(Xlist, ylist, Sigmalist, onlycoef=TRUE)  {
#
################################################################################
# FUNCTION TO COMPUTE THE GLS ESTIMATE + OPTIONAL INTERMEDIATE PRODUCTS
#
  Ulist <- lapply(Sigmalist,chol)
  invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
  invtUXlist <- mapply(function(invU,X) crossprod(invU,X),
    invUlist,Xlist,SIMPLIFY=FALSE)
  invtUylist <- mapply(function(invU,y) crossprod(invU,y),
    invUlist,ylist,SIMPLIFY=FALSE)
  invtUX <- do.call("rbind",invtUXlist)
  invtUy <- do.call("rbind",invtUylist)
  coef <- as.numeric(qr.solve(invtUX,invtUy))
#
  if(onlycoef) return(coef)
#
  list(coef=coef,Ulist=Ulist,invUlist=invUlist,invtUXlist=invtUXlist,
    invtUX=invtUX,invtUy=invtUy)
}
