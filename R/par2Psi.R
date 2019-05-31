###
### R routines for the R package mixmeta (c)
#
par2Psi <-
function(par, bscov, k, q, fix) {
#
################################################################################
# FUNCTION TO GENERATE THE MATRICES OF THE RANDOM PART FROM THE PARAMETERS
#
  # IF A SINGLE VECTOR STORING PARAMETERS FOR MULTIPLE LEVELS
  par <- getList(getPar(par,bscov,k,q))
#
  # DEFINE THE MATRICES
  Psi <- lapply(seq_along(bscov), function(i) {
    d <- k*q[i]
    switch(bscov[i],
      # IF UNSTRUCTURED, CROSSPRODUCT FROM CHOLESKY
      unstr = {
        L <- diag(0,d)
        L[lower.tri(L,diag=TRUE)] <- par[[i]]
        tcrossprod(L)
      },
      # DIAGONAL: THE EXPONENTIAL OF THE PARAMETERS
      diag = diag(exp(par[[i]]),d),
      # IDENTITY: THE EXPONENTIAL OF THE PARAMETER
      id = diag(d)*exp(par[[i]]),
      # COMPOUND SYMMETRY
      cs = {
        R <- matrix((exp(par[[i]][2])-1L/(d-1L))/(exp(par[[i]][2L])+1L),d,d)
        R[row(R) == col(R)] <- 1L
        exp(par[[i]][1L]*2L) * R
      },
      # HETEROSCEDASTIC COMPOUND SYMMETRY
      hcs = {
        R <- matrix((exp(par[[i]][d+1L])-1L/(d-1L))/(exp(par[[i]][d+1])+1L),d,d)
        R[row(R) == col(R)] <- 1L
        D <- diag(sqrt(exp(par[[i]][seq(d)]*2L)),d)
        D%*%R%*%D
      },
      # AUTOREGRESSIVE OF FIRST ORDER
      ar1 = {
        cor <- plogis(par[[i]][2L])*2L-1L
        R <- cor^abs(outer(seq(d),seq(d),"-"))
        D <- diag(sqrt(exp(par[[i]][1L])),d)
        D%*%R%*%D
      },
      # HETEROGENEOUS AUTOREGRESSIVE OF FIRST ORDER
      har1 = {
        cor <- plogis(par[[i]][d+1])*2L-1L
        R <- cor^abs(outer(seq(d),seq(d),"-"))
        D <- diag(sqrt(exp(par[[i]][seq(d)])),d)
        D%*%R%*%D
      },
      # PROPORTIONAL
      prop = {
        if(is.list(fix)) fix <- fix[[names(bscov)[[i]]]]
        exp(par[[i]])*fix
      },
      # KNOWN CORRELATION
      cor = {
        if(is.list(fix)) fix <- fix[[names(bscov)[[i]]]]
        inputcov(sqrt(exp(par[[i]])),cov2cor(fix))
      },
      # FIXED
      fixed = {
        if(is.list(fix)) fix <- fix[[names(bscov)[[i]]]]
        fix
      }
    )
  })
#
  # NAMES
  names(Psi) <- names(bscov)
#
  # DROP THE LIST STRUCTURE IF ONLY ONE COMPONENT
  dropList(Psi)
}
