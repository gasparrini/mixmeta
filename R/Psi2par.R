###
### R routines for the R package mixmeta (c)
#
Psi2par <-
function(Psi, bscov, k, q, fix) {
#
################################################################################
# FUNCTION TO EXTRACT THE PARAMETERS DEFINING THE RANDOM PART
#
  # IF ONLY ONE MATRIX, CREATE THE LIST
  Psi <- getList(Psi)
#
  # EXTRACT THE PARAMETERS
  par <- lapply(seq_along(bscov), function(i) {
    d <- k*q[i]
    switch(bscov[i],
      # IF UNSTRUCTURED, LOWER TRIANGULAR OF THE CHOLESKY DECOMPOSITION
      unstr = vechMat(t(chol(Psi[[i]]))),
      # DIAGONAL: THE LOGARITHM OF THE DIAGONAL ELEMENTS
      diag = log(diag(Psi[[i]])),
      # IDENTITY: THE LOGARITHM OF THE MEAN OF THE DIAGONAL ELEMENTS
      id = log(mean(diag(Psi[[i]]))),
      # COMPOUND SYMMETRY (WITH POSITIVE-DEFINITENESS CONSTRAINT)
      cs = {
        if(d<2L) stop("bscov='cs' only meaningful with more than 1 random term")
        cor <- mean(cov2cor(Psi[[i]])[row(Psi[[i]])!=col(Psi[[i]])])
        if(cor <= -1L/(d-1L)) cor <- -1L/d
        c(log(mean(diag(Psi[[i]])))/2L, log((cor+1L/(d-1L))/(1L-cor)))
      },
      # HETEROGENEOUS COMPOUND SYMMETRY (WITH POSITIVE-DEFINITENESS CONSTRAINT)
      hcs = {
        if(d<2L) stop("bscov='hcs' only meaningful with more than 1 random term")
        cor <- mean(cov2cor(Psi[[i]])[row(Psi[[i]])!=col(Psi[[i]])])
        if(cor <= -1L/(d-1L)) cor <- -1L/d
        c(log(diag(Psi[[i]]))/2L, log((cor+1L/(d-1L))/(1L-cor)))
      },
      # AUTOREGRESSIVE OF FIRST ORDER
      ar1 = {
        if(d<2L) stop("bscov='ar1' only meaningful with more than 1 random term")
        cor <- mean(cov2cor(Psi[[i]])[row(Psi[[i]])-col(Psi[[i]])==1L])
        c(log(mean(diag(Psi[[i]]))), (qlogis((cor+1L)/2L)))
      },
      # HETEROGENEOUS AUTOREGRESSIVE OF FIRST ORDER
      har1 = {
        if(d<2L) stop("bscov='har1' only meaningful with more than 1 random term")
        cor <- mean(cov2cor(Psi[[i]])[row(Psi[[i]])-col(Psi[[i]])==1L])
        c(log(diag(Psi[[i]])), (qlogis((cor+1L)/2L)))
      },
      # PROPORTIONAL: THE LOGARITHM OF THE MEAN RATIO OF THE DIAGONAL ELEMENTS
      prop = {
        if(is.list(fix)) fix <- fix[[names(bscov)[[i]]]]
        log(mean(diag(fix/Psi[[i]])))
      },
      # KNOWN CORRELATION: SAME AS DIAGONAL
      cor = log(diag(Psi[[i]])),
      # FIXED
      fixed = NULL
    )
  })
#
  # NAMES
  names(par) <- names(bscov)
#
  # DROP THE LIST STRUCTURE IF ONLY ONE COMPONENT
  dropList(par)
}
