###
### R routines for the R package mixmeta (c)
#
getQlist <-
function(Zlist, nalist, rep, k, q) {
#
################################################################################
#
  # IF Zlist IS NULL, CREATE LIST OF IDENTITY MATRICES (REMOVING MISSING ROWS)
  if(is.null(Zlist)) Zlist <- lapply(nalist, function(na)
    list(list(diag(length(na))[!na,,drop=FALSE])))
#
  # DESIGN MATRIX MAPPING THE PARAMETERS TO BE ESTIMATED
  # GOLDSTEIN COMP STAT & DATA ANAL 1992, PAGE 65
  # FIRST LOOP IN Zexp
  Qlist <- lapply(seq(Zlist), function(i) {
    # CBIND THE BLOCK-DIAGONAL PARTS
    Zexp <- do.call(cbind,lapply(seq(length(q)),function(j)
      bdiagMat(Zlist[[i]][[j]])))
    # SECOND LOOP: LEVEL OF RANDOM EFFECTS (THEN LISTS PUT TOGETHER)
    do.call(c,lapply(seq_along(q),function(j) {
      # DEFINE THE ROWS/COLS OF THE MATRIX, IDENTIFYING TERMS
      # DEFINE ALSO THE STARTING POINT IN Zexp CORRESPONDING TO THIS LEVEL
      rows <- vechMat(row(diag(q[j]*k)))
      cols <- vechMat(col(diag(q[j]*k)))
      start <- c(0,cumsum(q*k*rep[i,]))[j]
      # THIRD LOOP: SEQUENCE OF TERMS
      lapply(seq(rows),function(t) {
        # FOURTH LOOP: SEQUENCE OF REPETITIONS IN THAT LEVEL (THEN SUM AND VEC)
        sumList(lapply(seq(rep[i,j]),function(r) {
          # DEFINE INDICES IN Z
          ind1 <- start+(r-1)*(q[j]*k)+rows[t]
          ind2 <- start+(r-1)*(q[j]*k)+cols[t]
          # IF VARIANCE OR COVARIANCE
          if(ind1==ind2) tcrossprod(Zexp[,ind1]) else
            tcrossprod(Zexp[,ind1],Zexp[,ind2]) +
            tcrossprod(Zexp[,ind2],Zexp[,ind1])
        }))
      })
    }))
  })
#
  Qlist
}
