\name{vechMat}
\alias{vechMat}
\alias{xpndMat}


\title{ Vectorization and Expansion of Symmetric Matrices }

\description{
The function \code{vechMat} transforms a symmetric matrix in a vector containing its lower triangular elements, taken by column. The function \code{xpndMat} reverses this transformation.
}

\usage{
vechMat(mat, diag=TRUE)

xpndMat(vech)
}

\arguments{
  \item{mat }{ a square matrix.}
  \item{vech }{ a vector.}
  \item{diag }{ a logical switch indicating if the diagonal entries must be included.}
}

\value{
A vector for \code{vechMat}, a symmetric matrix for \code{xnpdMat}.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\seealso{ 
See functions \code{vech} and \code{xpnd} in package \pkg{MCMCpack}.
}

\examples{
# GENERATE A POSITIVE-DEFINITE MATRIX, VECTORIZE IT AND THEN RE-EXPAND
(M <- crossprod(matrix(rnorm(9),3)))
(v <- vechMat(M))
xpndMat(v)

# EXTRACT VECTORIZED S, EXPAND TO A LIST, AND RE-VECTORIZE
(S <- as.matrix(berkey98[5:7]))
(Slist <- lapply(seq(nrow(S)), function(i) xpndMat(S[i,])))
t(sapply(Slist,vechMat))
}

\keyword{manip}
