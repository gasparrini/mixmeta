\name{bdiagMat}
\alias{bdiagMat}


\title{ Block-Diagonal Expansion of a List of Matrices }

\description{
The function \code{bdiagMat} builds a single matrix with block-diagonal from a list of matrices.
}

\usage{
bdiagMat(x)
}

\arguments{
  \item{x }{ a list of matrices, or a single matrix.}
}

\value{
A matrix with block-diagonal form if \code{x} is a list, or otherwise \code{x} itself if a matrix.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\seealso{ 
See functions \code{bldiag} in package \pkg{metafor}.
}

\examples{
# GENERATE A LIST OF MATRICES, AND CREATE THE BLOCK-DIAGONAL MATRIX
(matlist <- list(matrix(1:4,2), matrix(1:8,2)))
bdiagMat(matlist)
}

\keyword{manip}
