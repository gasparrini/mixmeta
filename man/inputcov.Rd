\name{inputcov}
\alias{inputcov}


\title{ Input (Co)Variance Matrices }

\description{
This function inputs (co)variance matrices of a set of outcomes given the corresponding standard deviation and correlation values.
}

\usage{
inputcov(sd, cor=NULL)
}

\arguments{
  \item{sd }{ a \eqn{m \times k}{m x k} matrix of standard deviations for \eqn{k} outcomes in \eqn{m} matrices, or a vector for \eqn{k} outcomes in a single matrix.}
  \item{cor }{ either a vector of length 1, \eqn{m} or \eqn{k(k-1)/2}, or alternatively a \eqn{k \times k}{k x k} or \eqn{m \times k(k-1)/2}{m x k(k-1)/2} matrix. See Details.}
}

\details{
Depending the number of outcomes \eqn{k} and matrices \eqn{m}, the argument \code{cor} is interpreted as:
\itemize{
  \item if a vector of length 1 (a scalar), the same correlation for all the \eqn{k} outcomes for all the \eqn{m} matrices;
  \item if a vector of length \eqn{m}, the same correlation for all the \eqn{k} outcomes for each of the \eqn{m} matrices;
  \item if a vector of length \eqn{k(k-1)/2}, the lower triangular elements (without diagonal, taken by column) of the correlation matrix of the \eqn{k} outcomes, the same for all the \eqn{m} matrices;
  \item if a \eqn{k \times k}{k x k} matrix, the correlation matrix for the single matrix (only when \eqn{m}=1);
  \item if a \eqn{m \times k(k-1)/2}{m x k(k-1)/2} matrix, each row represents the lower triangular elements (without diagonal, taken by column) of the correlation matrix of the \eqn{k} outcomes for each of the \eqn{m} matrices.
}
}

\value{
For a single matrix, the (co)variance matrix itself. For multiple matrices, a \eqn{m \times k(k+1)/2}{m x k(k+1)/2} matrix, where each row represents the vectorized entries of the lower triangle (with diagonal, taken by column) of the related (co)variance matrix (see \code{\link{vechMat}}).
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\note{
This function is called internally by \code{\link{mixmeta}} for multivariate models to input the correlation(s) when only the within-unit variances are provided through the argument \code{S}. In this case, the correlation values are set through the argument \code{Scor} in the control list (see \code{\link{mixmeta.control}}).
}

\seealso{
See \code{\link{xpndMat}}. See \code{\link{mixmeta.control}}.
}

\examples{
# SOME RANDOM SD FOR A SINGLE MATRIX, WITH CONSTANT CORRELATION
(M <- inputcov(runif(4, 0.1, 3), 0.7))
# CHECK CORRELATION
cov2cor(M)

# NOW WITH A MORE COMPLEX CORRELATION STRUCTURE
(M <- inputcov(runif(3, 0.1, 3), c(0.7,0.2,0.4)))
cov2cor(M)

# MULTIPLE MATRICES
(V <- matrix(runif(5*3, 0.1, 3), 5, 3,
  dimnames=list(1:5, paste("V", 1:3, sep=""))))
inputcov(V, 0.6)

# WITH REAL DATA WHEN CORRELATIONS AVAILABLE
hyp
(S <- inputcov(hyp[c("sbp_se","dbp_se")], cor=hyp$rho))
# CHECK FIRST STUDY
cov2cor(xpndMat(S[1,]))

# USED INTERNALLY IN mixmeta
p53
inputcov(sqrt(p53[c("V1","V2")]), 0.5)
model <- mixmeta(cbind(y1,y2), S=cbind(V1,V2), data=p53, control=list(Scor=0.5))
model$S
}

\keyword{manip}
