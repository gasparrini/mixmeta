\name{inputna}
\alias{inputna}


\title{ Input Missing Values }

\description{
This function augment data by replacing missing values. It can be used internally in \code{\link{mixmeta}} through the \code{\link[=mixmeta.control]{control}} list.
}

\usage{
inputna(y, S, inputvar=10^4)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{n} units and \eqn{k} outcomes:

  \item{y }{ a \eqn{n}-dimensional vector (for univariate models) or \eqn{m \times k}{m x k} matrix (for multivariate models) of outcomes.}
  \item{S }{ series of within-unit variances (or (co)variance matrices for multivariate models) of the estimated outcome(s). For univariate models, this is usually a \eqn{n}-dimensional vector. For multivariate models, it can be provided as: a \eqn{m}-dimensional list of \eqn{k \times k}{k x k} matrices; a tri-dimensional \eqn{k \times k \times m}{k x k x n} array; a matrix or data frame with \eqn{n} rows and \eqn{k(k+1)/2} or \eqn{k} columns, depending on the availability of the within-unit correlations.}
  \item{inputvar }{ multiplier for inputting the missing variances in \code{S}.}
}

\details{
The function augments the data by replacing missing values in the outcomes and the associated (co)variances. Specifically, it replaces missing outcomes and missing covariances (if provided) with 0, and missing variances with the largest observed variance multiplied by \code{inputvar}. This value is expected to be very high, by default \eqn{10^4}{10^4}, so that the corresponding observation contributes only negligibly to the final estimate.
}

\value{
A matrix with the first \eqn{k} column corresponding to the augmented outcomes, and the remaining \eqn{k(k+1)/2} or \eqn{k} columns (depending on the availability of the within-study covariances) corresponding to vectorized entries of the lower triangle of the related (co)variance matrices.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. In submission.

Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Jackson D, Riley R, White IR (2011). Multivariate meta-analysis: Potential and promise. \emph{Statistics in Medicine}. \bold{30}(20);2481--2498.

White IR (2009). Multivariate random-effects meta-analysis. \emph{Stata Journal}. \bold{9}(1):40--56.

White IR (2011). Multivariate random-effects meta-regression: updates to mixmeta. \emph{Stata Journal}. \bold{11}(2):255-270.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\note{
Data augmentation used to be the approach to deal with missing values in the first implementation of \pkg{mixmeta}. The current algorithms directly account for missing.

Inputting missing values can be useful when two or more outcomes are never observed jointly, and the estimation is entirely based on indirect comparison. This method can be applied in network meta-analysis, also called indirect treatment comparison.

This approach can produce different results than standard methods, especially when the occurrence of missing is substantial. Preliminary analyses indicate that likelihood-based estimation methods do not seem to be affected, while non-iterative estimators such as method of moments and variance components are more sensitive. The user should be careful on the application of missing augmentation.
}

\seealso{
See \code{\link{inputcov}} for inputting (co)variance matrices.
}

\examples{
# INSPECT THE DATA
head(smoking)

# STANDARD APPROACH TO MISSING DATA
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
mod1 <- mixmeta(y, S)
summary(mod1)

# WITH DATA AUGMENTATION
augdata <- inputna(y, S)
y <- augdata[,1:3]
S <- augdata[,-c(1:3)]
mod2 <- mixmeta(y, S)
summary(mod2)
# NB: SAME PARAMETER ESTIMATES, BUT WRONG NYUMBER OF OBS

# USED INTERNALLY IN mixmeta
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
mod3 <- mixmeta(y, S, control=list(inputna=TRUE))
summary(mod3)
# NOW RIGHT NUMBER OF OBS
}

\keyword{manip}
