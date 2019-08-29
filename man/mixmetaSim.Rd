\name{mixmetaSim}
\alias{mixmetaSim}
\alias{simulate.mixmeta}


\title{ Simulating Responses for mixmeta Models }

\description{
These functions simulate sets of responses (either univariate or multivariate) for a group of units, in terms of their mean (expected) values and within and between-group (co)variances. These sets of outcomes can be used in meta-analytical models for simulation purposes.
}

\usage{
mixmetaSim(y, S, Psi, random, data, nsim=1, seed=NULL, ...)

\method{simulate}{mixmeta}(object, nsim=1, seed=NULL, \dots)
}

\arguments{
In order to simulate \eqn{k} outcomes for \eqn{n} units:

  \item{y }{ a \eqn{n}-dimensional vector (for simulating univariate responses) or \eqn{m \times k}{m x k} matrix (for simulating multivariate responses) of mean (expected) outcomes.}
  \item{S }{ series of within-unit variances (or (co)variance matrices for multivariate models) of the estimated outcome(s). For the list of accepted format, see the argument with the same name in \code{\link{mixmeta}}. Covariances or more complex error structures can be passed through additional arguments. See Details below.}
  \item{Psi }{ the random-effects (co)variance matrix (or a list of matrices for multilevel models) of the outcomes. Dimension must be consistent with the specification of the random-effects structure in \code{random}.}
  \item{random }{ a one-sided formula (or a list of formulae for multilevel models) offering a symbolic description of the linear predictor(s) and grouping structure for the random-effects part of the model. See the argument with the same name in \code{\link{mixmeta}}.}
  \item{data }{ an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame), optionally containing the variables in \code{y}, \code{S}, and \code{random}.}
  \item{nsim }{ number of simulation sets.}
  \item{seed }{ an object specifying if and how the random number generator should be initialized.}
  \item{object }{  an object of class \code{"mixmeta"}.}
  \item{\dots }{ further optional arguments.}
}

\details{
The set(s) of responses can be simulated either from a fitted model, using the method function \code{simulate} for objects of class \code{"mixmeta"}, or directly through the function \code{mixmetaSim}. In the former case, the fitted values from the model are used as mean (expected) outcomes, together with the within-unit and estimated random-effects (co)variance structure. In the latter option, this information need to be provided by the user in the correct dimensions and forms.

Additional arguments can be passed in \code{'...'}. Specifically, arguments \code{Scor} and \code{addSlist} can be added to input missing within-unit error covariances, or to specify more complex within-unit error structures, respectively. Another argument can be \code{checkPD} (logical), that checks the semi-positive definiteness of the matrices). See  \code{\link{mixmeta.control}} for details.

The functions simulate the responses for each study separately from a marginal multivariate normal distribution with mean equal to the expected values and (co)variance equal to the sum of the within-unit errors and random-effects components. The computation is identical to that implemented in the function \code{mvrnorm} of the package \pkg{MASS}, involving a eigen decomposition of the marginal (co)variance matrix. Numerical negative definiteness is checked, and positive semi-definiteness then forced by truncating the eigenvalues at a value close to zero (see \code{\link[=mixmeta.control]{control}}).
}

\value{
If \code{nsim=1}, a matrix or vector of simulated \eqn{k} outcomes for the \eqn{n} units. If more simulation sets are required (\code{nsim} higher than 1), a list of matrices or vectors.
}

\references{
Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
Studies with missing values in the fitted values or in the components of the within (co)variances are excluded by \code{simulate}. Missing values are instead not accepted in \code{metaSim}.
}

\seealso{
See \code{\link{simulate}} for the general method function. See \code{\link{inputcov}} for inputting correlations. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# RUN A MODEL
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98)

# SIMULATE A NEW SET OF OUTCOMES
simulate(model)

# SIMULATE FROM SCRATCH: 3 OUTCOMES, 8 STUDIES
(y <- matrix(0, 8, 3))
(S <- inputcov(matrix(runif(8*3, 0.1, 2), 8, 3, dimnames=list(NULL,
  c("V1","V2","V3"))), cor=c(0,0.5,0.7)))
(Psi <- inputcov(1:3, cor=0.3))
mixmetaSim(y, S, Psi)

# 2 SIMULATION SETS
mixmetaSim(y, S, Psi, nsim=2)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
