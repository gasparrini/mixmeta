\name{ml.newton}
\alias{ml.newton}
\alias{reml.newton}


\title{ Quasi-Newton Iterative Algorithms for mixmeta Models }

\description{
These functions implement quasi-Newton iterative algorithms for (restricted) maximum likelihood estimators for random-effects meta-analytical models. They are meant to be used internally and not directly run by the users.
}

\usage{

ml.newton(Psi, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const,
  bscov, fix, control)

reml.newton(Psi, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const,
  bscov, fix, control)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{n} units aggregated within \eqn{m} (outer-level) groups, \eqn{k} outcomes, \eqn{p} fixed-effects predictors, and \eqn{q} random-effects predictors:

  \item{Psi }{ a matrix (or a list of matrices for multilevel models) representing the initial estimate of the random-effects (co)variance matrix.}
  \item{Xlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{Zlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the random-effects part of the model. Each element of this list represents a list of matrices correponding to the (optionally multiple) grouping levels of random effects. In each matrix, rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of group-specific vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-group (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nalist }{ a \eqn{m}-dimensional list of group-specific logical vectors, identifying missing outcomes.}
  \item{rep }{ matrix with \eqn{m} rows where each column identifies the number of repetitions (number of groups) for each grouping level. The first column (outer level) is by definition a vector of 1's.}
  \item{k, q, nall }{ number of outcomes, number of random-effects predictors (including the intercept), total number of observations (excluding missing), respectively. While usually all are scalars, in the case of multilevel models \code{q} can be a numeric vector representing the number of predictors for each level.}
  \item{const }{ value of the constant to be included in the (restricted) likelihood, therefore not computed in the iterative algorithms.}
  \item{bscov }{ a character vector defining the structure of the (co)variance matrix for each level or random effects. See \code{\link{mixmeta}}.}
  \item{fix }{ a matrix (or optionally a list of matrices for multilevel models) defining the fixed components of the random-effects part of the model. See \code{\link{mixmeta.control}} for details.}
  \item{control }{ list of parameters for controlling the fitting process, usually internally set to default values by \code{\link{mixmeta.control}}.}
}

\details{
These functions are called internally by the fitting functions \code{\link{mixmeta.ml}} and \code{\link{mixmeta.reml}} to perform quasi-Newton iterative optimization algorithms for estimating random-effects meta-analytical models.

Starting values for the iterations are defined by \code{Psi}, representing a random-effects (co)variance matrix (or a list of matrices for multilevel models). In the default hybrid procedure (see \code{\link{mixmeta.ml}}), these are provided using few iterations of a \code{\link[=ml.igls]{(R)IGLS}} algorithm. If a full quasi-Newton method is used, the starting values are instead defined by default as a matrix (or matrices) with a diagonal form and 0.001 variances, or otherwise selected by the user in the \code{control} argument of \code{\link{mixmeta}} (see \code{\link{mixmeta.control}}).

The functions first re-define \code{Psi} as a set of random-effects parameters, depending on the chosen \code{\link[=mixmetaCovStruct]{structure(s)}}, using parameterizations that ensure the positive-definiteness of the estimated matrix (or matrices). Then, the function \code{\link{optim}} with \code{method="BFGS"} is called internally to perform the quasi-Newton optimization, using specific \code{\link[=ml.loglik.fn]{likelihood functions}} that compute the value of the (restricted) likelihood and (optionally) the vector of its first partial derivatives. The latter are used only in the case of basic random-effects structures, or otherwise the derivatives are computed numerically.

Some parameters of the optimization procedures are determined through \code{\link{mixmeta.control}}. Specifically, the user can obtain the Hessian matrix of the estimated parameters (appropriately transformed, see \code{\link{mixmetaCovStruct}}) in the optimization function by setting \code{hessian=TRUE}, and specific settings of the optimization process can be defined by the control list argument \code{optim}. These values are passed to the optimization function \code{\link{optim}}.
}

\value{
These functions return an intermediate list object, with components corresponding to the estimated random-effects (co)variance matrix (or list of matrices), the maximum (restricted) log-likelihood value, an indicator of convergence, the number of iterations, and optionally the Hessian matrix. These are then re-processed, with other components added later within other functions to finalize an object of class \code{"mixmeta"}. See \code{\link{mixmetaObject}}.
}

\references{
Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

Lindstrom MJ and Bates DM (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. \emph{Journal of the American Statistical Association}. \bold{83}(404):1014--1022.

Harville DA (1977) Maximum likelihood approaches to variance component estimation and to related problems. \emph{Journal of
the American Statistical Association}. \bold{72}(358):320--338.

Pinheiro JC and Bates DM (2000). \emph{Mixed-Effects Models in S and S-PLUS}. New York, Springer Verlag.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
As stated earlier, these functions are called internally by \code{\link{mixmeta.ml}} and \code{\link{mixmeta.reml}}, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. They are however exported in the namespace and documented for completeness.
}

\seealso{
See \code{\link{mixmeta.fit}} and \code{\link{mixmeta.ml}} for additional info on the fitting procedures. See \code{\link{mixmeta.control}} to determine specific parameters of the fitting procedures. See \code{\link{mixmetaCovStruct}} for (co)variance structures. See \code{\link{chol}} and \code{\link{qr}} for info on the Cholesky and QR decomposition.
See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
