\name{ml.loglik.fn}
\alias{ml.loglik.fn}
\alias{reml.loglik.fn}
\alias{ml.loglik.gr}
\alias{reml.loglik.gr}


\title{ Likelihood Functions for mixmeta Models }

\description{
These functions compute the value of the log-likelihood and the related vectors of first partial derivatives for random-effects meta-analytical models, in terms of model parameters. They are meant to be used internally and not directly run by the users.
}

\usage{

ml.loglik.fn(par, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const,
  bscov, fix)
ml.loglik.gr(par, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const,
  bscov, fix)

reml.loglik.fn(par, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const,
  bscov, fix)
reml.loglik.gr(par, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const,
  bscov, fix)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{n} units aggregated within \eqn{m} (outer-level) groups, \eqn{k} outcomes, \eqn{p} fixed-effects predictors, and \eqn{q} random-effects predictors:

  \item{par }{ a vector representing the random-effects parameters defining the random-effects (co)variance matrix (or multiple matrices for multilevel models).}
  \item{Xlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{Zlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the random-effects part of the model. Each element of this list represents a list of matrices correponding to the (optionally multiple) grouping levels of random effects. In each matrix, rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of group-specific vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-group (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nalist }{ a \eqn{m}-dimensional list of group-specific logical vectors, identifying missing outcomes.}
  \item{rep }{ matrix with \eqn{m} rows where each column identifies the number of repetitions (number of groups) for each grouping level. The first column (outer level) is by definition a vector of 1's.}
  \item{k, q, nall }{ number of outcomes, number of random-effects predictors (including the intercept), total number of observations (excluding missing), respectively. While usually all are scalars, in the case of multilevel models \code{q} can be a numeric vector representing the number of predictors for each level.}
  \item{const }{ value of the constant to be included in the (restricted) likelihood, therefore not computed in the \code{\link[=ml.igls]{iterative algorithms}}.}
  \item{bscov }{ a character vector defining the structure of the (co)variance matrix for each level or random effects. See \code{\link{mixmeta}}.}
  \item{fix }{ a matrix (or optionally a list of matrices for multilevel models) defining the fixed components of the random-effects part of the model. See \code{\link{mixmeta.control}} for details.}
}

\details{
These functions are called internally by fitting functions, in particular \code{\link{ml.newton}} and \code{\link{reml.newton}}, to compute the (restricted) log-likelihood and its first partial derivatives in terms of random-effects parameters for meta-analytical models.

These functions actually specify the \emph{profiled} version of the (restricted) likelihood, expressed only in terms of random-effects parameters, while the estimate of the fixed-effects coefficients is computed at each iteration using a generalized least squares estimator, based on the current value of the between-study (co)variance matrix. At convergence, the value of this profiled version is identical to the full (restricted) likelihood. This approach is computationally efficient, as it reduces the number of parameters in the optimization routine, especially for meta-regression models.

The random-effects parameters in \code{par} depends on the chosen \code{\link[=mixmetaCovStruct]{structure(s)}} for the random-effects (co)variance matrix (or multiple matrices for multilevel models). The parameterization ensures positive-definiteness. A Cholesky decomposition is then performed on the marginal (co)variance matrix in order to re-express the problem as standard least square equations, an approach which speeds up the computation of matrix inverses and determinants. These equations are finally solved through a QR decomposition, which guarantees stability. More details are provided in the references below.

Some parameters of the fitting procedures are determined through \code{\link{mixmeta.control}}. Specifically, the user can obtain the Hessian matrix of the estimated parameters (appropriately transformed, see \code{\link{mixmetaCovStruct}}) in the optimization function by setting \code{hessian=TRUE}, and specific settings of the optimization process can be defined by the control list argument \code{optim}. These values are passed to the optimization function \code{\link{optim}}.
}

\value{
\code{ml.loglik.fn} and \code{reml.loglik.fn} return the value of the (restricted) log-likelihood for a given set of parameters in \code{par}. \code{ml.loglik.gr} and \code{reml.loglik.gr} return instead the related vector of first partial derivatives.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Lindstrom MJ and Bates DM (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. \emph{Journal of the American Statistical Association}. \bold{83}(404):1014--1022.

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
