\name{ml.igls}
\alias{ml.igls}
\alias{reml.rigls}
\alias{igls.iter}
\alias{rigls.iter}


\title{ IGLS and RIGLS Iterative Algorithms for mixmeta Models }

\description{
These functions implements (restricted) iterative generalized least squares (IGLS and RIGLS) algorithms for (restricted) maximum likelihood estimators for random-effects meta-analytical models. They are meant to be used internally and not directly run by the users.
}

\usage{

ml.igls(Psi, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const, bscov,
  fix, control)

reml.rigls(Psi, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, nall, const, bscov,
  fix, control)

igls.iter(Psi, Qlist, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, bscov,
  fix, control)

rigls.iter(Psi, Qlist, Xlist, Zlist, ylist, Slist, nalist, rep, k, q, bscov,
  fix,  control)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{n} units aggregated within \eqn{m} (outer-level) groups, \eqn{k} outcomes, \eqn{p} fixed-effects predictors, and \eqn{q} random-effects predictors:

  \item{Psi }{ a matrix (or a list of matrices for multilevel models) representing the initial estimate of the random-effects (co)variance matrix.}
  \item{Xlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{Zlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the random-effects part of the model. Each element of this list represents a list of matrices correponding to the (optionally multiple) grouping levels of random effects. In each matrix, rows corresponding to missing outcomes have been excluded.}
  \item{Qlist }{ a \eqn{m}-dimensional list of group-specific design matrices mapping the random-effects parameters to be estimated in \code{Psi}. See references below for details.}
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
These functions are called internally by the fitting functions \code{\link{mixmeta.ml}} and \code{\link{mixmeta.reml}} to perform (R)IGLS optimization algorithms for estimating random-effects meta-analytical models.

These estimators are not sensitive to the choice of the starting values, and quickly converge to the vicinity of the (restricted) maximum likelihood. The starting values in \code{Psi} are therefore defined by default as a matrix (or matrices) with a diagonal form and 0.001 variances, or otherwise selected by the user in the \code{control} argument of \code{\link{mixmeta}} (see \code{\link{mixmeta.control}}).

The functions \code{ml.igls} and \code{reml.rigls} first produce a design matrix that maps the entries of \code{Psi}, and then call \code{iter.igls} and \code{iter.rigls}, respectively, to obtain updated results at each iteration following a (R)IGLS procedure described in Goldstein and colleagues (1992). Convergence is assessed as (lack of) changes in \code{Psi}. Positive semi-definiteness is forced by setting the negative eigenvalues of the estimated matrix to a value close to 0 at each iteration (see \code{\link[=mixmeta.control]{control}}).
}

\value{
The functions \code{ml.igls} and \code{reml.rigls} return an intermediate list object, with components corresponding to the estimated random-effects (co)variance matrix (or list of matrices), its parameters, the maximum (restricted) log-likelihood value, an indicator of convergence, and the number of iterations. These are then re-processed, with other components added later within other functions to finalize an object of class \code{"mixmeta"} (see \code{\link{mixmetaObject}}). The functions \code{iter.igls} and \code{iter.rigls} return an updated version of \code{Psi}.
}

\references{
Sera F, Armstrong B, Blangiardo M, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. 2019;38(29):5429-5444. [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

Goldstein H (1992). Efficient computational procedures for the estimation of parameters in multilevel models based on iterative generalized least squares. \emph{Computational Statistics & Data Analysis}. \bold{13}(1):63--71.

Goldstein H (1986). Multilevel mixed linear model analysis using iterative generalized least squares. \emph{Biometrika}. \bold{73}(1):43--56.

Goldstein H (1989). Restricted unbiased iterative generalized least-squares estimation. \emph{Biometrika}. \bold{76}(3):622--623.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
As stated earlier, these functions are called internally by \code{\link{mixmeta.ml}} and \code{\link{mixmeta.reml}}, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. They are however exported in the namespace and documented for completeness.
}

\seealso{
See \code{\link{mixmeta.fit}} and \code{\link{mixmeta.ml}} for additional info on the fitting procedures. See \code{\link{mixmeta.control}} to determine specific parameters of the fitting procedures. See \code{\link{mixmetaCovStruct}} for (co)variance structures.
See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
