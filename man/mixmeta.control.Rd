\name{mixmeta.control}
\alias{mixmeta.control}


\title{ Ancillary Parameters for Controlling the Fit in mixmeta Models }

\description{
This internal function sets the parameter options used for fitting meta-analytical models, commonly to pre-specified default values. It is usually internally called by \code{\link{mixmeta}}.
}

\usage{
mixmeta.control(optim=list(), showiter=FALSE, maxiter=100, initPsi=NULL, Psifix=NULL,
  Scor=NULL, addSlist=NULL, inputna=FALSE, inputvar=10^4, loglik.iter="hybrid",
  igls.inititer=10, hessian=FALSE, vc.adj=TRUE, reltol=sqrt(.Machine$double.eps),
  checkPD=NULL, set.negeigen=sqrt(.Machine$double.eps))
}

\arguments{

  \item{optim }{ list of parameters passed to the \code{control} argument of the function \code{\link{optim}}, which performs the \code{\link[=ml.newton]{quasi-Newton optimization}} in likelihood-based random-effects models. See \code{\link{optim}} for the list of arguments. See Details for additional info.}
  \item{showiter }{ logical. If \code{TRUE}, the progress of iterative optimization is shown.}
  \item{maxiter }{ positive interger value. Maximum number of iterations in methods involving optimization procedures.}
  \item{initPsi }{ either a matrix or a vector of its lower triangular elements (with diagonal, taken by column), or optionally a named list with one or more of such objects. Used as starting values of random-effects parameters in likelihood-based optimization routines. See Details.}
  \item{Psifix }{ either a matrix or a vector of its lower triangular elements (with diagonal, taken by column), or optionally a named list with one or more of such objects. Used to define fixed parts of the random-effects \code{\link[=mixmetaCovStruct]{(co)variance structures}}. See Details.}
  \item{Scor }{ either a scalar, vector or matrix representing the within-unit correlation(s) to be inputted when the covariances are not provided in multivariate models, and ignored if they are. See \code{\link{inputcov}}.}
  \item{addSlist }{ a list of \eqn{m} matrices for the (outer-level) groups of units defining the (known) error (co)variance structure, when this cannot be passed through the argument \code{S} of \code{\link{mixmeta}}. See Details.}
  \item{inputna }{ logical. If missing values must be internally inputted. To be used with caution. See \code{\link{inputna}}.}
  \item{inputvar }{ multiplier for inputting the missing variances, to be passed as an argument to \code{\link{inputna}}.}
  \item{loglik.iter }{ iterative scheme used in in likelihood-based optimization routines. Options are \code{"hybrid"}, \code{"newton"}, and \code{"igls"} or \code{"RIGLS"}. See \code{\link{mixmeta.ml}}.}
  \item{igls.inititer }{ number of iterations of the (restricted) iterative generalized least square algorithm when used in the initial phase of hybrid optimization procedure of likelihood-based estimators. See \code{\link{mixmeta.ml}}.}
  \item{hessian }{ logical. If \code{TRUE}, the Hessian matrix of the parameters estimated in the optimization process is computed and returned. Only applicable to likelihood-based estimation methods. For details, see the info provided in the help pages of the \code{\link[=ml.loglik.fn]{optimizations algorithms}} and \code{\link[=mixmetaCovStruct]{(co)variance structures}}.}
  \item{vc.adj }{ logical. If \code{TRUE}, an adjustement to the way the marginal variance part is computed in the (co)variance components estimator is applied in the variance components estimator. See \code{\link{mixmeta.vc}}.}
  \item{reltol }{ relative convergence tolerance in methods involving optimization procedures. The algorithm stops if it is unable to reduce the value by a factor of \code{reltol * (abs(val) + reltol)} at a step.}
  \item{checkPD }{ logical. Determines if the semi-positiveness of within-unit error or random-effects (co)variance matrices must be checked.}
  \item{set.negeigen }{ positive value. Value to which negative eigenvalues are to be set in estimators where such method is used to force semi-positive definiteness of the estimated between-study (co)variance matrix.}
}

\details{
This function has default values for most of the arguments, some of them set internally. Non-default values are passed through the control argument of \code{\link{mixmeta}}. Many arguments refer to specific fitting procedures. See the help page of the related estimator for details.

The function automatically sets non-default values for some control arguments for \code{\link{optim}}, unless explicitly set in the list passed to it. Specifically, the function selects \code{fnscale=-1}, \code{maxit=maxiter} and \code{reltol=reltol}, where the latter two are specified by other arguments of this function.

The arguments \code{initPsi} and \code{Psifix} are used to provide information for estimation procedures of the random-effects parameters in likelihood-based methods. Specifically, the former is used to choose non-default starting values (see \code{\link{mixmeta.ml}}), and the latter for defining the fixed (known) part of specific \code{\link[=mixmetaCovStruct]{(co)variance structures}}. In multilevel models, these arguments must be lists with named components referring to one or more levels of grouping defined by the argument \code{random} of \code{\link{mixmeta}}.

The argument \code{addSlist} can be used to define more complex (known) error structures of the outcome(s) that are usually provided through the argument \code{S} of \code{\link{mixmeta}} as within-unit variances (or (co)variance matrices for multivariate models). This can be useful when these error structures spans multiple units (rows), and the between-unit correlation cannot be defined through \code{S}, for instance in dose-response meta-analysis (see examples in \code{\link{mixmeta}}). Note that this information is passed internally after the data have be re-ordered following the grouping defined by \code{random} in \code{\link{mixmeta}}, and this should be consistent in \code{addSlist}. It is required to re-order the data before fitting the model, so to ensure consistency.
}

\value{
A list with components named as the arguments.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
The function is expected to be extended and/or modified at every release of the package \pkg{mixmeta}. It is strongly suggested to check the arguments of this function at every release.
}

\seealso{
See \code{\link{mixmeta}}. See also \code{\link{glm.control}}. See the help pages of the related fitting functions for details on each parameter. See \code{\link{mixmeta-package}} for an overview of this modelling framework.
}

\examples{
# PRINT THE ITERATIONS (SEE ?optim) AND CHANGE THE DEFAULT FOR STARTING VALUES
mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98,
  control=list(showiter=TRUE, igls.inititer=20))

# INPUT THE CORRELATION
mixmeta(cbind(y1,y2), S=cbind(V1,V2), data=p53, control=list(Scor=0.5))

# FIX (PARTS OF) THE RANDOM-EFFECTS (CO)VARIANCE MATRIX
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
mixmeta(y, S, bscov="prop", control=list(Psifix=diag(3)+1))
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
