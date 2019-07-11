\name{mixmeta.fixed}
\alias{mixmeta.fixed}


\title{ Fixed-Effects Estimator for mixmeta Models }

\description{
This function implements a generalized least square estimator for fixed-effects meta-analysis and meta-regression, including standard univariate models and non-standard multivariate extensions. It is meant to be used internally and not directly run by the users.
}

\usage{
mixmeta.fixed(Xlist, ylist, Slist, nall, control, ...)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{m} independent groups (usually studies) providing a single estimate (unit of analysis), \eqn{k} outcomes and \eqn{p} fixed-effects predictors:

  \item{Xlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of group-specific vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-group (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nall }{ numeric scalar with the total number of observations (excluding missing).}
  \item{control }{ list of parameters for controlling the fitting process, usually internally set to default values by \code{\link{mixmeta.control}}.}
  \item{\dots }{ further arguments passed to or from other methods. Currently not used.}
}

\details{
The estimation involves only the \eqn{kp} fixed-effects coefficients. Note that, in this fixed-effects estimator, each unit is assumed independent from the others, and therefore the number of groups (the length of the lists) is identical to the number of units (\code{m=n}). However, this is not important in fixed-effects models, where no random (and therefore grouping) structure is used.

The routine is based on a standard generalized least square (GLS) algorithm implemented in the internal function \code{glsfit}. The between-study (co)variance matrix is set to zero, so the marginal (co)variance matrix, composed only by elements of the within-unit component, is assumed as completely known. Similarly to the likelihood-based estimators implemented in  \code{\link{mixmeta.ml}} and \code{\link{mixmeta.reml}}, the computation involves Cholesky and and QR decompositions for computational stability and efficiency. The method is described in details in Gasparrini and collaborators (2012) (see references below).
}

\value{
These functions return an intermediate list object, with some components then processed and some others added later within \code{\link{mixmeta.fit}} and \code{\link{mixmeta}} to finalize an object of class \code{"mixmeta"}. See \code{\link{mixmetaObject}}.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. In submission.

Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Berkey CS, Anderson JJ, Hoaglin DC (1996). Multiple-outcome meta-analysis of clinical trials. \emph{Statistics in Medicine}. \bold{15}(5):537--547.

Berkey CS, et al. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}. \bold{17}(22):2537--2550.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
As stated earlier, this function is called internally by \code{\link{mixmeta.fit}}, and is not meant to be used directly. In particular, its code does not contain any check on the arguments provided, which are expected in specific formats. The function is however exported in the namespace and documented for completeness.

The arguments above are prepared by \code{\link{mixmeta.fit}} from its arguments \code{X}, \code{y} and \code{S}. The list structure, although requiring more elaborate coding, is computationally more efficient, as it avoids the specification of sparse block-diagonal matrices, especially for meta-analysis involving a large number of studies.

Some parameters of the fitting procedures are determined by the \code{control} argument, with default set by \code{\link{mixmeta.control}}. No missing values are accepted in the fitting functions. See details on \code{\link[=na.omit.data.frame.mixmeta]{missing values}}.
}

\seealso{
See \code{\link{mixmeta}} for the general usage of the functions. See \code{\link{mixmeta.control}} to determine specific parameters of the fitting procedures. Use the triple colon operator ('\code{:::}') to access the code of the internal functions, such as \code{glsfit}. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# UNIVARIATE FIXED-EFFECTS MODEL
mod1 <- mixmeta(yC, S=SCC, data=smoking, method="fixed")
summary(mod1)

# MULTIVARIATE FIXED-EFFECTS MODEL
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
mod2 <- mixmeta(y, S, method="fixed")
summary(mod2)

# MULTIVARIATE FIXED-EFFECTS MODEL: REPLICATE THE RESULTS IN BERKEY ET AL. 1998
mod3 <- mixmeta(cbind(PD,AL) ~ I(pubyear-1983), S=berkey98[5:7], data=berkey98,
  method="fixed")
summary(mod3)
}


\keyword{models}
\keyword{regression}
\keyword{multivariate}
