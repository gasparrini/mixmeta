\name{mixmeta.mm}
\alias{mixmeta.mm}


\title{ Method of Moments Estimator for mixmeta Models }

\description{
This function implements a method of moments estimator for multivariate and univariate  random-effects meta-analysis and meta-regression. It is meant to be used internally and not directly run by the users.
}

\usage{
mixmeta.mm(Xlist, ylist, Slist, nalist, k, m, p, nall, control, ...)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{m} independent groups (usually studies) providing a single estimate (unit of analysis), \eqn{k} outcomes and \eqn{p} fixed-effects predictors:

  \item{Xlist }{ a \eqn{m}-dimensional list of group-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of group-specific vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-group (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nalist }{ a \eqn{m}-dimensional list of group-specific logical vectors, identifying missing outcomes.}
  \item{k, m, p, nall }{ numeric scalars: number of outcomes, number of groups included in estimation (equal to the length of lists above), number of predictors (including the intercept), total number of observations (excluding missing).}
  \item{control }{ list of parameters for controlling the fitting process, usually internally set to default values by \code{\link{mixmeta.control}}.}
  \item{\dots }{ further arguments passed to or from other methods. Currently not used.}
}

\details{
In this method of moments estimator, only a basic random-effects structure is allowed, where each group (usually corresponding to an independent study) provides a single estimate (unit of analysis) for one or multiple outcomes. This implies that the number of groups (\emph{i.e.}, the length of the lists) is identical to the number of units (\code{m=n}). In addition, only an unstructured form for the(co)variance matrix of the single level of random effects is permitted. Therefore, the estimation involves \eqn{kp} fixed-effects coefficients and \eqn{k(k+1)/2} random-effects parameters, corresponding to the lower triangular entries of the between-study (co)variance matrix.

The method of moment estimator implemented here represents a multivariate extension of the traditional estimator proposed by DerSimonian and Laird (1986), and simplifies to the standard method in the univariate case. The estimator used here is described in Jackson and collaborators (2013) as a generalization of that developed by Chen and collaborators (2012). However, this general version is computationally more intensive, and may turn out to be slow when applied to meta-analysis of a relatively high number of studies. An alternative and computationally faster method of moment estimator was previously proposed by Jackson and collaborators (2010), although it is not invariant to reparameterization. This latter estimator is currently not implemented in \pkg{mixmeta}. See references below.

This method of moments estimator is not bounded to provide a positive semi-definite random-effects (co)variance matrix, as shown in the simulation study by Liu and colleagues (2009). Here positive semi-definiteness is forced by setting the negative eigenvalues of the estimated matrix to a positive value close to zero at each iteration (see \code{\link[=mixmeta.control]{control}}). Little is known about the impact of such constraint.
}

\value{
This function returns an intermediate list object, with some components then processed and some others added later within \code{\link{mixmeta.fit}} and \code{\link{mixmeta}} to finalize an object of class \code{"mixmeta"}. See \code{\link{mixmetaObject}}.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Jackson D, White IR, Riley RD (2013). A matrix based method of moments for fitting the multivariate random effects model for meta-analysis and meta-regression. \emph{Biometrical Journal}. \bold{55}(2):231--45.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\note{
As stated earlier, this function is called internally by \code{\link{mixmeta.fit}}, and is not meant to be used directly. In particular, its code does not contain any check on the arguments provided, which are expected in specific formats. The function is however exported in the namespace and documented for completeness.

The arguments above are prepared by \code{\link{mixmeta.fit}} from its arguments \code{X}, \code{y} and \code{S}. The list structure, although requiring more elaborate coding, is computationally more efficient, as it avoids the specification of sparse block-diagonal matrices, especially for meta-analysis involving a large number of studies.

Some parameters of the fitting procedures are determined by the \code{control} argument, with default set by \code{\link{mixmeta.control}}. No missing values are accepted in the fitting functions. See details on \code{\link[=na.omit.data.frame.mixmeta]{missing values}}.
}

\seealso{
See \code{\link{mixmeta}} for the general usage of the functions. See \code{\link{mixmeta.control}} to determine specific parameters of the fitting procedures. Use the triple colon operator ('\code{:::}') to access the code of the internal functions, such as \code{fbtr}. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# MM ESTIMATOR: UNIVARIATE MODEL
mod1 <- mixmeta(PD ~ pubyear, S=berkey98[,5], data=berkey98, method="mm")
summary(mod1)

# MULTIVARIATE MODEL: REPRODUCE THE RESULTS IN CHEN ET AL. (2012)
S <- as.matrix(hsls[5:10])
mod2 <- mixmeta(cbind(b1,b2,b3), S, data=hsls, method="mm")
summary(mod2)

# MULTIVARIATE MODEL: REPRODUCE THE RESULTS IN JACKSON ET AL. (2013)
S <- inputcov(hyp[c("sbp_se","dbp_se")], cor=hyp$rho)
mod3 <- mixmeta(cbind(sbp,dbp), S=S, data=hyp, method="mm")
summary(mod3)
}


\keyword{models}
\keyword{regression}
\keyword{multivariate}
