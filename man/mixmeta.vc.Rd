\name{mixmeta.vc}
\alias{mixmeta.vc}


\title{ Variance Components Estimator for mixmeta Models }

\description{
This function implements a variance components estimator for multivariate and univariate  random-effects meta-analysis and meta-regression. It is meant to be used internally and not directly run by the users.
}

\usage{
mixmeta.vc(Xlist, ylist, Slist, nalist, k, m, p, nall, control, ...)
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
In this variance components estimator, only a basic random-effects structure is allowed, where each group (usually corresponding to an independent study) provides a single estimate (unit of analysis) for one or multiple outcomes. This implies that the number of groups (\emph{i.e.}, the length of the lists) is identical to the number of units (\code{m=n}). In addition, only an unstructured form for the(co)variance matrix of the single level of random effects is permitted. Therefore, the estimation involves \eqn{kp} fixed-effects coefficients and \eqn{k(k+1)/2} random-effects parameters, corresponding to the lower triangular entries of the between-study (co)variance matrix.

The procedure is based on the estimate of the between-group (co)variance as the difference between the marginal (co)variance and the average within-group (co)variance. This in turn requires the estimate of the marginal (co)variance, obtained by the residuals of the fitted model. The procedure is iterative, with the current estimate of the between-group (co)variance plugged into a generalized least square (GLS) routine. Starting values are provided by a fixed-effects estimator (see \code{\link{mixmeta.fixed}}). The algorithm is fast and generally converges with few iterations.

Similar versions of this estimator has been previously proposed. Berkey and collaborators (1998) simply called it GLS method, and a non-iterative approach was proposed by Ritz and collaborators (2008), referred to as MVEE3 in their article. A non-iterative version for univariate models is discussed in Sidik and Jonkman (2007). The results from Berkey and collaborators (1998) are reproduced in the example below.

In the original approach, the estimate of the marginal (co)variance is obtained from the sum of the residual components using a denominator equal to \code{m-p}. Following the development proposed by Kauermann and Carroll (2001) and Fay and Graubard (2001) in the context of sandwich (co)variance estimators, then discussed by Lu and collaborators (2007), an adjusted denominator can be computed as a quantity derived from the hat matrix. This method is expected to perform better in the presence of missing values and small data sets. This alternative adjustment is chosen by default by setting \code{vc.adj=TRUE} in the \code{control} argument.

This variance component estimator is not bounded to provide a positive semi-definite between-study (co)variance matrix, as shown in the simulation study by Liu and colleagues (2009). Here positive semi-definiteness is forced by setting the negative eigenvalues of the estimated matrix to a value close to zero at each iteration (see \code{\link[=mixmeta.control]{control}}). Little is known about the impact of such constraint.
}

\value{
This function returns an intermediate list object, with some components then processed and some others added later within \code{\link{mixmeta.fit}} and \code{\link{mixmeta}} to finalize an object of class \code{"mixmeta"}. See \code{\link{mixmetaObject}}.
}

\references{
Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].

Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Ritz J, Demidenko E, Spiegelman G (2008). Multivariate meta-analysis for data consortia, individual patient meta-analysis, and pooling projects. \emph{Journal of Statistical Planning and Inference}. \bold{139}(7):1919--1933.

Berkey CS, et al. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}. \bold{17}(22):2537--2550.

Liu Q, et al (2009). A two-stage hierarchical regression model for meta-analysis of epidemiologic nonlinear dose-response data. \emph{Computational Statistics and Data Analysis}. \bold{53}(12):4157--4167

Sidik K, Jonkman JN (2007). A comparison of heterogeneity variance estimators in combining results of studies. \emph{Statistics in  Medicine}. \bold{26}(9):1964--81.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\note{
As stated earlier, this function is called internally by \code{\link{mixmeta.fit}}, and is not meant to be used directly. In particular, its code does not contain any check on the arguments provided, which are expected in specific formats. The function is however exported in the namespace and documented for completeness.

The arguments above are prepared by \code{\link{mixmeta.fit}} from its arguments \code{X}, \code{y} and \code{S}. The list structure, although requiring more elaborate coding, is computationally more efficient, as it avoids the specification of sparse block-diagonal matrices, especially for meta-analysis involving a large number of studies.

Some parameters of the fitting procedures are determined by the \code{control} argument, with default set by \code{\link{mixmeta.control}}. No missing values are accepted in the fitting functions. See details on \code{\link[=na.omit.data.frame.mixmeta]{missing values}}.
}

\seealso{
See \code{\link{mixmeta}} for the general usage of the functions. See \code{\link{mixmeta.control}} to determine specific parameters of the fitting procedures. Use the triple colon operator ('\code{:::}') to access the code of the internal functions, such as \code{sumlist}. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# VC ESTIMATOR: UNIVARIATE MODEL
mod1 <- mixmeta(PD ~ pubyear, S=berkey98[,5], data=berkey98, method="vc")
summary(mod1)

# VC ESTIMATOR: MULTIVARIATE MODEL
mod2 <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98,
  method="vc")
summary(mod2)

# VC ESTIMATOR: NON-ITERATIVE VERSION
mod3 <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98,
  method="vc", control=list(maxiter=1))
summary(mod3)

# VARIANCE COMPONENTS ESTIMATOR: REPLICATE THE RESULTS IN BERKEY ET AL. (1998)
mod4 <- mixmeta(cbind(PD,AL) ~ I(pubyear-1983), S=berkey98[5:7], data=berkey98,
  method="vc", control=list(vc.adj=FALSE))
summary(mod4)
}


\keyword{models}
\keyword{regression}
\keyword{multivariate}
