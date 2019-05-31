\name{qtest.mixmeta}
\alias{qtest.mixmeta}
\alias{print.qtest.mixmeta}


\title{ Cochran Q Test of Heterogeneity for mixmeta Models }

\description{
This method function performs a Cochran Q test of (residual) heterogeneity on fitted meta-analytical models represented in objects of class \code{"mixmeta"}.
}

\usage{
\method{qtest}{mixmeta}(object, \dots)

\method{print}{qtest.mixmeta}(x, digits=3, \dots)
}

\arguments{

  \item{object, x}{ objects of classes \code{"mixmeta"} and \code{"qtest.mixmeta"}, respectively.}
  \item{digits }{  an integer specifying the number of digits to which printed results must be rounded.}
  \item{\dots}{ further arguments passed to or from other methods.}

}

\details{
In multivariate models, the test assesses the null hypothesis that the variability in the multivariate distribution of the outcomes is explained only in terms of estimation error in each unit, measured by the within-unit (co)variance matrices stored in the component \code{S} of \code{mixmeta} objects. This is equal to test the hypothesis that the random-effects (co)variance matrix is a zero matrix, and there is no random deviation in unit-specific estimates. Tests for single outcome parameters, comparable to estimates from multiple univariate meta-analysis, are also reported. This test reduces to the standard Q test in univariate models.

The function compute the statistics by actually fitting the related fixed-effects model, re-evaluating the call of the model with \code{method} changed to \code{"fixed"}.
}

\value{
A list object of class \code{"qtest.mixmeta"} with the following components:

  \item{Q }{ the vector of test statistics for overall and outcome-specific tests, distributed under the null hypothesis as a Chi-square with degrees of freedom \code{df}.}
  \item{df }{ the vector of degrees of freedom of the null distribution for overall and outcome-specific tests. For the overall test, equal to the number of observations used for estimation minus the number of coefficients in the fixed part of the model. For outcome-specific test, equal to number of observed values minus the number of coefficients.}
  \item{pvalue }{ the vector of p-values for overall and outcome-specific tests.}
  \item{residual }{ logical switch indicating if a meta-regression model is assessed, meaning that the tested heterogeneity is residual.}
  \item{k }{ dimensionality of the overall test, that is the number of outcome parameters in the model.}

As usual, the \code{print} method function for class \code{"qtest.mixmeta"} does not return any value.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].

Berkey, et al. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}. \bold{17}(22):2537--2550.

Ritz J, Demidenkob E, Spiegelman G (2008). Multivariate meta-analysis for data consortia, individual patient meta-analysis, and pooling projects. \emph{Journal of Statistical Planning and Inference}. \bold{139}(7):1919--1933.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\note{
Tests on single outcome parameters are performed by extracting the related estimates and variances, but they do not account for the correlation between them, which nevertheless has been considered in estimation. These tests are not therefore comparable with those performed by running a univariate model on each outcome parameter.
}

\seealso{
See \code{\link{qtest}} for the generic method function. See \code{\link{mixmeta-package}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL
model <- mixmeta(cbind(PD,AL) ~ 1, S=berkey98[5:7], data=berkey98)

# MULTIVARIATE COCHRAN Q TEST FOR HETEROGENEITY
test <- qtest(model)
print(test, digits=2)
unclass(test)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
\keyword{htest}
