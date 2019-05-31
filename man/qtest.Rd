\name{qtest}
\alias{qtest}


\title{ Cochran Q Test of Heterogeneity }

\description{
This is a generic function to perform a Cochran Q test of (residual) heterogeneity. The function invokes particular \code{\link{methods}} which depend on the \code{\link{class}} of the first argument.
}

\usage{
qtest(object, \dots)
}

\arguments{

  \item{object }{ an object for which the test is desired}
  \item{\dots}{ further arguments passed to specific methods.}
}

\details{
A specific method function \code{\link{qtest.mixmeta}} has been provided for objects of class \code{"mixmeta"}.
}

\value{
Returned values depend on the specific class. Usually, the results of the test.
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}>}

\seealso{
See \code{\link{qtest.mixmeta}} for the specific method function applied to objects of class \code{"mixmeta"}.
}

\keyword{methods}
\keyword{htest}
