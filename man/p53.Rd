\name{p53}
\alias{p53}
\docType{data}

\title{ Mutant p53 Gene and Squamous Cell Carcinoma }

\description{
The dataset includes studies providing evidence about whether the presence of mutant p53 tumour suppressor gene is a prognostic factor for patients presenting with squamous cell carcinoma arising from the oropharynx cavity. Unadjusted estimates of log hazard ratios of mutant p53 to normal p53 for disease-free and overall survival, together with the associated variances, are collected from 6 observational studies.
}

\usage{p53}

\format{
  A data frame with 6 observations on the following 5 variables:
  \itemize{
    \item{\code{study}: study ID.}
    \item{\code{y1, V1}: estimate and associated variance of the log hazard ratio for disease-free survival.}
    \item{\code{y2, V2}: estimate and associated variance of the log hazard ratio for overall survival.}
  }
}

\details{
Only 3 studies provide estimates for disease-free survival. The within-study correlations are not reported in the original studies but are expected to be highly positively correlated. The original data are described in Tandon and colleagues (2010) and used as an example by Jackson and colleagues (2011).
}

\note{
The data provide an example of application of multivariate meta-analysis when the within-study correlations are not known. These correlations can be inputted directly in the \code{\link{mixmeta}} function through the \code{control} argument. See \code{\link{mixmeta.control}} for details.
}

\source{
Jackson D, Riley R, White IR (2011). Multivariate meta-analysis: Potential and promise. \emph{Statistics in Medicine}. \bold{30}(20);2481--2498.

Tandon S, Tudur-Smith C, Riley RD, et al. (2010). A systematic review of p53 as a prognostic factor of survival in squamous cell carcinoma of the four main anatomical subsites of the head and neck. \emph{Cancer Epidemiology, Biomarkers and Prevention}. \bold{19}(2):574--587.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\examples{
### REPRODUCE THE RESULTS OF EXAMPLE 3 IN JACKSON ET AL. (2011)

# INSPECT THE DATA
p53

# REML MODEL WITH INPUTTED CORRELATION EQUAL TO 0.95
model <- mixmeta(cbind(y1,y2), cbind(V1,V2), data=p53, control=list(Scor=0.95))
print(summary(model), digits=2)
}

\keyword{datasets}
