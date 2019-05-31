\name{fibrinogen}
\alias{fibrinogen}
\docType{data}

\title{ Fibrinogen Studies Collaboration }

\description{
The Fibrinogen Studies Collaboration is a meta-analysis of individual data on 154,012 adults from 31 prospective cohort studies with information on plasma fibrinogen and major disease outcomes. The dataset reports a subset of the results of a first-stage analysis consisting of the log-hazard ratio of coronary heart disease for categories of levels of fibrinogen versus a baseline category.
}

\usage{fibrinogen}

\format{
  A data frame with 31 observations on the following 15 variables:
  \itemize{
    \item{\code{cohort}: study ID.}
    \item{\code{b2, b3, b4, b5}: estimated log-hazard ratios for the second to fifth categories versus the baseline category.}
    \item{\code{V_2_2, V_3_3, V_4_4, V_5_5}: variances of the estimated log-hazard ratios.}
    \item{\code{V_2_3, V_2_4, V_2_5, V_3_4, V_3_5, V_4_5}: covariances of the estimated log-hazard ratios.}
  }
}

\details{
The published analysis adopted a fixed-effects model on 10 categories of fibrinogen (Fibrinogen Studies Collaboration 2004, 2005). Here a subset of the results of the first-stage analysis is reported, namely the log-hazard ratio for 4 categories and associated (co)variance terms, ordered as the lower triangular elements of the (co)variance matrix taken by column. Details on the first-stage model and the second-stage meta-analysis are provided in White (2009) and Jackson and colleagues (2010).
}

\note{
The data provide an example of application of multivariate meta-analysis for multi-parameter association, where a relationship is defined by functions specified by several coefficients. In this case, the coefficients refer to log-hazard ratio for strata of the original variable versus a baseline category. A general overview of the application of multivariate meta-analysis in this setting is provided by Gasparrini and colleagues (2012).
}

\source{
Fibrinogen Studies Collaboration (2004). Collaborative meta-analysis of prospective studies of plasma fibrinogen and cardiovascular disease. \emph{European Journal of Cardiovascular Prevention and Rehabilitation}. \bold{11}:9--17.

Fibrinogen Studies Collaboration (2005). Plasma fibrinogen level and the risk of major cardiovascular diseases and nonvascular mortality: an individual participant meta-analysis. \emph{Journal of the American Medical Association}. \bold{294}:1799--1809.

White IR (2009). Multivariate random-effects meta-analysis. \emph{Stata Journal}. \bold{9}(1):40--56.

Jackson D, White IR, Thompson SG (2010). Extending DerSimonian and Laird's methodology to perform multivariate random effects meta-analyses. \emph{Statistics in Medicine}. \bold{29}(12):1282--1297.

Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].
}

\examples{
### REPRODUCE THE RESULTS IN WHITE (2009) AND JACKSON ET AL. (2010)

# INSPECT THE DATA
head(fibrinogen)

# REML MODEL
y <- as.matrix(fibrinogen[2:5])
S <- as.matrix(fibrinogen[6:15])
model <- mixmeta(y, S)

# SUMMARIZE THE RESULTS
print(summary(model), digits=3)
round(model$Psi, 3)
}

\keyword{datasets}
