\name{bcg}
\alias{bcg}
\docType{data}

\title{ Efficacy of BCG Vaccine in the Prevention of Tuberculosis }

\description{
The dataset contains the data on 13 prospective clinical trials that compared the rates of tubercolosis in groups vaccinated with the Bacillus Calmette-Guerin (BCG) vaccine and non-vaccinated control populations. The outcome here is reported as both relative risk (RR) and odds ratio (OR), with associated uncertanty.
}

\usage{bcg}

\format{
  A data frame with 13 observations on the following 13 variables:
  \itemize{
    \item{\code{trial}: sequence identifying the trial.}
    \item{\code{author}: label identifying the author(s).}
    \item{\code{year}: year of publication.}
    \item{\code{tpos, tneg}: number of positive and negative TB cases in the treated (vaccinated) group.}
    \item{\code{cpos, cneg}: number of positive and negative TB cases in the control (non-vaccinated) group.}
    \item{\code{ablat}: absolute latitude of the study location (in degrees).}
    \item{\code{alloc}: method of treatment allocation (random, alternate, or systematic assignment).}
  }
}

\source{
van Houwelingen HC, Arends LR, Stijnen T (2002). Advanced methods in meta-analysis: multivariate approach and meta-regression. \emph{Statistics in Medicine}. \bold{21}(4):589--624.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\note{
The data provide an example of application of standard univariate meta-analysis and meta-regression, with independent studies providing a single estimate of a single effect size. Interestingly, the data can be analyzed also as a multivariate meta-analysis, using a bivariate outcome where risks or odds of TB can be measured separaterly in treatment and control groups. Results can be compared with those reported van Houwelingen, Arends, and Stijnen (2002).

The dataset is also available in the same format in the dataframe \code{dat.colditz1994} of the package \pkg{metafor}.
}

\examples{
### REPRODUCE THE RESULTS IN VAN HOUWELINGEN ET AL (2002)

# FIXED-EFFECTS META-ANALYSIS (SECTION 3.1.1)
unifix <- mixmeta(logor, logorvar, data=bcg, method="fixed")
print(summary(unifix), digits=3)

# RANDOM-EFFECTS META-ANALYSIS WITH MAXIMUM LIKELIHOOD (SECTION 3.1.2)
uniran <- mixmeta(logor, logorvar, data=bcg, method="ml")
print(summary(uniran), digits=3, report="var")

# ORIGINAL ESTIMATES AND BEST-LINEAR UNBIASED PREDICTIONS (FIGURE 3)
pred <- with(bcg, cbind(logor, logor-1.96*sqrt(logorvar),
  logor+1.96*sqrt(logorvar)))
blup <- blup(uniran, pi=TRUE)
plot(pred[,1], rev(bcg$trial)+0.2, xlim=c(-3,3), ylim=c(0,14), pch=18,
  axes=FALSE, xlab="Log odds ratio", ylab="Trial", main="Forest plot")
axis(1)
axis(2, at=bcg$trial, labels=rev(bcg$trial), lty=0, las=1)
abline(v=coef(uniran))
segments(pred[,2], rev(bcg$trial)+0.2, pred[,3], rev(bcg$trial)+0.2, lty=5)
points(blup[,1], rev(bcg$trial)-0.2, pch=19)
segments(blup[,2], rev(bcg$trial)-0.2, blup[,3], rev(bcg$trial)-0.2)

# COMPUTE THE OUTCOME SEPARATELY FOR TREATMENT AND CONTROL GROUPS
y <- with(bcg, log(cbind(tpos/tneg, cpos/cneg)))
S <- with(bcg, cbind(1/tpos+1/tneg, 1/cpos+1/cneg))

# BIVARIATE RANDOM-EFFECTS META-ANALYSIS (SECTION 4)
mvran <- mixmeta(y, S, method="ml")
print(summary(mvran), digits=3, report="var")

# META-REGRESSION (SECTION 5)
uniranlat <- update(uniran, .~. + ablat)
print(summary(uniranlat), digits=3, report="var")
drop1(uniranlat, test="Chisq")
}

\keyword{datasets}
