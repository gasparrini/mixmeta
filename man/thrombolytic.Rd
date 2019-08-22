\name{thrombolytic}
\alias{thrombolytic}
\docType{data}

\title{ Randomized Trials of Thrombolytic Therapy }

\description{
The dataset contains the data on 20 randomized trials of thrombolytic therapy, which evaluated effect on short-term mortality after a myocardial infarction (up to 35 days) in 50,246 patients in relation to treatment delay. The hypothesis is that the thrombolytic therapy reduces the mortality risk following the myocardial infarction, and that the benefit is particularly substantial for very early treatment. Some of the trials report separate results according to treatment delay, generating 38 observations from full trials or subgroups of trials. Effect sizes were reported as absolute risk reduction computed as the difference between treated and control groups in each trial or subgroup.
}

\usage{thrombolytic}

\format{
  A data frame with 38 observations on the following 10 variables:
  \itemize{
    \item{\code{trial}: label identifying the trial.}
    \item{\code{time2treat}: treatment delay after the onset of the symptoms of a myocardial infaction, reported in hours.}
    \item{\code{dtreat, ntreat}: number of deaths and total patients in the treated group, respectively.}
    \item{\code{dcontr, ncontr}: number of deaths and total patients in the control group, respectively.}
    \item{\code{risktreat, riskcontr}: risk of death in the treatment and control groups, respectively.}
    \item{\code{absrisk}: absolute risk difference of death between the treatment and control groups. See Details.}
    \item{\code{var}: variance of the absolute risk difference. See Details.}
  }
}

\details{
The absolute risk is simply the difference in risk, which is computed empirically as ratio of the number of deaths and the number of total patients in treated and control groups (\eqn{p_1=d_1/N_1}{p1=d1/N1} and \eqn{p_0=d_0/N_0}{p0=d0/N0}, respectively). The variance of the absolute risk difference is computed as \eqn{p_0(1-p_0)/N_0+p_1(1-p_1)/N_1}{p0(1-p0)/N0+p1(1-p1)/N1}. See Thompson and colleagues (2001) for details.
}

\source{
Thompson SG, Turner RM, Warn DE (2001). Multilevel models for meta-analysis, and their application to absolute risk differences. \emph{Statistical Methods in Medical Research}. \bold{10}(6):375--392.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\note{
The data provide an example of application of multilevel meta-analysis with repeated observations in an inner level within an outer level, corresponding here to treatment subgroups within each trial. This more complex correlation structure is modelled by two levels of random effects, including meta-predictors that can explain part of the heterogeneity at each level. Results can be compared with those reported by Thompson and colleagues (2001).
}

\examples{
### REPRODUCE THE RESULTS IN THOMPSON ET AL (2001), TABLES 2, 3, AND 4

# STANDARD FIXED-EFFECTS META-ANALYSIS
mod1 <- mixmeta(absrisk, var, data=thrombolytic, method="fixed")
print(summary(mod1), digits=5)

# STANDARD RANDOM-EFFECTS META-ANALYSIS
subtrial <- seq(nrow(thrombolytic))
mod2 <- mixmeta(absrisk, var, random= ~ 1|subtrial, data=thrombolytic)
print(summary(mod2), digits=5)

# TWO-LEVEL RANDOM-EFFECTS META-ANALYSIS
mod3 <- mixmeta(absrisk, var, random= ~ 1|trial/subtrial, data=thrombolytic)
print(summary(mod3), digits=5)

# TWO-LEVEL RANDOM-EFFECTS META-REGRESSION
mod4 <- mixmeta(absrisk~time2treat, var, random= ~ 1|trial/subtrial,
  data=thrombolytic)
print(summary(mod4), digits=5)

# TWO-LEVEL RANDOM-EFFECTS META-REGRESSION WITH NON-LINEAR TERM
mod5 <- mixmeta(absrisk ~ time2treat + I(1/time2treat), var,
  random= ~ 1|trial/subtrial, data=thrombolytic)
print(summary(mod5), digits=5)

### SEE help(school) FOR A COMPLEMENTARY EXAMPLE
}

\keyword{datasets}
