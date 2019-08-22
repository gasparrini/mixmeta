\name{school}
\alias{school}
\docType{data}

\title{ Studies on Modified School Calendar and Student Achievement }

\description{
The dataset contains the results of 56 studies that evaluate the effect of a modified school calendar on student achievement. The studies assessed students from grade 1 to 9 and reported standardized reading achievement differences between schools that follow a year-round versus the traditional nine-month calendar. The studies were performed in separate school districts, with at least three studies in each district.
}

\usage{school}

\format{
  A data frame with 56 observations on the following 5 variables:
  \itemize{
    \item{\code{district, study}: numbers identifying the school district and study, respectively.}
    \item{\code{effect}: estimated standardized effect, reported as difference in reading achievement expressed in standard deviation units.}
    \item{\code{var}: within-study variance of the estimated effects.}
    \item{\code{year}: year when the study was performed.}
  }
}

\source{
Kostantopoulos S (2011). Fixed effects and variance components estimation in three-level meta-analysis. \emph{Research Synthesis Methods}. \bold{2}(1):61--76.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\note{
The data provide an example of application of multilevel meta-analysis with multiple nested random-effects levels, where effect sizes are correlated between studies within school district. This more complex correlation structure is modelled by two levels of random effects. Results can be compared with the so-called three-level model in Kostantopoulos (2011), that is defined as a two-level meta-analysis here.
}

\examples{
### REPRODUCE THE RESULTS IN KOSTANTOPOULOS (2011), TABLES 4 AND 5

# STANDARD META-ANALYSIS (NB: random NOT STRICTLY NEEDED HERE)
mod1 <- mixmeta(effect, var, random= ~ 1|study, data=school, method="ml")
print(summary(mod1), digits=3, report="var")

# STANDARD META-REGRESSION
yearcen <- school$year - mean(school$year)
mod2 <- mixmeta(effect ~ yearcen, var, random= ~ 1|study, data=school,
  method="ml")
print(summary(mod2), digits=3, report="var")

# TWO-LEVEL META-ANALYSIS
mod3 <- mixmeta(effect, var, random= ~ 1|district/study, data=school,
  method="ml")
print(summary(mod3), digits=3, report="var")

# TWO-LEVEL META-REGRESSION
yearcen2 <- with(school, year - mean(tapply(year, district, mean)))
mod4 <- mixmeta(effect ~ yearcen2, var, random= ~ 1|district/study, data=school,
  method="ml")
print(summary(mod4), digits=3, report="var")

### SEE help(thrombolytic) FOR A COMPLEMENTARY EXAMPLE
}

\keyword{datasets}
