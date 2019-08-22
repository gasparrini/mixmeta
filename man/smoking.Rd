\name{smoking}
\alias{smoking}
\docType{data}

\title{ Meta-Analysis of Interventions to Promote Smoking Cessation }

\description{
The dataset contains the results of 24 trials comparing four alternative interventions to promote smoking cessation. The trials have different designs, comparing two or three different interventions. The data consist of the number of successes out of the total participants, and the estimated log-odds ratio for arms B (self-help), C (individual counselling), and D (group counselling) relative to arm A (no contact), as well as the (co)variance matrix of these three estimates.
}

\usage{smoking}

\format{
  A data frame with 24 observations on the following 19 variables:
  \itemize{
    \item{\code{study}: study ID.}
    \item{\code{design}: design of the trial, reporting the interventions being compared.}
    \item{\code{dA, dB, dC, dD}: number of successes for each intervention.}
    \item{\code{nA, nB, nC, nD}: number of participants for each intervention.}
    \item{\code{yB, yC, yD}: estimated log-odds ratios for interventions B, C and D versus intervention A.}
    \item{\code{SBB, SBC, SBD, SCC, SCD, SDD}: variances and co-variances of the estimated log-odds ratios for interventions B, C and D versus intervention A. The order corresponds to the lower triangular elements of the (co)variance matrix taken by column.}
  }
}

\details{
Intervention A is chosen as the reference category. Trials without an arm A (trials 2 and 21-24) are augmented with an arm A with 0.01 individuals and 0.001 successes. Trials containing zero cells (trials 9 and 20) have 1 individual with 0.5 successes added to each intervention. Details on the data augmentation and estimation of (co)variances of the log-odds ratios are provided by White (2011).
}

\note{
The data provide an example of application of network meta-analysis, also referred to as indirect mixed-treatment comparison. Additional information using examples based on these data are provided by Lu and Ades (2006), White (2011) and Higgins and colleagues (2012).
}

\source{
Lu G and Ades AE (2006). Assessing evidence inconsistency in mixed treatment comparisons. \emph{Journal of the American Statistical Association}. \bold{101}:447--459.

Higgins JPT, et al. (2012). Consistency and inconsistency in network meta-analysis: concepts and models for multi-arm studies. \emph{Research Synthesis Methods}. \bold{3}(2):98--110.

White IR (2011). Multivariate random-effects meta-regression. \emph{The Stata Journal}. \bold{11}:255--270.

Sera F, Gasparrini A (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. DOI: 10.1002/sim.8362 [Freely available \href{http://www.ag-myresearch.com/2019_sera_statmed.html}{\bold{here}}].
}

\examples{
### REPRODUCE THE RESULTS IN WHITE (2011)

# INSPECT THE DATA
head(smoking)
names(smoking)

# CONSISTENCY MODEL, UNSTRUCTURED BETWEEN-STUDY (CO)VARIANCE
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
mod1 <- mixmeta(y, S)
summary(mod1)

# CONSISTENCY MODEL, STRUCTURED BETWEEN-STUDY (CO)VARIANCE (PROPORTIONAL)
mod2 <- mixmeta(y, S, bscov="prop", control=list(Psifix=diag(3)+1))
summary(mod2)

# TRANSFORM IN LONG FORMAT, WITH S AS LIST (EXCLUDING MISSING)
long <- na.omit(reshape(smoking[,c(1,2,11:13)], varying=list(3:5), idvar="study", 
  v.names="y", timevar="outcome", times=colnames(y), direction="long"))
Slist <- lapply(lapply(seq(nrow(S)), function(i) xpndMat(S[i,])), function(x)
  x[!is.na(diag(x)), !is.na(diag(x)), drop=FALSE])

# THE MODELS ABOVE CAN BE REPLICATED IN THE LONG FORMAT
mod2b <- mixmeta(y ~ 0 + factor(outcome), random= ~ 0 + factor(outcome)|study,
  data=long, bscov="prop", control=list(addS=Slist, Psifix=diag(3)+1))
summary(mod2b)

# DEFINE AND ADD INDICATORS FOR OUTCOME AND DESIGN
dummy <- cbind(model.matrix(~0+outcome, long), model.matrix(~0+design, long))
colnames(dummy) <- c(levels(factor(long$outcome)), levels(long$design))
long <- cbind(long, data.frame(dummy))

# INCONSISTENCY MODEL (SPECIAL PARAMETERIZATION OF OUTCOME-BY-DESIGN INTERACTION)
formula <- y ~ 0 + yB + yC + yC:acd + yC:bc + yC:bcd + yD + yD:acd + yD:bcd + 
  yD:bd + yD:cd
mod3 <- update(mod2b, formula=formula)
summary(mod3)
}

\keyword{datasets}
