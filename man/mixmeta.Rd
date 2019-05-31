\name{mixmeta}
\alias{mixmeta}
\alias{mixmeta.fit}


\title{ Fitting Standard and Extended Meta-Analysis and Meta-Regression Models }

\description{
The function \code{mixmeta} performs various meta-analytical models under a common mixed-effects framework, including standard univariate fixed and random-effects meta-analysis and meta-regression, and non-standard extensions such as multivariate, multilevel, longitudinal, and dose-response models. The function \code{mixmeta.fit} is a wrapper for actual fitting functions based on different estimation methods, usually called internally. See \code{\link{mixmeta-package}} for an overview.
}

\usage{
mixmeta(formula, S, data, random, method="reml", bscov="unstr", offset, subset,
  contrasts=NULL, na.action, model=TRUE, control=list())

mixmeta.fit(X, Z, y, S, groups, method, bscov, control)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{n} units aggregated within \eqn{m} (outer-level) groups, \eqn{k} outcomes, \eqn{p} fixed-effects predictors, and \eqn{q} random-effects predictors:

  \item{formula }{ an object of class \code{"\link{formula}"} (or one that can be coerced to that class) offering a symbolic description of the linear predictor for the fixed-effects part of the model. Alternatively, for meta-analysis with no fixed-effects predictors, a single  vector (for univariate models) or matrix-type object (for multivariate models). Terms in \code{formula} must be vector or matrix-type objects, optionally provided in the \code{data} argument below. See \code{\link{mixmetaFormula}}.}
  \item{S }{ series of within-unit variances (or (co)variance matrices for multivariate models) of the estimated outcome(s). For univariate models, this is usually a \eqn{n}-dimensional vector. For multivariate models, it can be provided as: a \eqn{m}-dimensional list of \eqn{k \times k}{k x k} matrices; a tri-dimensional \eqn{k \times k \times m}{k x k x n} array; a matrix or data frame with \eqn{n} rows and \eqn{k(k+1)/2} or \eqn{k} columns, depending on the availability of the within-unit correlations. \code{mixmeta.fit} accepts only the last option. Optionally, for more complex error structures, this argument can be omitted and passed through \code{addSlist} in \code{control}. See Details below.}
  \item{data }{ an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in \code{formula} and \code{random}. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{mixmeta} is called.}
  \item{random }{ a one-sided formula (or a list of formulae for multilevel models) offering a symbolic description of the linear predictor(s) and grouping structure for the random-effects part of the model. The usual form is \code{~ z1 + ... + zq | g}, with the grouping factor separated from the linear predictor by the symbol \code{'|'}. Multiple levels with the same linear predictor can be defined by separating multiple grouping factors using the symbol \code{'/'}. Alternatively, in a list form the grouping factors can also be provided as list names. In both cases, the levels are considered nested (from outer to inner following the order). See \code{\link{mixmetaFormula}} and Details below.}
  \item{method }{ estimation method: \code{"fixed"} for fixed-effects models, \code{"ml"} or \code{"reml"} for random-effects models fitted through (restricted) maximum likelihood, \code{"mm"} for random-effects models fitted through method of moments, and \code{"vc"} for random-effects models fitted through variance components. See Details below. If \code{"model.frame"}, the model frame is returned, as in \code{\link{lm}} or \code{\link{glm}}.}
  \item{bscov }{ a character vector defining the structure of the random-effects (co)variance matrices. Default to \code{"unstr"} (unstructured). Names corresponding to grouping factors (see \code{random} above) can be used to refer to specific random-effects levels for non-default values. If unnamed, the values can be recycled. Among various \code{\link[=mixmetaCovStruct]{(co)variance structures}}, the user can select \code{"diag"} (diagonal), \code{"cs"} (compound symmetry), \code{"hcs"} (heterogeneous compound symmetry), \code{"ar1"} (autoregressive of first order), or \code{"fixed"} (fixed). See also Details.}
  \item{offset }{ optionally, a \eqn{n}-dimensional numeric vector used to specify an a priori known component in the linear predictor. One or more \code{\link{offset}} terms can be included in the formula instead or as well. See \code{\link{model.offset}}.}
  \item{subset }{ an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{contrasts }{ an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix}}.}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s. Default to \code{na.action} setting of \code{\link{options}}, usually \code{\link[=na.omit.data.frame.mixmeta]{na.omit}}. \code{\link[=na.exclude.data.frame.mixmeta]{na.exclude}} can be useful. See details on \code{\link[=na.omit.data.frame.mixmeta]{missing values}} in \code{mixmeta}.}
  \item{model }{ a logical value indicating whether the model frame should be included as a component of the returned value. See the \code{\link[=model.frame.mixmeta]{model.frame}} method function.}
  \item{control }{ list of parameters for controlling the fitting process. These are passed to  \code{\link{mixmeta.control}} to replace otherwise selected default values.}
  \item{X }{ a \eqn{n \times p}{n x p} design matrix containing the \eqn{p} fixed-effects predictors, appropriately ordered by \code{groups}. Usually produced internally by \code{mixmeta} from \code{formula} above.}
  \item{Z }{ a \eqn{n \times q}{n x q} design matrix (or a list of design matrices for multilevel models) containing the \eqn{q} random-effects predictors, appropriately ordered by \code{groups}. Usually produced internally by \code{mixmeta} from \code{random} above.}
  \item{y }{ a \eqn{n}-dimensional vector (for univariate models) or \eqn{m \times k}{m x k} matrix (for multivariate models) of outcomes, appropriately ordered by \code{groups}. Usually produced internally by \code{mixmeta} from \code{formula} above.}
  \item{groups }{ matrix with \eqn{n} rows, with each column identifying the groups for each level of random-effects, appropriately ordered. Usually produced internally by \code{mixmeta} from \code{random} above.}
}

\details{
The function \code{mixmeta} resembles standard regression functions in \R. See \code{\link[nlme]{lme}} in particular, or \code{\link{lm}} or \code{\link{glm}}, for information on most of the arguments. Internally, this function assembles the data components, defines the (potentially multiple) grouping levels and re-order the data accordingly, and then pass them to \code{mixmeta.fit}. This is a wrapper for actual fitting functions that are automatically selected. Functions other than \code{mixmeta} are not expected to be called directly for model fitting.

Fixed or random-effects models for meta-analysis are simply defined using \code{y ~ 1} in \code{formula}, where \code{y} is a response vector optionally stored in \code{data}. In  meta-regression models, other terms are added in the right-hand side of the formula as \code{y ~ x1 + ... + xp}, defining the linear meta-predictor. Factors, variable transformations and interactions are allowed, following the usual formula specification (see \code{\link{mixmetaFormula}}).

In this standard usage, each of the \eqn{n} rows is assumed to represent a single estimate of an outcome from a set of independent studies. In random-effects models, the grouping structure is automatically derived by assigning a group to each row of data (with \eqn{m=n}). Extensions to multivariate models (\eqn{k>1}) are straightforward, and only require using a matrix in the left-hand side, where each of the \eqn{k} columns represents a different outcome, or the form \code{cbind(y1, ..., yk) ~ 1}. See \code{\link{mixmetaFormula}}.

Non-standard random-effects models can be specified through the optional argument \code{random}. This is commonly represented by a one-sided formula, whose basic random-intercept form is \code{~ 1 | g}, where \code{g} is a grouping factor. A more complex linear meta-predictor for the random-effects part can be also specified by \code{~ z1 + ... + zq | g}.  The argument \code{random} also accepts a list of one-sided formulae referring to multiple random-effects levels (see \code{\link{mixmetaFormula}}). The use of \code{random} extends the standard meta-analytical setting by relaxing the assumption of independence between units, allowing multiple estimates from the same group (with \eqn{m<n}) and multiple nested grouping levels. This provides the possibility to fit longitudinal, multilevel, and dose-response meta-analysis, among other extensions. See the examples below.

The argument \code{bscov} allows the definition of specific structures for the random-effects (co)variance matrices corresponding the each level. The default unstructured form requires \code{kq(kq+1)/2} parameters for a single-level meta-analysis. The choice of other structures reduces the number of parameters, although requiring stronger assumptions. More information and complete list of options is available at a specific help page (see \code{\link{mixmetaCovStruct}}).

The within-unit (co)variances are provided through the argument \code{S}, usually as a matrix. If the correlations are available, each of the \eqn{m} row represents the \eqn{k(k+1)/2} vectorized entries of the lower triangle of the related (co)variance matrix, taken by column (see \code{\link{xpndMat}}). If correlations are not available, each row represents the \eqn{k} variances, and the correlations are inputted internally through the argument \code{Scor} of the \code{control} list (see \code{\link{inputcov}}). For more complex error structures that span multiple units, the argument \code{S} can be omitted and passed through \code{addSlist} in \code{control}, although requiring the observations to be re-ordered accordingly to groups (see \code{\link{mixmeta.control}}).

Different estimator are available in the package \code{mixmeta} and chosen through the argument \code{method}, with related fitting functions called internally. In the current version, the options are:

\itemize{
  \item \code{method="fixed"}: \code{\link[=mixmeta.fixed]{Fixed-effects estimator}}
  \item \code{method="ml"}: \code{\link[=mixmeta.ml]{Maximum likelihood (ML) estimator}}
  \item \code{method="reml"}: \code{\link[=mixmeta.reml]{Restricted maximum likelihood (REML) estimator}}
  \item \code{method="mm"}: \code{\link[=mixmeta.mm]{Method of moments estimator}}
  \item \code{method="vc"}: \code{\link[=mixmeta.vc]{Variance components estimator}}
}

Note that non-standard random-effects models and the use of structured (co)variance matrices are only available for \code{"ml"} and \code{"reml"} methods. See their help pages for further details on the estimation procedures, following the links above.

Missing values are allowed in both sides of \code{formula}. In the case of missing predictors (right-hand side of \code{formula}), the related unit is entirely excluded from estimation. In contrast, a unit still contributes to estimation if at least outcome is non-missing. This behaviour is different from standard regression functions such as \code{\link{lm}} or \code{\link{glm}}. Before the call to \code{mixmeta.fit}, units matching such stricter missing definition are removed from the the model frame. The missing pattern in \code{S} must be consistent with that in \code{y}. See further details on handling \code{\link[=na.omit.data.frame.mixmeta]{missing values}} in \code{mixmeta}.

The fitting procedure can be controlled through the additional terms specified in \code{control}, which are passed to the function \code{\link{mixmeta.control}}.
}

\value{
The \code{mixmeta} function typically returns a list object of class \code{"mixmeta"} representing the meta-analytical model fit, as described in \code{\link{mixmetaObject}}. When \code{method="data.frame"}, the model is not fitted and the model frame is returned, namely a data frame with special attributes (see the default method \code{\link{model.frame}}) and, in this case, the additional class \code{"data.frame.mixmeta"}.

The wrapper function \code{mixmeta.fit} is usually called internally in \code{mixmeta}, and returns an intermediate list object with some of the components expected in the \code{"mixmeta"} class.

Several method functions for regression objects are available, either default or specifically written for the \code{"mixmeta"} class. See \code{\link{mixmetaObject}} for a complete list.
}

\references{
Sera F, Gasparrini A. (2019). An extended mixed-effects framework for meta-analysis.\emph{Statistics in Medicine}. In submission.
  
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/2012_gasparrini_statmed.html}{\bold{here}}].
}

\author{Antonio Gasparrini <\email{antonio.gasparrini@lshtm.ac.uk}> and Francesco Sera <\email{francesco.sera@lshtm.ac.uk}>}

\seealso{
See additional info on the estimation procedures at the related page of the fitting functions See \code{\link{mixmetaFormula}} for the use of formulae to define the fixed and random parts of the model. See alternative \code{\link[=mixmetaCovStruct]{(co)variance structures}} for likelihood-based estimation methods. See handling of \code{\link[=na.omit.data.frame.mixmeta]{missing values}} in \code{mixmeta}. See \code{\link[nlme]{lme}}, \code{\link{lm}} or \code{\link{glm}} for standard regression functions. See \code{\link{mixmeta-package}} for an overview of this modelling framework.
}

\examples{
### STANDARD MODELS

# RANDOM-EFFECTS META-ANALYSIS, ESTIMATED WITH REML
model <- mixmeta(logor, logorvar, data=bcg)
summary(model)

# RANDOM-EFFECTS META-REGRESSION, ESTIMATED WITH ML
model <- mixmeta(logor~ablat, logorvar, data=bcg, method="ml")
summary(model)


### MAIN METHOD FUNCTIONS

# COEFFICIENTS AND (CO)VARIANCE MATRIX
coef(model)
vcov(model)

# RESIDUALS AND FITTED VALUES
residuals(model)
fitted(model)

# MODEL FRAME AND MODEL MATRIX
model.frame(model)
model.matrix(model)

# LOG-LIKELIHOOD AND AIC VALUE
logLik(model)
AIC(model)

# COCHRAN Q TEST FOR RESIDUAL HETEROGENEITY
qtest(model)


### PREDICTIONS

# PREDICTED EFFECTS
predict(model)
predict(model, se=TRUE)
predict(model, newdata=data.frame(ablat=2:5*10), ci=TRUE)

# BEST LINEAR UNBIASED PREDICTION
blup(model)
blup(model, pi=TRUE)

# SEE help(predict.mixmeta) AND help(BLUP.mixmeta) FOR MORE INFO


### MULTIVARIATE MODELS

### BIVARIATE MODELS
model <- mixmeta(cbind(PD,AL) ~ pubyear, S=berkey98[5:7], data=berkey98)
summary(model)
residuals(model)

### MULTIVARIATE META-ANALYSIS WITH MORE THAN 2 OUTCOMES
y <- as.matrix(fibrinogen[2:5])
S <- as.matrix(fibrinogen[6:15])
model <- mixmeta(y, S)
summary(model)
predict(model, se=TRUE)
predict(model, se=TRUE, aggregate="outcome")


### OTHER EXTENSIONS

# MULTILEVEL META-ANALYSIS
model <- mixmeta(effect, var, random= ~ 1|district/study, data=school)
summary(model)
# SEE help(school) AND help(thrombolytic) FOR MORE EXAMPLES

# DOSE-RESPONSE META-ANALYSIS (SIMPLIFIED)
model <- mixmeta(logrr ~ 0 + dose, S=se^2, random= ~ 0 + dose|id, data=alcohol,
 subset=!is.na(se))
summary(model)
# SEE help(alcohol) FOR MORE EXAMPLES

# LONGITUDINAL META-ANALYSIS
model <- mixmeta(logOR~time, S=logORvar, random=~I(time-15)|study, data=gliomas)
summary(model)
# SEE help(gliomas) AND help(dbs) FOR MORE EXAMPLES


### FIXED-EFFECTS MODELS AND ALTERNATIVE ESTIMATORS

# FIXED-EFFECTS MODEL
model <- mixmeta(sbp~ish, S=sbp_se^2, data=hyp, method="fixed")
summary(model)

# METHOD OF MOMENTS
S <- as.matrix(hsls[5:10])
model <- mixmeta(cbind(b1,b2,b3), S, data=hsls, method="mm")
summary(model)

# VARIANCE COMPONENTS ESTIMATOR
model <- mixmeta(cbind(PD,AL)~pubyear, S=berkey98[5:7], data=berkey98,
  method="vc")
summary(model)


### IN THE PRESENCE OF MISSING VALUES

# RUN THE MODEL
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
model <- mixmeta(y, S)
summary(model)
model.frame(model)

# SEE help(na.omit.data.frame.mixmeta) FOR MORE EXAMPLES


### WHEN WITHIN-STUDY COVIARIANCES ARE NOT AVAILABLE AND/OR NEED TO BE INPUTTED

# GENERATE S
(S <- inputcov(hyp[c("sbp_se","dbp_se")], cor=hyp$rho))

# RUN THE MODEL
model <- mixmeta(cbind(sbp,dbp), S=S, data=hyp)

# INPUTTING THE CORRELATION DIRECTLY IN THE MODEL
model <- mixmeta(cbind(y1,y2), cbind(V1,V2), data=p53, control=list(Scor=0.95))
summary(model)

# SEE help(hyp) AND help(p53) FOR MORE EXAMPLES


### STRUCTURING THE BETWEEN-STUDY (CO)VARIANCE

# DIAGONAL
S <- as.matrix(hsls[5:10])
model <- mixmeta(cbind(b1,b2,b3), S, data=hsls, bscov="diag")
summary(model)
model$Psi

# COMPOUND SYMMETRY
model <- mixmeta(cbind(b1,b2,b3), S, data=hsls, bscov="cs")
summary(model)
model$Psi

# SEE help(mixmetaCovStruct) FOR DETAILS AND ADDITIONAL EXAMPLES


### USE OF THE CONTROL LIST

# PRINT THE ITERATIONS AND CHANGE THE DEFAULT FOR STARTING VALUES
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
model <- mixmeta(y, S, control=list(showiter=TRUE, igls.inititer=20))

# SEE help(mixmeta.control) FOR FURTHER DETAILS
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}

