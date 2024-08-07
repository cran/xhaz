% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIC.mexhazLT.R
\name{BIC.mexhazLT}
\alias{BIC.mexhazLT}
\title{Bayesian Information Criterion for excess hazard model from mexhazLT function}
\usage{
\method{BIC}{mexhazLT}(object, ...)
}
\arguments{
\item{object}{a fitted model object obtained from \code{mexhazLT} function}

\item{...}{optionally more fitted model objects obtained from \code{mexhazLT} function}
}
\value{
the value corresponds to the BIC calculated from the total
log-likelihood of the fitted model if just one object is provided.
If multiple objects are provided, a data.frame with columns corresponding to the
objects and rows representing the number of parameters in the model (df) and the BIC.
}
\description{
Calculates the Bayesian Information Criterion' for fitted
models from \code{mexhazLT}.
}
\examples{
\donttest{
library("xhaz")

data("breast")
# load the data sets 'breast'.

 # Flexible mexhaz model: baseline excess hazard with cubic B-splines
 # assumption on the life table available :
 # other cause mortality in the cohort is comparable to the mortality
 # observed in the general population with the same characteristics.

# The life table to be used is survexp.us. Note that SEX is coded 2 instead of female in survexp.us.
breast$sexe <- "female"

fit.haz <- exphaz(
                  formula = Surv(temps, statut) ~ 1,
                  data = breast, ratetable = survexp.us,
                  only_ehazard = FALSE,
                  rmap = list(age = 'age', sex = 'sexe', year = 'date'))

breast$expected <- fit.haz$ehazard
breast$expectedCum <- fit.haz$ehazardInt

mod.bs <- mexhazLT(formula = Surv(temps, statut) ~ agecr + armt,
                  data = breast,
                  ratetable = survexp.us, degree = 3,
                  knots=quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3)),
                  expected = "expected",expectedCum = "expectedCum",
                  base = "exp.bs", pophaz = "classic")

mod.bs

BIC(mod.bs)
}

}
\keyword{internal}
