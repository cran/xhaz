% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qbs.R
\name{qbs}
\alias{qbs}
\title{qbs function}
\usage{
qbs(x)
}
\arguments{
\item{x}{a covariate to be considered in the \code{xhaz} formula with a
time-dependant effect. Quadratic B-splines with two interior knots are used.}
}
\value{
No return value, called for side effects.
}
\description{
a function indicating which covariates have a time-dependent
effect in the formula.
}
\examples{

\donttest{

library("xhaz")
library("numDeriv")
library("survexp.fr")
library("splines")

fit.tdphBS <- xhaz(
              formula = Surv(obs_time_year, event) ~ ageCentre + qbs(immuno_trt),
              data = dataCancer, ratetable = survexp.fr,
              interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
              rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
              baseline = "bsplines", pophaz = "classic")

print(fit.tdphBS)
}

}
\keyword{qbs}
