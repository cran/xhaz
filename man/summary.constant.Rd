% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.xhaz.constant.R
\name{summary.constant}
\alias{summary.constant}
\title{A summary.constant Function used to print a object of class \code{xhaz.constant}}
\usage{
\method{summary}{constant}(object, ci_type = "lognormal", ...)
}
\arguments{
\item{object}{an object of class xhaz.constant}

\item{ci_type}{method for confidence intervals calculation}

\item{...}{additionnal parameters which can be used in the \code{print}
function}
}
\value{
Estimated parameters of the model in different scales for interpretation purposes.
}
\description{
This function present the estimated coefficients for the excess
hazard baseline coefficient and for the covariate effects
}
\examples{

library("xhaz")
library("numDeriv")
data("simuData", package = "xhaz")    # load the data sets 'simuData'

# Esteve et al. model: baseline excess hazard is a piecewise function
#                      linear and proportional effects for the covariates on
#                      baseline excess hazard.


set.seed(1980)
simuData2 <- simuData[sample(nrow(simuData), size = 500), ]

fit.estv2 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
                  data = simuData2,
                  ratetable = survexp.us,
                  interval = c(0, NA, NA, NA, NA, NA, 6),
                  rmap = list(age = 'age', sex = 'sex', year = 'date'),
                  baseline = "constant", pophaz = "classic")





summary(fit.estv2)

}
\seealso{
\code{\link{xhaz}}, \code{\link{summary.constant}}, \code{\link{print.bsplines}}
}
\keyword{summary.constant}
