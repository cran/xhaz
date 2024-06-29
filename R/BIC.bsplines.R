#' @title Bayesian Information Criterion for excess hazard model with
#' baseline hazard following a B-splines functions
#'
#' @description Calculates the Bayesian Information Criterion' for fitted
#' models from `xhaz`.
#'
#' @param object a fitted model object obtained from `xhaz` function
#'
#' @param ... optionally more fitted model objects obtained from `xhaz` function
#'
#' @return the value corresponds to the BIC calculated from the total
#' log-likelihood of the fitted model if just one object is provided.
#' If multiple objects are provided, a data.frame with columns corresponding to the
#' objects and rows representing the number of parameters in the model (df) and the BIC.
#'
#' @examples
#' \donttest{
#' library("xhaz")
#'
#' #Giorgi et al model: baseline excess hazard is a quadratic Bsplines
#' #                    function with two interior knots and allow here a
#' #                    linear and proportional effects for the covariates on
#' #                    baseline excess hazard.
#'
#' fitphBS <- xhaz(formula = Surv(time_year, status) ~ agec + race,
#'                 data = simuData,
#'                 ratetable = survexp.us,
#'                 interval = c(0, NA, NA, 6),
#'                 rmap = list(age = 'age', sex = 'sex', year = 'date'),
#'                 baseline = "bsplines", pophaz = "classic")
#'
#' fitphBS
#' BIC(fitphBS)
#' }
#'
#' @keywords internal
#' @export
BIC.bsplines <- function(object, ...) {

  dots.object <- list(...)
  if (length(dots.object) == 0) {
    if (inherits(object, "bsplines")) {
      df <- length(object$coefficients)
      val <- (log((object$n)) * length(object$coefficients) - 2 * (object$loglik))[2]
    } else{
      stop("object must be a xhaz function output")
    }
    return(val)
  } else{
    object <- list(object, ...)
    aic_bis <- function(i) {
      if (inherits(object[[i]], "bsplines")) {
        df <- length(object[[i]]$coefficients)

        val <-
          (log((object[[i]]$n)) * length(object[[i]]$coefficients) - 2 * (object[[i]]$loglik))[2]
        resval <- data.frame(df = df, AIC = val)
        return(resval)
      } else{
        stop("object must be a xhaz function output")
      }
    }

    val <- sapply(1:length(object), aic_bis)
    Call <- match.call()
    colnames(val) <- as.character(Call[-1])
    return(val)
  }

}
