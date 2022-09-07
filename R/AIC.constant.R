#' @title Akaike's An Information Criterion for excess hazard model with
#' baseline hazard following a piecewise constant function
#'
#' @description Calculates the Akaike's ‘An Information Criterion’ for fitted
#' models from `xhaz`.
#'
#' @param object a fitted model object obtained from `xhaz` function
#'
#' @param ... optionally more fitted model objects obtained from `xhaz` function
#'
#' @param k numeric, the penalty per parameter to be used; the default \code{k = 2}
#' is the classical AIC.
#'
#' @return the value corresponds to the AIC calculated from the total
#' log-likelihood of the fitted model if just one object is provided.
#' If multiple objects are provided, a data.frame with columns corresponding to the
#' objects and rows representing the number of parameters in the model (df) and the AIC
#'
#' @examples
#' library("xhaz")
#'
#'# Esteve et al. model: baseline excess hazard is a piecewise function
#'#                      linear and proportional effects for the covariates on
#'#                      baseline excess hazard.
#'
#' levels(simuData$sex) <- c("male", "female")
#'
#' set.seed(1980)
#' simuData2 <- simuData[sample(nrow(simuData), size = 500), ]
#'
#' fit.estv2 <- xhaz(formula = Surv(time_year, status) ~ agec + race,
#'                   data = simuData2,
#'                   ratetable = survexp.us,
#'                   interval = c(0, NA, NA, NA, NA, NA, 6),
#'                   rmap = list(age = 'age', sex = 'sex', year = 'date'),
#'                   baseline = "constant", pophaz = "classic")
#'
#'
#' fit.estv2
#'
#' AIC(fit.estv2)
#'
#'
#' @export
AIC.constant <- function(object, ..., k = 2) {
  dots.object <- list(...)
  if (length(dots.object) == 0) {
    if (inherits(object, "constant")) {
      df <- length(object$coefficients)
      val <-
        (k * length(object$coefficients) - 2 * (object$loglik))[2]
    } else{
      stop("object must be a xhaz function output")
    }
    return(val)
  } else{
    object <- list(object, ...)
    aic_bis <- function(i) {
      if (inherits(object[[i]], "constant")) {
        df <- length(object[[i]]$coefficients)

        val <-
          (k * length(object[[i]]$coefficients) - 2 * (object[[i]]$loglik))[2]
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






