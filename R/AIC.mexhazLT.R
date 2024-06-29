#' @title Akaike's Information Criterion for excess hazard model from mexhazLT function
#'
#' @description Calculates the Akaike's Information Criterion' for fitted models from `mexhazLT`.
#'
#' @param object a fitted model object obtained from `mexhazLT` function
#'
#' @param ... optionally more fitted model objects obtained from `mexhazLT` function
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
#' \donttest{
#' library("xhaz")
#'
#' data("breast")
#' # load the data sets 'breast'.
#'
#'  # Flexible mexhaz model: baseline excess hazard with cubic B-splines
#'  # assumption on the life table available :
#'  # other cause mortality in the cohort is comparable to the mortality
#'  # observed in the general population with the same characteristics.
#'
#' # The life table to be used is survexp.us. Note that SEX is coded 2 instead of female in survexp.us.
#' breast$sexe <- "female"
#'
#' fit.haz <- exphaz(
#'                   formula = Surv(temps, statut) ~ 1,
#'                   data = breast, ratetable = survexp.us,
#'                   only_ehazard = FALSE,
#'                   rmap = list(age = 'age', sex = 'sexe', year = 'date'))
#'
#' breast$expected <- fit.haz$ehazard
#' breast$expectedCum <- fit.haz$ehazardInt
#'
#' mod.bs <- mexhazLT(formula = Surv(temps, statut) ~ agecr + armt,
#'                   data = breast,
#'                   ratetable = survexp.us, degree = 3,
#'                   knots=quantile(breast[breast$statut==1,]$temps, probs=c(1:2/3)),
#'                   expected = "expected",expectedCum = "expectedCum",
#'                   base = "exp.bs", pophaz = "classic")
#'
#' mod.bs
#'
#' AIC(mod.bs)
#' }
#' @keywords internal
#' @export
AIC.mexhazLT <- function(object, ..., k = 2) {
  dots.object <- list(...)
  if (length(dots.object) == 0) {
    if (inherits(object, "mexhazLT")) {
      df <- object$n.par
      val <- (k * object$n.par - 2 * (object$loglik))
    } else{
      stop("object must be a mexhazLT function output")
    }
    return(val)
  } else{
    object <- list(object, ...)
    aic_bis <- function(i) {
      if (inherits(object[[i]], "mexhazLT")) {
        df <- object[[i]]$n.par

        val <- (k * object[[i]]$n.par - 2 * (object[[i]]$loglik))
        resval <- data.frame(df = df, AIC = val)
        return(resval)
      } else{
        stop("object must be a mexhazLT function output")
      }
    }

    val <- sapply(1:length(object), aic_bis)
    Call <- match.call()
    colnames(val) <- as.character(Call[-1])
    return(val)
  }
}






