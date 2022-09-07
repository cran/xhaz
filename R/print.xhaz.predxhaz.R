#' @title A print.predxhaz Function used to print a object of class predxhaz
#'
#' @description  This function present the print of the predict function
#'
#' @param x an object of class predxhaz
#'
#' @param ... other parameters used for print function
#'
#' @return an object of class data.frame containing the following components:
#'
#'
#' \item{times.pts}{The time at which the estimations of excess hazard and net
#' survival are predicted}
#'
#' \item{hazard}{the predicted excess hazard at the fixed times}
#'
#' \item{survival}{the predicted net survival at the fixed times}
#'
#' @keywords print.predxhaz
#'
#' @examples
#'
#' \donttest{
#'
#' library("xhaz")
#' library("survexp.fr")
#' library("splines")
#'
#' data("dataCancer", package = "xhaz")   # load the data set in the package
#'
#' fit.phBS <- xhaz(
#'         formula = Surv(obs_time_year, event) ~ ageCentre + immuno_trt,
#'         data = dataCancer, ratetable = survexp.fr,
#'         interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
#'         rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
#'         baseline = "bsplines", pophaz = "classic")
#'
#'
#' fit.phBS
#'
#'
#' predicted <- predict(object = fit.phBS,
#'                      new.data = dataCancer[1:10,],
#'                      times.pts = c(seq(0,10,1)),
#'                      baseline = TRUE)
#'
#'
#'
#' #a list of predicted hazard and survival at different time points
#' print(predicted)
#'
#'
#' #predicted hazard and survival at time points 10 years
#' print(predicted[[10]])
#' }
#'
#' @export

print.predxhaz <- function(x, ...)
{
  if (any(class(x) == "predxhaz")) {
    cl <- try(attributes(x)$call)
    if (!is.null(cl)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }
    attributes(x) <- NULL
    print(x, ...)
    invisible()
  }

}

