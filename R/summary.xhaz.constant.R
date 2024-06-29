#' @title A summary.constant Function used to print a object of class \code{xhaz.constant}
#'
#' @description This function present the estimated coefficients for the excess
#' hazard baseline coefficient and for the covariate effects
#'
#' @param object an object of class xhaz.constant
#'
#' @param ci_type method for confidence intervals calculation
#'
#'
#' @param ... additionnal parameters which can be used in the \code{print}
#' function
#'
#'
#' @return Estimated parameters of the model in different scales for interpretation purposes.
#'
#'
#' @keywords summary.constant
#'
#' @seealso \code{\link{xhaz}}, \code{\link{summary.constant}}, \code{\link{print.bsplines}}
#'
#' @examples
#'
#' library("xhaz")
#' library("numDeriv")
#' data("simuData", package = "xhaz")    # load the data sets 'simuData'
#'
#' # Esteve et al. model: baseline excess hazard is a piecewise function
#' #                      linear and proportional effects for the covariates on
#' #                      baseline excess hazard.
#'
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
#'
#'
#'
#' summary(fit.estv2)
#'
#' @importFrom stats printCoefmat
#' @export



summary.constant <-
  function(object,
           ci_type = "lognormal",
           ...)
  {
    digits <- max(options()$digits - 4, 3)

    if (is.null(object$coefficients)) {
      return(object)
    }
    cl <- try(object$call)
    if (!is.null(cl)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }
    if (!is.null(object$fail)) {
      cat(" Esteveph failed.", object$fail, "\n")
      return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    if (!is.null(object$add.rmap)) {
      nalpha <- nlevels(object$add.rmap)
      if (nalpha == 1) {
        indxAlpha <- which(stringr::str_detect(names(object$coefficients),
                                               pattern = "alpha"))
        names(object$coef)[indxAlpha] <- "log(alpha)"
      }
      else{
        indxAlpha <- which(stringr::str_detect(names(object$coefficients),
                                               pattern = "alpha"))
        if (object$add.rmap.cut$breakpoint == TRUE) {
          initnames <- names(object$coef)[c(indxAlpha)]
          names(object$coef)[c(indxAlpha)] <- paste("log(", initnames,
                                                    ")",
                                                    sep = "")
        } else {
          names(object$coef)[c(indxAlpha)] <- paste("log(", paste0('alpha.',
                                                                   levels(object$add.rmap)),
                                                    ")",
                                                    sep = "")
        }

        nalpha <- length(indxAlpha)

      }
    } else{
      nalpha <- 0
    }

    nstrata <- ifelse(is.null(attr(object$terms, "nstrata")),
                      1,
                      attr(object$terms, "nstrata"))
    nvar <-
      length(object$coef) - nstrata * (length(object$interval) - 1) - nalpha
    coef <- object$coef
    object$var <- object$varcov
    se <- numeric(length(coef))
    if (nvar > 0) {
      if (nvar != 1) {
        se[1:nvar] <- sqrt(diag(object$var[1:nvar, 1:nvar]))
      }
      else{
        se[1:nvar] <- sqrt(object$var[1])
      }
    }

    if (ci_type == "delta.method") {
      if (!is.null(object$add.rmap)) {
        se[(nvar + 1):(length(object$coef) - length(indxAlpha))] <- c(sqrt(exp(object$coef[(nvar + 1):(length(object$coef) - length(indxAlpha))]) %*%
                                                                             object$var[(nvar + 1):(length(object$coef) - length(indxAlpha)),
                                                                                   (nvar + 1):(length(object$coef) - length(indxAlpha))] %*%
                                                                             exp(object$coef[(nvar + 1):(length(object$coef) - length(indxAlpha))])))

        coef[(nvar + 1):(length(object$coef) - length(indxAlpha))] <-
          exp(coef[(nvar + 1):(length(object$coef) - length(indxAlpha))])
        if (length(indxAlpha) > 1) {
          se[c(indxAlpha)] <-
            sqrt(diag(object$var[c(indxAlpha), c(indxAlpha)]))

        } else{
          se[c(indxAlpha)] <- (sqrt(object$var[c(indxAlpha), c(indxAlpha)]))
        }
      } else{
        se[(nvar + 1):length(object$coef)] <- sqrt(exp(object$coef[(nvar + 1):length(object$coef)]) %*%
                                                     object$var[(nvar + 1):length(object$coef),
                                                           (nvar + 1):length(object$coef)] %*%
                                                     exp(object$coef[(nvar + 1):length(object$coef)]))
        coef[(nvar + 1):length(object$coef)] <-
          exp(coef[(nvar + 1):length(object$coef)])
      }

    } else if (ci_type == "lognormal") {
      #lognormal based CI
      if (!is.null(object$add.rmap)) {
        se[(nvar + 1):(length(object$coef) - length(indxAlpha))] <-
          c(sqrt(diag((
            exp(2 * object$coef[(nvar + 1):(length(object$coef) - length(indxAlpha))] +
                  object$var[(nvar + 1):(length(object$coef) - length(indxAlpha)),
                        (nvar + 1):(length(object$coef) - length(indxAlpha))])
          ) *
            (
              exp(object$var[(nvar + 1):(length(object$coef) - length(indxAlpha)),
                             (nvar + 1):(length(object$coef) - length(indxAlpha))]) - 1
            ))))

        coef[(nvar + 1):(length(object$coef) - length(indxAlpha))] <-
          exp(coef[(nvar + 1):(length(object$coef) - length(indxAlpha))] + 1 / 2 *
                diag(object$var[(nvar + 1):(length(object$coef) - length(indxAlpha)), (nvar + 1):(length(object$coef) - length(indxAlpha))]))
        if (length(indxAlpha) > 1) {
          se[c(indxAlpha)] <-
            (sqrt(diag(object$var[c(indxAlpha), c(indxAlpha)])))

        } else{
          se[c(indxAlpha)] <- (sqrt(object$var[c(indxAlpha), c(indxAlpha)]))
        }
      } else{
        se[(nvar + 1):length(object$coef)] <- sqrt(diag(exp(2 * object$coef[(nvar + 1):length(object$coef)] + (object$var[(nvar + 1):length(object$coef),
                                                                                                                     (nvar + 1):length(object$coef)])) * (exp(object$var[(nvar + 1):length(object$coef),
                                                                                                                                                                         (nvar + 1):length(object$coef)]) - 1)))
        coef[(nvar + 1):length(object$coef)] <-
          exp(coef[(nvar + 1):length(object$coef)] + 1 / 2 * diag(object$var[(nvar + 1):length(object$coef), (nvar + 1):length(object$coef)]))
      }

    }



    if (is.null(coef) | is.null(se))
      stop("Input is not valid")
    tmp <-
      cbind(
        coef,
        se,
        coef - abs(qnorm((1 - object$level) / 2)) * se,
        coef + abs(qnorm((1 - object$level) / 2)) * se,
        coef / se,
        signif(1 - pchisq((coef / se) ^ 2, 1), digits - 1)
      )
    dimnames(tmp) <-
      list(names(coef),
           c(
             "coef",
             "se(coef)",
             paste("lower", object$level, sep = " "),
             paste("upper", object$level, sep = " "),
             "z",
             "Pr(>|z|)"
           ))
    cat("\n")
    printCoefmat(
      tmp,
      P.values = TRUE,
      digits = digits,
      signif.stars = TRUE,
      na.print = "NA",
      ...
    )

    if (nvar > 0) {
      if (!is.null(object$add.rmap)) {
        index <- c(1:nvar, indxAlpha)
        nalpha <- nlevels(object$add.rmap)
        if (nalpha == 1) {
          names(object$coef)[indxAlpha] <- "alpha"
        }
        else{
          if (object$add.rmap.cut$breakpoint == FALSE) {
            names(object$coef)[c(indxAlpha)] <-
              paste0('alpha.', levels(object$add.rmap))

          } else {
            names(object$coef)[c(indxAlpha)] <- initnames

          }
        }
      }
      coef <- object$coef



      if (!is.null(object$add.rmap)) {
        if (ci_type == "delta.method") {
          se_alpha <- sapply(1:length(indxAlpha), function(i)
            (matrix(exp(
              object$coef[c(indxAlpha[i])]
            )) %*% sqrt(object$var[indxAlpha[i], indxAlpha[i]])))
          coef_alpha <- exp(coef[indxAlpha])

          mlevel <- abs(qnorm((1 - object$level) / 2))
          if (inherits(attributes(object$terms)$factors, "matrix")) {
            tmp_new <- cbind((exp(coef[index])),
                             c(exp(coef[1:nvar] - abs(qnorm((1 - object$level) / 2
                             )) * se[1:nvar]),
                             c(coef_alpha - mlevel * c(se_alpha))),
                             c(exp(coef[1:nvar] + abs(qnorm((1 - object$level) / 2
                             )) * se[1:nvar]),
                             coef_alpha + mlevel * c(se_alpha)))

            }

        } else {
          #lognormal based
          if (length(indxAlpha) > 1) {
            se_alpha <- sqrt(diag(object$var[indxAlpha, indxAlpha]))

          } else{
            se_alpha <- sqrt((object$var[indxAlpha, indxAlpha]))

          }

          coef_alpha <- exp(coef[indxAlpha])

          mlevel <- abs(qnorm((1 - object$level) / 2))
          if (inherits(attributes(object$terms)$factors, "matrix")) {
            tmp_new <- cbind(c(exp(coef[1:nvar]), coef_alpha),
                             c(exp(coef[1:nvar] - abs(qnorm((1 - object$level) / 2
                             )) *
                               se[1:nvar]),
                             exp(c(
                               coef[indxAlpha] - mlevel * c(se_alpha)
                             ))),
                             c(exp(coef[1:nvar] + abs(qnorm((1 - object$level) / 2
                             )) *
                               se[1:nvar]),
                             exp(coef[indxAlpha] + mlevel * c(se_alpha))))
          }

        }


      } else{

        if (inherits(attributes(object$terms)$factors, "matrix")) {
          index <- c(1:nvar)

          tmp_new <- cbind(exp(coef[1:nvar]),
                           exp(coef[1:nvar] - abs(qnorm((1 - object$level) / 2
                           )) * se[1:nvar]),
                           exp(coef[1:nvar] + abs(qnorm((1 - object$level) / 2
                           )) * se[1:nvar]))
        }

      }

      if (inherits(attributes(object$terms)$factors, "matrix")) {

        dimnames(tmp_new) <- list(names(coef[index]),
                                  c(
                                    "exp(coef)",
                                    paste("lower", object$level, sep = " "),
                                    paste("upper", object$level, sep = " ")
                                  ))
      }




      cat("\n")
      if (object$pophaz != "classic") {
        if (object$pophaz == "rescaled") {
          cat("\n")
          cat(
            "Excess hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
          )
          if (inherits(attributes(object$terms)$factors, "matrix")) {
          lines_al <- which(stringr::str_detect(rownames(tmp_new),
                                                pattern = "alpha"))
          print(tmp_new[-c(lines_al), ], digits = digits + 1)
          cat("\n")
          cat("and rescaled parameter on population hazard \n")

          tmp_new_alpha <- matrix(tmp_new[c(lines_al), ], nrow = 1)
          dimnames(tmp_new_alpha)[[1]] <- "alpha"
          dimnames(tmp_new_alpha)[[2]] <- dimnames(tmp_new)[[2]]
          print(tmp_new_alpha, digits = digits + 1)
          }

        } else if (object$pophaz == "corrected" &
                   object$add.rmap.cut$breakpoint == FALSE) {
          cat("\n")
          cat(
            "Excess hazard hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
          )
          if (inherits(attributes(object$terms)$factors, "matrix")) {
          lines_al <- which(stringr::str_detect(rownames(tmp_new),
                                                pattern = "alpha"))
          print(tmp_new[-c(lines_al), ], digits = digits + 1)
          cat("\n")
          cat("and corrected scale parameters on population hazard \n")
          print(tmp_new[c(lines_al), ], digits = digits + 1)
          }
        } else if (object$pophaz == "corrected" &
                   object$add.rmap.cut$breakpoint == TRUE) {
          cat("\n")
          cat(
            "Excess hazard hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
          )
          if (inherits(attributes(object$terms)$factors, "matrix")) {
          lines_al <- which(stringr::str_detect(rownames(tmp_new),
                                                pattern = "alpha"))
          print(tmp_new[-c(lines_al), ], digits = digits + 1)
          cat("\n")
          cat(
            "and corrected scale parameters on population hazard \n (non proportional correction using breakpoint approach)\n"
          )
          cat("\n")

          print(tmp_new[c(lines_al), ], digits = digits + 1)
          }


          if (!is.na(object$add.rmap.cut$cut[1])) {
            n_break <- length(object$add.rmap.cut$cut)
            n_int <- 1 + n_break

          }

          cat("\n")

        }

      } else{
        cat("\n")
        cat(
          "Excess hazard hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
        )
        if (inherits(attributes(object$terms)$factors, "matrix")) {
        print(tmp_new, digits = digits + 1)
        }
      }


      logtest <- (-2 * (object$loglik[1] - object$loglik[2]))
      df <- length(object$coef)
      cat("\n")
      cat("number of observations:",
          paste0(format(object$n), "; "),
          "number of events:",
          object$n.events)
      cat("\n")
      cat(
        "Likelihood ratio test: ",
        format(round(logtest, 2)),
        "  on ",
        df,
        " degree(s) of freedom,",
        " p=",
        format(1 - pchisq(logtest, df)),
        sep = ""
      )
      cat("\n")
    }

    if (sum(object$cov.test) > 0) {
      cat("\n")
      cat("Results of tests for '",
          names(coef)[grep("T", as.character(object$cov.test))], "' equal to 0")
      cat("\n")
      cat(
        "  Likelihood ratio test=",
        format(round(object$loglik.test, 2)),
        "  on ",
        object$cov.df,
        " degree(s) of freedom,",
        " p=",
        format(1 - pchisq(object$loglik.test, object$cov.df)),
        sep = ""
      )
      cat("\n")
      cat(
        "  Wald test=",
        format(round(object$wald.test, 2)),
        "  on ",
        object$cov.df,
        " degree(s) of freedom,",
        " p=",
        format(1 - pchisq(object$wald.test, object$cov.df)),
        sep = ""
      )
      cat("\n")
      cat(
        "  Score test=",
        format(round(object$score.test, 2)),
        "  on ",
        object$cov.df,
        " degree(s) of freedom,",
        " p=",
        format(1 - pchisq(object$score.test, object$cov.df)),
        sep = ""
      )
      cat("\n")
      cat("number of observations:",
          format(object$n),
          "number of events:",
          object$n.events)
      cat("\n")



    }
    if (inherits(attributes(object$terms)$factors, "matrix")) {
      if (any(tmp_new[, "lower 0.95"] < 0, na.rm = TRUE)) {
        warning("\nlower 0.95 CI approximation may be incorrect")
      }
      }

    invisible()
  }



