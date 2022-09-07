#' @title A print.constant Function used to print a object of class constant
#'
#' @description This function present the estimated coefficients for the excess
#' hazard baseline coefficient and for the covariate effects
#'
#' @param x an object of class xhaz.constant
#'
#' @param digits minimal number of significant digits.
#'
#' @param ci_type method for confidence intervals calculation
#'
#' @param ... additionnal parameters which can be used in the \code{print}
#' function
#'
#' @return Estimated parameters of the model in different scales for interpretation purposes.
#'
#'
#' @keywords print.constant
#'
#' @seealso \code{\link{xhaz}}, \code{\link{summary.constant}}, \code{\link{print.bsplines}}
#'
#' @examples
#'
#' library("numDeriv")
#' library("survexp.fr")
#'
#' data("simuData","rescaledData", "dataCancer")
#' # load the data sets 'simuData', 'rescaledData' and 'dataCancer'.
#'
#' # Esteve et al. model: baseline excess hazard is a piecewise function
#' #                      linear and proportional effects for the covariates on
#' #                      baseline excess hazard.
#'
#' levels(simuData$sex) <- c("male", "female")
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
#' print(fit.estv2)
#'
#'
#' @importFrom stats printCoefmat
#'
#'
#' @export



print.constant <-
  function(x,
           ci_type = "lognormal",
           digits = max(options()$digits - 4, 3),
           ...)
  {
    if (is.null(x$coefficients)) {
      return(x)
    }
    cl <- try(x$call)
    if (!is.null(cl)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }
    if (!is.null(x$fail)) {
      cat(" Esteveph failed.", x$fail, "\n")
      return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    if (!is.null(x$add.rmap)) {
      nalpha <- nlevels(x$add.rmap)
      if (nalpha == 1) {
        indxAlpha <- which(stringr::str_detect(names(x$coefficients),
                                               pattern = "alpha"))
        names(x$coef)[indxAlpha] <- "log(alpha)"
      }
      else{
        indxAlpha <- which(stringr::str_detect(names(x$coefficients),
                                               pattern = "alpha"))
        if (x$add.rmap.cut$breakpoint == TRUE) {
          initnames <- names(x$coef)[c(indxAlpha)]
          names(x$coef)[c(indxAlpha)] <- paste("log(", initnames,
                                               ")",
                                               sep = "")
        } else {
          names(x$coef)[c(indxAlpha)] <- paste("log(", paste0('alpha.',
                                                              levels(x$add.rmap)),
                                               ")",
                                               sep = "")
        }

        nalpha <- length(indxAlpha)

      }
    } else{
      nalpha <- 0
    }

    nstrata <- ifelse(is.null(attr(x$terms, "nstrata")),
                      1,
                      attr(x$terms, "nstrata"))
    nvar <-
      length(x$coef) - nstrata * (length(x$interval) - 1) - nalpha
    coef <- x$coef
    x$var <- x$varcov
    se <- numeric(length(coef))
    if (nvar > 0) {
      if (nvar != 1) {
        se[1:nvar] <- sqrt(diag(x$var[1:nvar, 1:nvar]))
      }
      else{
        se[1:nvar] <- sqrt(x$var[1])
      }
    }

    if (ci_type == "delta.method") {
      if (!is.null(x$add.rmap)) {
        se[(nvar + 1):(length(x$coef) - length(indxAlpha))] <- c(sqrt(exp(x$coef[(nvar + 1):(length(x$coef) - length(indxAlpha))]) %*%
                                                                        x$var[(nvar + 1):(length(x$coef) - length(indxAlpha)),
                                                                              (nvar + 1):(length(x$coef) - length(indxAlpha))] %*%
                                                                        exp(x$coef[(nvar + 1):(length(x$coef) - length(indxAlpha))])))

        coef[(nvar + 1):(length(x$coef) - length(indxAlpha))] <-
          exp(coef[(nvar + 1):(length(x$coef) - length(indxAlpha))])
        if (length(indxAlpha) > 1) {
          se[c(indxAlpha)] <- sqrt(diag(x$var[c(indxAlpha), c(indxAlpha)]))

        } else{
          se[c(indxAlpha)] <- (sqrt(x$var[c(indxAlpha), c(indxAlpha)]))
        }
      } else{
        se[(nvar + 1):length(x$coef)] <- sqrt(exp(x$coef[(nvar + 1):length(x$coef)]) %*%
                                                x$var[(nvar + 1):length(x$coef),
                                                      (nvar + 1):length(x$coef)] %*%
                                                exp(x$coef[(nvar + 1):length(x$coef)]))
        coef[(nvar + 1):length(x$coef)] <-
          exp(coef[(nvar + 1):length(x$coef)])
      }

    } else if (ci_type == "lognormal") {
      #lognormal based CI
      if (!is.null(x$add.rmap)) {
        se[(nvar + 1):(length(x$coef) - length(indxAlpha))] <- c(sqrt(diag((
          exp(2 * x$coef[(nvar + 1):(length(x$coef) - length(indxAlpha))] +
                x$var[(nvar + 1):(length(x$coef) - length(indxAlpha)),
                      (nvar + 1):(length(x$coef) - length(indxAlpha))])
        ) *
          (
            exp(x$var[(nvar + 1):(length(x$coef) - length(indxAlpha)),
                      (nvar + 1):(length(x$coef) - length(indxAlpha))]) - 1
          ))))

        coef[(nvar + 1):(length(x$coef) - length(indxAlpha))] <-
          exp(coef[(nvar + 1):(length(x$coef) - length(indxAlpha))] + 1 / 2 * diag(x$var[(nvar + 1):(length(x$coef) - length(indxAlpha)), (nvar + 1):(length(x$coef) - length(indxAlpha))]))
        if (length(indxAlpha) > 1) {
          se[c(indxAlpha)] <- (sqrt(diag(x$var[c(indxAlpha), c(indxAlpha)])))

        } else{
          se[c(indxAlpha)] <- (sqrt(x$var[c(indxAlpha), c(indxAlpha)]))
        }
      } else{
        se[(nvar + 1):length(x$coef)] <- sqrt(diag(exp(2 * x$coef[(nvar + 1):length(x$coef)] + (x$var[(nvar + 1):length(x$coef),
                                                                                                      (nvar + 1):length(x$coef)])) * (exp(x$var[(nvar + 1):length(x$coef),
                                                                                                                                                (nvar + 1):length(x$coef)]) - 1)))
        coef[(nvar + 1):length(x$coef)] <-
          exp(coef[(nvar + 1):length(x$coef)] + 1 / 2 * diag(x$var[(nvar + 1):length(x$coef), (nvar + 1):length(x$coef)]))
      }

    }



    if (is.null(coef) | is.null(se))
      stop("Input is not valid")
    tmp <-
      cbind(
        coef,
        se,
        coef - abs(qnorm((1 - x$level) / 2)) * se,
        coef + abs(qnorm((1 - x$level) / 2)) * se,
        coef / se,
        signif(1 - pchisq((coef / se) ^ 2, 1), digits - 1)
      )
    dimnames(tmp) <-
      list(names(coef),
           c(
             "coef",
             "se(coef)",
             paste("lower", x$level, sep = " "),
             paste("upper", x$level, sep = " "),
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
      if (!is.null(x$add.rmap)) {
        index <- c(1:nvar, indxAlpha)
        nalpha <- nlevels(x$add.rmap)
        if (nalpha == 1) {
          names(x$coef)[indxAlpha] <- "alpha"
        }
        else{
          if (x$add.rmap.cut$breakpoint == FALSE) {
            names(x$coef)[c(indxAlpha)] <- paste0('alpha.', levels(x$add.rmap))

          } else {
            names(x$coef)[c(indxAlpha)] <- initnames

          }
        }
      }
      coef <- x$coef



      if (!is.null(x$add.rmap)) {
        if (ci_type == "delta.method") {
          se_alpha <- sapply(1:length(indxAlpha), function(i)
            (matrix(exp(x$coef[c(indxAlpha[i])])) %*% sqrt(x$var[indxAlpha[i], indxAlpha[i]])))
          coef_alpha <- exp(coef[indxAlpha])

          mlevel <- abs(qnorm((1 - x$level) / 2))

          tmp_new <- cbind((exp(coef[index])),
                           c(exp(coef[1:nvar] - abs(qnorm((1 - x$level) / 2
                           )) * se[1:nvar]),
                           c(coef_alpha - mlevel * c(se_alpha))),
                           c(exp(coef[1:nvar] + abs(qnorm((1 - x$level) / 2
                           )) * se[1:nvar]),
                           coef_alpha + mlevel * c(se_alpha)))
        } else {
          #lognormal based
          if (length(indxAlpha) > 1) {
            se_alpha <- sqrt(diag(x$var[indxAlpha, indxAlpha]))

          } else{
            se_alpha <- sqrt((x$var[indxAlpha, indxAlpha]))

          }

          coef_alpha <- exp(coef[indxAlpha])

          mlevel <- abs(qnorm((1 - x$level) / 2))

          tmp_new <- cbind(c(exp(coef[1:nvar]), coef_alpha),
                           c(exp(coef[1:nvar] - abs(qnorm((1 - x$level) / 2
                           )) *
                             se[1:nvar]),
                           exp(c(
                             coef[indxAlpha] - mlevel * c(se_alpha)
                           ))),
                           c(exp(coef[1:nvar] + abs(qnorm((1 - x$level) / 2
                           )) *
                             se[1:nvar]),
                           exp(coef[indxAlpha] + mlevel * c(se_alpha))))
        }


      } else{
        index <- c(1:nvar)

        tmp_new <- cbind(exp(coef[1:nvar]),
                         exp(coef[1:nvar] - abs(qnorm((1 - x$level) / 2
                         )) * se[1:nvar]),
                         exp(coef[1:nvar] + abs(qnorm((1 - x$level) / 2
                         )) * se[1:nvar]))
      }
      dimnames(tmp_new) <- list(names(coef[index]),
                                c(
                                  "exp(coef)",
                                  paste("lower", x$level, sep = " "),
                                  paste("upper", x$level, sep = " ")
                                ))



      cat("\n")
      if (x$pophaz != "classic") {
        if (x$pophaz == "rescaled") {
          cat("\n")
          cat(
            "Excess hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
          )
          lines_al <-
            which(stringr::str_detect(rownames(tmp_new), pattern = "alpha"))
          print(tmp_new[-c(lines_al), ], digits = digits + 1)
          cat("\n")
          cat("and rescaled parameter on population hazard \n")

          tmp_new_alpha <- matrix(tmp_new[c(lines_al), ], nrow = 1)
          dimnames(tmp_new_alpha)[[1]] <- "alpha"
          dimnames(tmp_new_alpha)[[2]] <- dimnames(tmp_new)[[2]]
          print(tmp_new_alpha, digits = digits + 1)


        } else if (x$pophaz == "corrected" &
                   x$add.rmap.cut$breakpoint == FALSE) {
          cat("\n")
          cat(
            "Excess hazard hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
          )
          lines_al <-
            which(stringr::str_detect(rownames(tmp_new), pattern = "alpha"))
          print(tmp_new[-c(lines_al), ], digits = digits + 1)
          cat("\n")
          cat("and corrected scale parameters on population hazard \n")
          print(tmp_new[c(lines_al), ], digits = digits + 1)
        } else if (x$pophaz == "corrected" &
                   x$add.rmap.cut$breakpoint == TRUE) {
          cat("\n")
          cat(
            "Excess hazard hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
          )
          lines_al <-
            which(stringr::str_detect(rownames(tmp_new), pattern = "alpha"))
          print(tmp_new[-c(lines_al), ], digits = digits + 1)
          cat("\n")
          cat(
            "and corrected scale parameters on population hazard \n (non proportional correction using breakpoint approach)\n"
          )
          cat("\n")

          print(tmp_new[c(lines_al), ], digits = digits + 1)

          if (!is.na(x$add.rmap.cut$cut[1])) {
            n_break <- length(x$add.rmap.cut$cut)
            n_int <- 1 + n_break

          }

          cat("\n")

        }

      } else{
        cat("\n")
        cat(
          "Excess hazard hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
        )
        print(tmp_new, digits = digits + 1)
      }


      logtest <- (-2 * (x$loglik[1] - x$loglik[2]))
      df <- length(x$coef)
      cat("\n")
      cat("number of observations:",
          paste0(format(x$n), "; "),
          "number of events:",
          x$n.events)
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

    if (sum(x$cov.test) > 0) {
      cat("\n")
      cat("Results of tests for '",
          names(coef)[grep("T", as.character(x$cov.test))], "' equal to 0")
      cat("\n")
      cat(
        "  Likelihood ratio test=",
        format(round(x$loglik.test, 2)),
        "  on ",
        x$cov.df,
        " degree(s) of freedom,",
        " p=",
        format(1 - pchisq(x$loglik.test, x$cov.df)),
        sep = ""
      )
      cat("\n")
      cat(
        "  Wald test=",
        format(round(x$wald.test, 2)),
        "  on ",
        x$cov.df,
        " degree(s) of freedom,",
        " p=",
        format(1 - pchisq(x$wald.test, x$cov.df)),
        sep = ""
      )
      cat("\n")
      cat(
        "  Score test=",
        format(round(x$score.test, 2)),
        "  on ",
        x$cov.df,
        " degree(s) of freedom,",
        " p=",
        format(1 - pchisq(x$score.test, x$cov.df)),
        sep = ""
      )
      cat("\n")
      cat("number of observations:",
          format(x$n),
          "number of events:",
          x$n.events)
      cat("\n")



    }
    if (any(tmp_new[, "lower 0.95"] < 0)) {
      warning("\nlower 0.95 CI approximation may be incorrect")
    }
    invisible()
  }
