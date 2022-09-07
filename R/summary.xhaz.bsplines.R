#' @title A summary.bsplines Function used to print a object of class \code{bsplines}
#'
#' @description This function present the estimated coefficients for the excess
#'  hazard baseline coefficient and for the covariate effects
#'
#' @param object an object of class \code{bsplines}
#'
#'
#' @param ... additionnal parameters which can be used in the \code{summary}
#' function
#'
#' @return Estimated parameters of the model in different scales for interpretation purposes.
#'
#' @keywords summary.bsplines
#'
#' @seealso \code{\link{xhaz}}, \code{\link{summary.bsplines}}, \code{\link{plot.bsplines}}
#'
#' @references Goungounga JA, Touraine C, Grafféo N, Giorgi R;
#' CENSUR working survival group. Correcting for misclassification
#' and selection effects in estimating net survival in clinical trials.
#'  BMC Med Res Methodol. 2019 May 16;19(1):104.
#'   doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#' Touraine C, Grafféo N, Giorgi R; CENSUR working survival group.
#' More accurate cancer-related excess mortality through correcting
#' background mortality for extra variables.
#'  Stat Methods Med Res. 2020 Jan;29(1):122-136.
#'  doi: 10.1177/0962280218823234. Epub 2019 Jan 23. PMID: 30674229.
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/30674229/}{PubMed})
#'
#'  Mba RD, Goungounga JA, Grafféo N, Giorgi R; CENSUR working survival group.
#'   Correcting inaccurate background mortality in excess hazard models
#'   through breakpoints. BMC Med Res Methodol. 2020 Oct 29;20(1):268.
#'   doi: 10.1186/s12874-020-01139-z. PMID: 33121436; PMCID: PMC7596976.
#'   (\href{https://pubmed.ncbi.nlm.nih.gov/33121436/}{PubMed})
#'
#'
#' @examples
#'
#' \donttest{
#' library("xhaz")
#' library("survival")
#' library("numDeriv")
#' library("survexp.fr")
#' library("splines")
#'
#' data("dataCancer", package = "xhaz")   # load the data set in the package
#'
#' fit.phBS <- xhaz(
#'               formula = Surv(obs_time_year, event) ~ ageCentre + immuno_trt,
#'               data = dataCancer, ratetable = survexp.fr,
#'               interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
#'               rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
#'               baseline = "bsplines", pophaz = "classic")
#'
#' summary(fit.phBS)
#' }
#'
#' @importFrom stats printCoefmat
#'
#' @export

summary.bsplines <- function(object, ...)
{
  digits <- max(options()$digits - 4, 3)
  cl <- try(object$call)
  if (!is.null(cl)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  if (!is.null(object$fail)) {
    cat(" xhaz.bsplines failed.", object$fail, "\n")
    return()
  }

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  if (!is.null(object$add.rmap)) {
    nalpha <- nlevels(object$add.rmap)
  } else{
    nalpha <- 0
  }
  if (nalpha) {
    if (nalpha == 1) {
      indxAlpha <- which(stringr::str_detect(names(object$coefficients),
                                             pattern = "alpha"))
      object$mycoef <- object$coefficients
      names(object$mycoef)[indxAlpha]  <- "alpha"
      names(object$coefficients)[indxAlpha] <- "log(alpha)"

    }
    else{
      object$mycoef <- object$coefficients

      indxAlpha <-
        which(stringr::str_detect(names(object$coefficients),
                                  pattern = "alpha"))
      names(object$coefficients)[c(indxAlpha)] <- paste("log(",
                                                        paste0('alpha.',
                                                               levels(object$add.rmap)),
                                                        ")",
                                                        sep = "")
      names(object$mycoef)[c(indxAlpha)] <-
        paste0('alpha.', levels(object$add.rmap))

    }
    coef <- c(object$coefficients)
  } else{
    coef <- c((object$coefficients))
  }
  object$var <- object$varcov

  # WARNING: works for one single alpha
  if (nalpha) {
    se <- sqrt((c(diag(object$var[1:(length(diag(object$var)) - 1),
                                  1:(length(diag(object$var)) - 1)]),
                  (object$var[length(diag(object$var)),
                              length(diag(object$var))]))))
  }
  else{
    se <- sqrt((diag(object$var)))
  }

  if (is.null(coef) | is.null(se))
    stop("Input is not valid")

  if (nalpha) {
    tmp <- cbind(coef,
                 se,
                 (c((object$coef[1:length(object$coef) - 1]),
                    (object$coef[length(object$coef)])) -
                    abs(qnorm((1 - object$level) / 2
                    )) * sqrt((
                      diag(object$var)
                    ))),
                 (c((object$coef[1:length(object$coef) - 1]),
                    (object$coef[length(object$coef)])) +
                    abs(qnorm((1 - object$level) / 2
                    )) * sqrt((
                      diag(object$var)
                    ))),
                 coef / se,
                 signif(1 - pchisq((coef / se) ^ 2, 1), digits - 1))

    dimnames(tmp) <- list(names(coef),
                          c(
                            "coef",
                            "se(coef)",
                            paste("lower", object$level, sep = " "),
                            paste("upper", object$level, sep = " "),
                            "z",
                            "Pr(>|z|)"
                          ))
  } else{
    tmp <-
      cbind(coef,
            se,
            (coef - abs(qnorm((1 - object$level) / 2
            )) * se),
            (coef + abs(qnorm((1 - object$level) / 2
            )) * se),
            coef / se,
            signif(1 - pchisq((coef / se) ^ 2, 1), digits - 1))
    dimnames(tmp) <- list(names(coef),
                          c(
                            "coef",
                            "se(coef)",
                            paste("lower", object$level, sep = " "),
                            paste("upper", object$level, sep = " "),
                            "z",
                            "Pr(>|z|)"
                          ))
  }

  cat("\n")

  printCoefmat(
    tmp,
    P.values = TRUE,
    digits = digits,
    signif.stars = TRUE,
    na.print = "NA",
    ...
  )


  nPH <- object$nPH
  nTD <- object$nTD
  if (nPH > 0) {
    tmp <-
      exp(tmp)[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH + nalpha), c(1, 3, 4)]


    if (is.null(dim(tmp))) {
      tmp <- t(as.matrix(tmp))
    }
    dimnames(tmp) <-
      list(c(names(coef[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]),
             names(object$mycoef[(5 + 5 * nTD + nPH + 1):(5 + 5 * nTD + nPH + nalpha)])),
           c(
             "exp(coef)",
             paste("lower", object$level, sep = " "),
             paste("upper", object$level, sep = " ")
           ))
    cat("\n")

    if (object$pophaz != "classic") {
      if (object$pophaz == "rescaled") {
        cat(
          "Excess hazard and expected hazard ratio(s) \n(proportional effect variable(s) for exess hazard ratio(s) )\n"
        )
        print(tmp, digits = digits + 1)


        cat("\n")
        cat(
          "Excess hazard ratio(s)\n(proportional effect variable(s) for exess hazard ratio(s))\n"
        )
        lines_al <-
          which(stringr::str_detect(rownames(tmp), pattern = "alpha"))
        print(tmp[-c(lines_al), ], digits = digits + 1)
        cat("\n")
        cat("and rescaled parameter on population hazard \n")

        tmp_new_alpha <- matrix(tmp[c(lines_al), ], nrow = 1)
        dimnames(tmp_new_alpha)[[1]] <- "alpha"
        dimnames(tmp_new_alpha)[[2]] <- dimnames(tmp)[[2]]
        print(tmp_new_alpha, digits = digits + 1)

      }

    } else{
      cat(
        "Excess hazard ratio(s) \n(proportional effect variable(s) for exess hazard ratio(s) )\n"
      )
      print(tmp, digits = digits + 1)
    }

  }

  cat("\n")
  cat(
    "number of observations:",
    paste0(format(object$n), "; "),
    "number of events:",
    object$n.events
  )
  cat("\n")
  cat(
    "log-likelihood: ",
    format(object$loglik),
    " (for ",
    length(object$coef),
    " degree(s) of freedom)"
  )
  cat("\n")



  if (sum(object$cov.test) > 0) {
    cat("\n")
    cat(
      "  Likelihood ratio test of PH effect for '",
      gsub("\\(|\\)", "",
           as.character(stringr::str_remove((
             attr(object$terms, "term.labels")[object$cov.test]
           ),
           "qbs"))),
      "'=",
      format(round(object$loglik.test, 2)),
      ",\n on ",
      object$cov.df,
      " degree(s) of freedom,"
    )
    cat(" p=",
        format(1 - pchisq(object$loglik.test, object$cov.df)),  sep =
          "")
  }
  invisible()
}
