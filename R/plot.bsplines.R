#' @title plot.bsplines
#'
#' @description to plot the log hazard ratio functions for non-proportional
#'  hazards model
#'
#'
#' @param x An object of class xhaz
#'
#' @param cov specify covariates for which a plot is required.
#'
#' @param conf.int a vector of logical values indicating whether (if TRUE)
#' confidence intervals will be plotted. The default is to do so if the plot
#' concerns only one curve.
#'
#' @param baseline a vector of logical values indicating whether (if \code{baseline = TRUE})
#' to plot the curve for the baseline group. Default is FALSE, except if cov
#' is unspecified.
#'
#' @param xrange vector indicating the minimum and the maximum values of the
#' x axis. By default, these values are automatically calculated for the first
#'    plot (i.e before the use of add argument).
#'
#' @param  yrange vector indicating the minimum and the maximum values of the y
#'  axis. By default, these values are automatically calculated for the
#'  first plot (i.e before the use of add argument).
#'
#' @param  xlegend value indicating the location of the legend over x axis.
#'  By default, location at the left of the plot.
#'
#' @param ylegend value indicating the location of the legend over y axis.
#'  By default, location at the top of the plot
#'
#' @param glegend vectors of names attributed to each lines of the excess hazard
#' to be displayed in the plot. If (\code{baseline = TRUE}), glegend is \code{"baseline"}.
#'
#' @param xaxs the x axis style, as listed in 'par'. Survival curves are
#' traditionally drawn with the curve touching the bounding box on the left
#' edge, but not touching it on the right edge. This corresponds to neither
#' of the two standard S axis styles of "e" (neither touches) or "i" (both touch).
#' If xaxis is missing or NULL the internal axis style is used (xaxs= i) but
#' only after the right endpoint has been extended.
#'
#' @param add a logical value indicating whether to add the survival curves to the
#' current plot (if \code{add = TRUE}). Default is FALSE.
#'
#' @param col a vector of integers specifying colors for each curve. The default
#' value is 1.
#'
#' @param lty a vector of integers specifying line types for each curve. The
#' default value is fixed by the number of covariates (plus 1 if \code{baseline = TRUE}).
#'
#' @param lwd a vector of numeric values for line widths. The default value is 1.
#'
#' @param ... additional arguments affecting the plot function
#'
#' @keywords plot.bsplines
#'
#' @return The return of this function produce graphics of log hazard ratio
#' functions for non-proportional hazards model
#'
#' @author Juste Goungounga, Robert Darlin Mba, Nathalie Grafféo and Roch Giorgi
#' @export
#'
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
#' Mba RD, Goungounga JA, Grafféo N, Giorgi R; CENSUR working survival group.
#' Correcting inaccurate background mortality in excess hazard models
#' through breakpoints. BMC Med Res Methodol. 2020 Oct 29;20(1):268.
#' doi: 10.1186/s12874-020-01139-z. PMID: 33121436; PMCID: PMC7596976.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/33121436/}{PubMed})
#'
#'
#' Giorgi R, Abrahamowicz M, Quantin C, Bolard P, Esteve J, Gouvernet J,
#' Faivre J. A relative survival regression model using B-spline functions
#' to model non-proportional hazards.
#' Statistics in Medicine 2003; 22: 2767-84.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/12939785/}{PubMed})
#'
#' @examples
#' \donttest{
#' # load the data set in the package
#' library("xhaz")
#' library("survexp.fr")
#'
#'data("dataCancer", package = "xhaz")   # load the data set in the package
#'
#' fit.nphBS <- xhaz(
#'       formula = Surv(obs_time_year, event) ~ ageCentre + qbs(immuno_trt),
#'       data = dataCancer,
#'       ratetable = survexp.fr,
#'       interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
#'       rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
#'       baseline = "bsplines", pophaz  = "classic")
#'
#'  plot(fit.nphBS, cov = "immuno_trt", col = "blue", baseline = FALSE)
#' }
#' @importFrom graphics plot lines grid legend

plot.bsplines <- function(x,
                          cov,
                          conf.int = TRUE,
                          baseline = FALSE,
                          xrange,
                          yrange,
                          xlegend,
                          ylegend,
                          glegend,
                          xaxs = NULL,
                          add = FALSE,
                          col = 1,
                          lty = 1,
                          lwd = 1,
                          ...) {
  rsfit <- x
  covlev <- levels(rsfit$data[, cov])
  if (inherits(rsfit, "bsplines")) {
    int <- rsfit$int
    n.tau <- 5
    id_coefvov <- (which(stringr::str_detect(names(
      rsfit$coefficients
    ),
    pattern = cov)))
    nTD <- length(id_coefvov) / 5
    if (missing(cov)) {
      if (ncol(attr(rsfit$terms, "factors")) != sum(rsfit$bsplines))
        stop(
          "At least one covariate is not TD. You must have to specify the TD covariates for which plot is required."
        )
      coef <- rsfit$coef
      ncov <-
        (length(coef) - n.tau) / 5	#/5 because there are five parameters for one covariate TD
      baseline <- TRUE
      npar <- ncov + 1		#All the tested variables plus the baseline
      cov <- (n.tau + 1):(length(coef))		#Location of the variables
      names.cov <- names(coef)[cov]
    }
    else {
      coef <- rsfit$coef
      if (is.character(cov)) {
        dum <- 0
        for (i in 1:length(cov))
          dum <- sum(length(grep(cov[i],
                                 as.character(
                                   names(rsfit$coef)
                                 ))) / n.tau + dum)
        #if (dum != length(cov))
        if (rsfit$nTD == 0)
          stop(
            "At least one covariate is not TD. You must have to specify the TD covariates for which plot is required."
          )

        npar <-
          length(cov)		#Used to take account when the baseline is needed.
        covv <- c(1:length(cov))

        for (i in 1:length(cov)) {
          covv[i] <-
            grep(stringr::str_sub(cov[i], start = 1, end = nchar(cov[i]) - 1),
                 as.character(attr(rsfit$terms, "term.labels")))
        }
        names.cov <- attr(rsfit$terms, "term.labels")[covv]
        id_cov <- which(stringr::str_detect(names(rsfit$coefficients), pattern = cov))
        cov <- names(rsfit$coefficients)[id_cov]
      }
      else{
        if (length(grep("T", as.character(rsfit$bsplines[cov]))) != length(rsfit$bsplines[cov]))
          stop(
            "At least one covariate is not TD. You must have to specify the TD covariates for which plot is required."
          )

        npar <-
          length(cov)		#Used to take account when the baseline is needed.
        names.cov <- attr(rsfit$terms, "term.labels")[cov]
        cov <-
          c(sapply(1:npar, function(i, names.cov, coef)
            (grep(names.cov[i], as.character(names(
              coef
            )))), names.cov = names.cov, coef = rsfit$coef))
      }
      ncov <-
        length(cov) / 5	#There are 5 coefficients for one covariable
      coef <- c(coef[1:n.tau], coef[cov])
      var <-
        rsfit$var[c(1:(n.tau + sum(length(id_coefvov) / 5) * 5)), c(1:(n.tau +
                                                                         sum(length(id_coefvov) / 5) * 5))]

      if (baseline == TRUE)
        npar <-
        length(cov) + 1		#used to take account when the baseline is needed.
      if (any(is.na(cov)) ||
          length(cov) > ncov * 5 || length(cov) < 1)
        stop("Invalid variable requested")
    }



    if (missing(conf.int)) {
      if (npar > 1) {
        conf.int <- FALSE
      } else {
        conf.int <- TRUE
      }
    } else
      conf.int <- TRUE

    if (length(conf.int) < npar)
      conf.int <- c(conf.int, rep(FALSE, npar - length(conf.int)))
    else if (length(conf.int) > npar)
      conf.int <- conf.int[seq(npar)]

    rsfit.coef <-
      c(rsfit$coef[cov], rsfit$coef[(ncol(attr(rsfit$terms, "factors")) + 1):(length(rsfit$coef))])

    x.time <- seq(0, int[4], 0.01) #c(0:int[4])
    knots <- c(0, 0, int[1], int[2], int[3], int[4], int[4], int[4])
    splinesbase <- splines::spline.des(knots, x.time, 3)$design

    logHRbeta <- matrix(rep(0, length(x.time) * ncov), ncol = ncov)
    for (n in 1:ncov) {
      dum <- t(coef[((n - 1) * 5 + 6):((n - 1) * 5 + 10)] * t(splinesbase))
      logHRbeta[, n] <- rowSums(dum)
    }

    if (baseline == FALSE) {
      logHRi <- logHRbeta
    }
    else{
      dum <- t(coef[1:n.tau] * t(splinesbase))
      logHRtau <- rowSums(dum)
      logHRi <- cbind(logHRtau, logHRbeta)
    }

    if (missing(xrange))
      xrange <- c(0, max(x.time))

    if (is.null(xaxs)) {
      xrange <- 1.04 * xrange
      xaxs <- "i"
    }
    ncov2 <- ncol(logHRi)
    if (length(col) != npar) {
      col <- rep(col, length = ncov2)
    }

    if (missing(lty))
      lty <- seq(ncov2)
    else if (length(lty) != ncov2)
      lty <- rep(lty, length = ncov2)
    if (length(lwd) != npar)
      lwd <- rep(lwd, length = ncov2)



    if (any(conf.int)) {
      lcibeta <- matrix(rep(0, length(x.time) * ncol(logHRi)), ncol = ncov)
      lcubeta <-
        matrix(rep(0, length(x.time) * ncol(logHRi)), ncol = ncov)

      varB <- diag(var)[(n.tau + 1):(n.tau + n.tau * nTD)]
      VcovB <-
        var[(n.tau + 1):(n.tau + n.tau * nTD), (n.tau + 1):(n.tau + n.tau * nTD)]

      for (n in 1:ncov) {
        varBB <-
          c(varB[n], varB[n + nTD], varB[n + (2 * nTD)], varB[n + (3 * nTD)], varB[n +
                                                                                     (4 * nTD)])
        w1 <- t(varBB * t((splinesbase ^ 2)))
        w1 <- rowSums(w1)
        w2 <- {
          0
        }
        for (i in 1:(n.tau - 1)) {
          for (j in (i + 1):n.tau) {
            w2 <-
              w2 + 2 * ((splinesbase[, i] * splinesbase[, j]) * VcovB[((i - 1) * ncov +
                                                                         n), ((j - 1) * ncov + n)])
          }
          lcibeta[, n] <-
            (logHRbeta[, n] - abs(qnorm((
              1 - rsfit$level
            ) / 2)) * sqrt(w1 + w2))
          lcubeta[, n] <-
            (logHRbeta[, n] + abs(qnorm((
              1 - rsfit$level
            ) / 2)) * sqrt(w1 + w2))
        }
      }
      if (baseline == FALSE) {
        lcii <- lcibeta
        lcui <- lcubeta
      }
      else{
        vartau <- diag(var)[1:n.tau]
        Vcovtau <- var[1:n.tau, 1:n.tau]
        w1t <- t(vartau * t(splinesbase ^ 2))
        w1t <- rowSums(w1t)
        w2t	<-	{
          0
        }
        for (i in 1:(n.tau - 1)) {
          for (j in (i + 1):n.tau) {
            w2t <- w2t + 2 * ((splinesbase[, i] * splinesbase[, j]) * Vcovtau[j, i])
          }
        }
        lcitau <-
          (logHRtau - abs(qnorm((
            1 - rsfit$level
          ) / 2)) * sqrt(w1t + w2t))
        lcutau <-
          (logHRtau + abs(qnorm((
            1 - rsfit$level
          ) / 2)) * sqrt(w1t + w2t))
        lcii <- cbind(lcitau, lcibeta)
        lcui <- cbind(lcutau, lcubeta)
      }
      if (!add) {
        if (missing(yrange))
          yrange <-
            c(min(logHRi, lcii, lcui) * 0.96, c(max(logHRi, lcii, lcui)) * 0.96)
        plot(xrange, yrange, type = "n", xaxs = xaxs, ...)
      }
    }
    else if (!add) {
      if (missing(yrange))
        yrange <- c(min(logHRi) * 0.96, c(max(logHRi)) * 0.96)
      plot(xrange,  yrange, type = "n", xaxs = xaxs,  ...)
    }

    if (any(conf.int)) {
      sapply(1:ncol(logHRi), function(i) {
        lines(x.time,
              logHRi[, i],
              lty = lty[i],
              col = col[i],
              lwd = lwd[i])
        lines(
          x.time,
          lcii[, i],
          lty = ifelse(ncov == 1,
                       lty[i] +
                         ifelse(add, 0, 1),
                       lty[i]),
          col = col[i],
          lwd = lwd[i]
        )

        lines(
          x.time,
          lcui[, i],
          lty = ifelse(ncov == 1,
                       lty[i] +
                         ifelse(add, 0, 1), lty[i]),
          col = col[i],
          lwd = lwd[i]
        )
      })
    } else {
      sapply(1:ncol(logHRi), function(i) {
        lines(x.time,
              logHRi[, i],
              lty = lty[i],
              col = col[i],
              lwd = lwd[i])
      })
    }


    if (missing(glegend)) {
      glegend <- paste(names.cov,
                       c(covlev)[2:(length(covlev))])
    }

    if (baseline == TRUE) {
      names.cov <- c(c("baseline"), glegend)
    } else {
      names.cov <- glegend
    }

    if (missing(ylegend) & missing(xlegend)) {
      legend(
        "bottomleft",
        legend = names.cov[1:ncol(logHRi)],
        lty = lty,
        lwd = lwd,
        col = col,
        bty = "n"
      )
    } else{
      legend(
        x = xlegend,
        y = ylegend,
        legend = names.cov[1:ncol(logHRi)],
        lty = lty,
        lwd = lwd,
        col = col,
        bty = "n"
      )
    }

  } else{
    stop("only implemented for time-dependent covariate effect")
  }
  invisible()

}
