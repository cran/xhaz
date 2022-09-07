
testM <- function(X,
                  Y,
                  ehazard,
                  ehazardInt,
                  int,
                  covtest,
                  bsplines = bsplines,
                  init,
                  control, event, Terms, strats,
                  add.rmap, add.rmap.cut, ageDiag, ageDC,
                  optim, trace, speedy, data) {

  if (min(add.rmap.cut$cut) < min(c(data$age, data$age + data$time))) {
    if (max(add.rmap.cut$cut) <= max(c(data$age, data$age + data$time))) {
      stop("Breakpoint(s) is (are) smaller than the minimum age")
    } else
      stop(
        "Breakpoint(s) is (are) smaller than the minimum age and breakpoint(s) greater than the maximum age"
      )
  } else{
    if (max(add.rmap.cut$cut) > max(c(data$age, data$age + data$time)))
      stop("Breakpoint(s) is (are) greater than the maximum age")
  }
  fitter <- get("esteve.ph.fit")


  res <- try(fitter(X,
                    Y,
                    ehazard,
                    ehazardInt,
                    int = int,
                    covtest,
                    bsplines =  bsplines,
                    init,
                    control, event, Terms, strats,
                    add.rmap, add.rmap.cut, ageDiag, ageDC,
                    optim, trace, speedy), TRUE)

  if (class(res)[[1]] != "try-error") {
    res$AIC <- (2 * length(res$coefficients) - 2 * (res$loglik))[2]
    res$BIC <- (log(length(Y)) * length(res$coefficients) - 2 * (res$loglik))[2]
  }

  res
}
