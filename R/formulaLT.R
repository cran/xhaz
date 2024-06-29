formulaLT <- function (formula, data, expected, expectedCum,
                       pophaz, base, degree, knots, bound, n.gleg, init,
                       random, n.aghq, recurrent, fnoptim, verbose, method,
                       iterlim, numHess, print.level, exactGradHess,
                       gradtol, testInit, keep.data, ...) {
  to_rescale <-  to_rescale2()
  if (pophaz == "classic") {

    res <- to_rescale(formula = formula, data = data,
                                   expected = expected,expectedCum =expectedCum,
                                   base = base, degree = degree, knots = knots,
                                   bound = bound, n.gleg = n.gleg, init = init,
                                   random = random, withAlpha = FALSE,
                                   n.aghq = n.aghq,
                                   verbose = verbose,iterlim = iterlim,
                                   print.level = print.level, gradtol = gradtol,
                                   keep.data = keep.data)
  }
  else if (pophaz == "rescaled") {
    res <- to_rescale(formula = formula, data = data,
                                   expected = expected,expectedCum =expectedCum,
                                   base = base,degree = degree, knots = knots,
                                   bound = bound, n.gleg = n.gleg, init = init,
                                   random = random, withAlpha = TRUE,
                                   n.aghq = n.aghq,
                                   verbose = verbose,iterlim = iterlim,
                                   print.level = print.level, gradtol = gradtol,
                                    keep.data = keep.data)
  }
return(res)

}
