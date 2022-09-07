#' @import statmod

giorgi.tdph.fit <- function(x, y, ehazard, ehazardInt, int, covtest, bsplines,
                            init, control, event, Terms, strats, add.rmap,
                            add.rmap.cut, ageDiag, ageDC, optim, trace, speedy,
                            nghq = nghq) {
  k <- 3
  nrowx <- nrow(x)
  ehazard <- c(ehazard)
  event <- c(event)
  y[, 1] <- c(y[, 1])# / 12)
  cst <- 1
  cpti <- 1
  cptj <- 1
  int <- int# / 12
  alpha <- NULL
  #Coefficients for the basis functions
  int23 <- int[2] * int[3]
  int4232 <- (int[4] - int[2]) * (int[3] - int[2])
  int4243 <- (int[4] - int[2]) * (int[4] - int[3])
  int424353 <- int4243 * (int[4] - int[3])
  int5343 <- (int[4] - int[3]) ^ 2
  #Basis functions
  spline1 <- c(1, -2 / int[2], 1 / (int[2] ^ 2),
               0, 0, 0, 0, 0, 0)

  spline2 <- c(0,
               (1 / int[2]) + (int[3] / int23),
               (-1 / int[2] ^ 2) - (1 / int23),
               (int[3] / (int[3] - int[2])),
               -2 / (int[3] - int[2]),
               1 / (int[3] * (int[3] - int[2])),
               0, 0, 0)

  spline3 <- c(0, 0,
               1 / int23,
               -((int[4] * int[2]) / int4232),
               (1 / (int[3] - int[2])) +
                 (int[4] / int4232) +
                 (int[2] / int4232),
               -1 / (int[3] * (int[3] - int[2])) - 1 / int4232,
               int[4] ^ 2 / int4243,
               -2 * int[4] / int4243,
               1 / int4243)

  spline4 <- c(0, 0, 0,
               int[2] ^ 2 / int4232,
               -2 * int[2] / int4232,
               1 / int4232,
               -(int[2] * int[4] * (int[4] - int[3]) +
                   int[3] * int[4] * (int[4] - int[2])) / int424353,
               (2 * int[4] * int[4] - 2 * int[2] * int[3]) / int424353,
               (int[2] + int[3] - 2 * int[4]) / int424353)

  spline5 <- c(0, 0, 0, 0, 0, 0,
               int[3] ^ 2 / int5343,
               -2 * int[3] / int5343,
               1 / int5343)

  #The resulting B-spline function
  p <- matrix(c(spline1, spline2, spline3, spline4, spline5),
              nrow = 5,
              byrow = T)
  knot <- c(int[c(-1,-length(int))])
  delta <- sort(c(rep(c(int[1], int[length(int)]), k), knot))

  P <- splines::splineDesign(knots = delta, x = y[, 1], ord = k)
  #Boundaries of the integral
  boundmax <- matrix(0, nrow(y), k)
  boundmin <- matrix(0, nrow(y), k)

  for (i in 1:k) {
    boundmax[, i] <- (y[, 1] > int[i + 1]) * int[i + 1] +
      (y[, 1] <= int[i + 1] & y[, 1] >= int[i]) * y[, 1]
    boundmin[, i] <- (y[, 1] >= int[i]) * int[i]
  }
  if (!is.null(add.rmap)) {
    if (!is.factor(add.rmap)) {
      stop("The alpha argument must be a factor.")
    } else{
      nalpha <- nlevels(add.rmap)
    }
  } else{
    nalpha <- 0
  }
  #Reorganize 'x' with TD covariates following by PH covariates.
  theta0 <- c(rep(0, nalpha + 5 + 5 * (length(bsplines))))
  thetaPH <- NULL
  nTD <- ncol(as.matrix(x[, grep("TRUE", as.character(bsplines))]))
  nPH <- (length(bsplines) - nTD)

  if (nalpha) {
    #Initialization of theta0
    if (!is.null(init)) {
      if (length(unlist(init)) != (5 + nTD * 5 + nPH + nalpha))
        stop("The number of initials values must the same as
               \nthe number of parameters to estimate.")
    }
  } else{
    #Initialization of theta0
    if (!is.null(init)) {
      if (length(unlist(init)) != (5 + nTD * 5 + nPH))
        stop("The number of initials values must the same as
               \nthe number of parameters to estimate.")
    }
  }
x_new <- x
  if (nTD != length(bsplines)) {
    if (!is.null(init)) {
      if (nPH != 0) {
        dummyF <- dimnames(x)[[2]][grep("FALSE", as.character(bsplines))]
        dummyF <- sapply(1:nPH,
                         function(i, init, dummyF)
                           unlist(init[grep(dummyF[i], names(init))]),
                         init = init,
                         dummyF = dummyF)
      }
      else {
        dummyF <- NULL
      }

      if (nTD != 0) {
        dummyT <- dimnames(x)[[2]][grep("TRUE", as.character(bsplines))]
        dummyT <- sapply(1:nTD,
                         function(i, init, dummyT)
                           unlist(init[grep(dummyT[i], names(init))]),
                         init = init,
                         dummyT = dummyT)
      }
      else{
        dummyT <- NULL
      }
indxAlpha <- (5 + nTD * 5 + nPH + 1) : (5 + nTD * 5 + nPH + nalpha)
      if (nalpha) {
        if (nTD != 0 & nPH != 0) {
          theta0 <- c(unlist(init[1][1:5]),
                      dummyT,
                      dummyF,
                      unlist(init[length(init)]))
        }
        else{
          if (nTD == 0 &
              nPH != 0) {
            theta0 <-
              c(unlist(init[1][1:5]), dummyF, unlist(init[indxAlpha#length(init)
                                                          ]))
          }
          if (nTD != 0 &
              nPH == 0) {
            theta0 <-
              c(unlist(init[1][1:5]), dummyT, unlist(init[indxAlpha#length(init)
                                                          ]))
          }
        }
      } else{
        if (nTD != 0 & nPH != 0) {
          theta0 <- c(unlist(init[1][1:5]), dummyT, dummyF)
        } else{
          if (nTD == 0 & nPH != 0) {
            theta0 <- c(unlist(init[1][1:5]), dummyF)
          }
          if (nTD != 0 &
              nPH == 0) {
            theta0 <- c(unlist(init[1][1:5]), dummyT)
          }
        }
      }

      names(theta0) <- NULL
      thetaPH <- theta0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
    }
    else{
      theta0 <- c(rep(0, 5 + 5 * nTD + nPH + nalpha))
      thetaPH <- theta0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
    }

    if (nTD != 0 & nPH != 0) {
      namesx <- c(dimnames(x)[[2]][grep("TRUE", as.character(bsplines))],
                  dimnames(x)[[2]][grep("FALSE", as.character(bsplines))])
      x <- cbind(as.matrix(x[, grep("TRUE", as.character(bsplines))]),
                 as.matrix(x[, grep("FALSE", as.character(bsplines))]))
    }
    else{
      if (nTD == 0 & nPH != 0) {
        namesx <- c(dimnames(x)[[2]][grep("FALSE", as.character(bsplines))])
        x <-  as.matrix(x[, grep("FALSE", as.character(bsplines))])
      }
      if (nTD != 0 & nPH == 0) {
        namesx <- c(dimnames(x)[[2]][grep("TRUE", as.character(bsplines))])
        x <-  as.matrix(x[, grep("TRUE", as.character(bsplines))])
      }
    }
    dimnames(x)[[2]] <- namesx
  }
  else{
    if (!is.null(init)) {
      theta0 <- c(unlist(init))
      names(theta0) <- NULL
      if (nalpha) {
        thet0 <- theta0[5 + 5 * nTD + nPH]
        if (nPH != 0) {
          thetaPH <- thet0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
        }
        else{
          thetaPH <- 0
        }
      }
      else{
        if (nPH != 0) {
          thetaPH <- theta0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
        } else{
          thetaPH <- 0
        }
      }
    } else{
      if (nalpha) {
        theta0 <- c(rep(0, 5 + 5 * nTD + nPH + nalpha))
        thet0 <- theta0[5 + 5 * nTD + nPH]
        thetaPH <- thet0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
      }
      else{
        theta0 <- c(rep(0, 5 + 5 * nTD + nPH))
        thetaPH <- theta0[(5 + 5 * nTD + 1):(5 + 5 * nTD + nPH)]
      }
    }

    namesx <- c(dimnames(x)[[2]][grep("TRUE", as.character(bsplines))],
                dimnames(x)[[2]][grep("FALSE", as.character(bsplines))])
    x <- cbind(as.matrix(x[, grep("TRUE", as.character(bsplines))]),
               as.matrix(x[, grep("FALSE", as.character(bsplines))]))

    dimnames(x)[[2]] <- namesx
  }





#Integral function using Gauss-Legendre quadrature
IntGL <- function(f, bound, cst, cpti, cptj, theta, x, nTD, p, nghq = nghq) {
  # default nghq=12
  GL <- gauss.quad(n = nghq, kind = "legendre")
  e <- matrix(GL$nodes,ncol = 1)
  w <- matrix(GL$weights,ncol = 1)
  dif <- t(bound[1, ]) - t(bound[2, ])
  addbound <- t(bound[2, ]) + t(bound[1, ])

  nodes2 <- 0.5 * (matrix(rep(addbound, nghq), ncol = nghq) + t(e %*% dif))
  weights2 <- t(as.vector(w) %*% (dif))
  GL$nodes2 <- c(nodes2)
  GL$weights2 <- weights2
  fx <- matrix(c(rep(0, nrow(nodes2) * ncol(nodes2))), ncol = ncol(nodes2)) +
    f(nodes2, cst, cpti, cptj, theta, x, nTD, p)

  fx <- colSums(t(as.matrix((fx))*weights2))
  fx
}

  #Used for the integral calculus of the disease-related mortality
  #hazard function ("corrected hazard")
  Int.chazard <- function(u, cst, cpti, cptj, theta, x, nTD, p) {
    tau.bs <- 0
    beta.bs <- 0
    nghq <- dim(u)[2]

    beta.bs <- matrix(0, nrow(x), nghq)

    for (i in 1:5) {
      #
      tau.bs <- tau.bs +
        theta[i] * (p[i, ((cst - 1) * 3 + 1)] +
                      p[i, ((cst - 1) * 3 + 2)] * u +
                      p[i, ((cst - 1) * 3 + 3)] * u ^ 2)
    }
    #The number of col is 12 or nghq because the shape of IntGL function
    if (nTD != 0) {
      for (niv in 1:nTD) {
        for (h in 1:5) {
          beta.bs <- beta.bs +
            theta[6 + ((h - 1) * nTD) + niv - 1] *
            x[, niv] * (p[h, ((cst - 1) * 3 + 1)] +
                          p[h, ((cst - 1) * 3 + 2)] * u +
                          p[h, ((cst - 1) * 3 + 3)] * u ^ 2)
        }
      }
    }
    exp(tau.bs + beta.bs)
  }
  #Used for the integral calculus of the first derivative
  #of the baseline function
  Int.FDbase <- function(u, cst, cpti, cptj, theta, x, nTD, p) {
    tau.bs <- 0
    beta.bs <- 0
    nghq <- dim(u)[2]
    beta.bs <- matrix(0, nrow(x), nghq)

    for (i in 1:5) {
      tau.bs <- tau.bs +
        theta[i] *
        (p[i, ((cst - 1) * 3 + 1)] +
           p[i, ((cst - 1) * 3 + 2)] * u +
           p[i, ((cst - 1) * 3 + 3)] * u ^ 2)
    }
    if (nTD != 0) {
      for (niv in 1:nTD) {
        for (h in 1:5) {
          beta.bs <- beta.bs +
            theta[6 + ((h - 1) * nTD) + niv - 1] *
            x[, niv] * (p[h, ((cst - 1) * 3 + 1)] +
                          p[h, ((cst - 1) * 3 + 2)] * u +
                          p[h, ((cst - 1) * 3 + 3)] * u ^ 2)
        }
      }
    }
    ((p[cpti, ((cst - 1) * 3 + 1)] +
        p[cpti, ((cst - 1) * 3 + 2)] * u +
        p[cpti, ((cst - 1) * 3 + 3)] * u ^ 2) *
        exp(tau.bs + beta.bs))
  }
  #Used for the integral calculus of the second derivative of the baseline function
  Int.SDbase <- function(u, cst, cpti, cptj, theta, x, nTD, p) {
    tau.bs <- 0
    beta.bs <- 0
    nghq <- dim(u)[2]
    beta.bs <- matrix(0, nrow(x), nghq)

    for (i in 1:5) {
      tau.bs <- tau.bs +
        theta[i] *
        (p[i, ((cst - 1) * 3 + 1)] +
           p[i, ((cst - 1) * 3 + 2)] * u +
           p[i, ((cst - 1) * 3 + 3)] * u ^ 2)
    }
    if (nTD != 0) {
      for (niv in 1:nTD) {
        for (h in 1:5) {
          beta.bs <- beta.bs +
            theta[6 + ((h - 1) * nTD) + niv - 1] *
            x[, niv] *
            (p[h, ((cst - 1) * 3 + 1)] +
               p[h, ((cst - 1) * 3 + 2)] * u +
               p[h, ((cst - 1) * 3 + 3)] * u ^ 2)
        }
      }
    }
    (p[cpti, ((cst - 1) * 3 + 1)] +
        p[cpti, ((cst - 1) * 3 + 2)] * u +
        p[cpti, ((cst - 1) * 3 + 3)] * u ^ 2) *
      (p[cptj, ((cst - 1) * 3 + 1)] +
         p[cptj, ((cst - 1) * 3 + 2)] * u +
         p[cptj, ((cst - 1) * 3 + 3)] * u ^ 2) *
      exp(tau.bs + beta.bs)
  }

  #Full model
  if (optim) {
    Fmodel <- giorgi.tdph.optim.maxim(x,
                                      theta0,
                                      nTD,
                                      nPH,
                                      event,
                                      ehazard,
                                      ehazardInt,
                                      P,
                                      p,
                                      k,
                                      nrowx,
                                      IntGL,
                                      cpti,
                                      cptj,
                                      cst,
                                      boundmin,
                                      boundmax,
                                      Int.chazard,
                                      Int.FDbase,
                                      Int.SDbase,
                                      int,
                                      control$iter.max,
                                      control$eps,
                                      Terms,
                                      add.rmap,
                                      trace,
                                      speedy,
                                      nghq = nghq)
    }
  else{
    Fmodel <- giorgi.tdph.maxim(x,
                                theta0,
                                nTD,
                                nPH,
                                event,
                                ehazard,
                                ehazardInt,
                                P,
                                p,
                                k,
                                nrowx,
                                IntGL,
                                cpti,
                                cptj,
                                cst,
                                boundmin,
                                boundmax,
                                Int.chazard,
                                Int.FDbase,
                                Int.SDbase,
                                int,
                                control, nghq = nghq)
  }
  #Tested model: if the likelihood ratio test of PH is required
  if (sum(covtest) > 0) {
    nPH.Full <- nPH
    nTD.Full <- nTD

    cov.test <- covtest
    theta0 <- Fmodel$coefficients
    nPH <- sum(covtest) + nPH
    nTD <- length(covtest) - nPH
    #Reorganize 'x' with TD covariates following by the new PH (tested)
    #and old PH covariates.
    dummy <- rep(1:length(bsplines), 1)
    attr(Terms, "term.labels") <- colnames(x)
    xTDarg <- ((bsplines == TRUE)*(covtest == FALSE))*dummy
    xPH1arg <- ((bsplines == TRUE)*(covtest == TRUE))*dummy
    xPH2arg <- ((bsplines == FALSE)*(covtest == FALSE))*dummy


    namesx <- c(attr(Terms, "term.labels")[xTDarg],
                attr(Terms, "term.labels")[xPH1arg],
                attr(Terms, "term.labels")[xPH2arg])

    xTD <- matrix(x[, attr(Terms, "term.labels")[xTDarg]])
    xPH1 <- matrix(x[, attr(Terms, "term.labels")[xPH1arg]])
    xPH2 <- matrix(x[, attr(Terms, "term.labels")[xPH2arg]])

    if (nrow(xTD) != 0) {
      x <- cbind(xTD, xPH1)
      if (nrow(xPH2) != 0) {
        x <- cbind(x, xPH2)
      }
    }
    else {
      x <- cbind(xPH1)
      if (nrow(xPH2) != 0) {
        x <- cbind(x, xPH2)
      }
    }
    dimnames(x)[[2]] <- namesx

    if (nPH.Full != 0)
      thetaPH <- theta0[(5 + 5 * nTD.Full + 1):length(theta0)]
    thetaPH <- c(rep(0, ncol(xPH1)), thetaPH)

    if (nTD >= 1) {
      vec <- c(rep(0, 5 * nTD))
      for (i in 1:nTD) {
        vectnTD <- ((i - 1) * 5 + 1):((i - 1) * 5 + 5)
        vec[vectnTD] <- grep(dimnames(x)[[2]][i], names(theta0))
      }
      thetaTD <- theta0[vec]
      if (nalpha) {
        theta0 <- c(theta0[1:5], thetaTD, thetaPH, alpha)
      } else{
        theta0 <- c(theta0[1:5], thetaTD, thetaPH)

      }
      names(theta0) <- NULL
    }
    else{
      theta0 <- c(theta0[1:5], thetaPH)
    }


    if (optim) {
      Tmodel <- giorgi.tdph.optim.maxim(x,
                                        theta0,
                                        nTD,
                                        nPH,
                                        event,
                                        ehazard,
                                        ehazardInt,
                                        P,
                                        p,
                                        k,
                                        nrowx,
                                        IntGL,
                                        cpti,
                                        cptj,
                                        cst,
                                        boundmin,
                                        boundmax,
                                        Int.chazard,
                                        Int.FDbase,
                                        Int.SDbase,
                                        int,
                                        control$iter.max,
                                        control$eps,
                                        Terms,
                                        add.rmap,
                                        trace,
                                        speedy,
                                        nghq = nghq)
    }
    else{
      Tmodel <- giorgi.tdph.maxim(x,
                                  theta0,
                                  nTD,
                                  nPH,
                                  event,
                                  ehazard,
                                  ehazardInt,
                                  P,
                                  p,
                                  k,
                                  nrowx,
                                  IntGL,
                                  nghq,
                                  cpti,
                                  cptj,
                                  cst,
                                  boundmin,
                                  boundmax,
                                  Int.chazard,
                                  Int.FDbase,
                                  Int.SDbase,
                                  int,
                                  control)
      }

    list(
      coefficients = Fmodel$coefficients,
      var = Fmodel$var,
      loglik = c(Tmodel$loglik, Fmodel$loglik),
      loglik.test = (-2 * (Tmodel$loglik - Fmodel$loglik)),
      iterations = Fmodel$iterations,
      cov.test = cov.test,
      cov.df = (4 * abs(nPH.Full - nPH)),
      message = Fmodel$message,
      convergence = Fmodel$convergence,
      p = p,
      nTD = nTD.Full,
      nPH = nPH.Full,
      nalpha = nalpha
    )
  }
  else{
    list(
      coefficients = Fmodel$coefficients,
      var = Fmodel$var,
      loglik = Fmodel$loglik,
      iterations = Fmodel$iterations,
      cov.test = FALSE,
      message = Fmodel$message,
      convergence = Fmodel$convergence,
      nTD = nTD,
      nPH = nPH,
      nalpha = nalpha,
      p = p
    )
  }

}
