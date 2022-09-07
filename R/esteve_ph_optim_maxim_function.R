
esteve.ph.optim.maxim <- function(x, y,
                                  theta0,
                                  nvar,
                                  k,
                                  indic,
                                  event,
                                  integr1,
                                  ehazard,
                                  ehazardInt,
                                  control.iter.max,
                                  control.eps,
                                  Terms,
                                  strats,
                                  nstrata,
                                  add.rmap,
                                  add.rmap.cut,
                                  ageDiag,
                                  trace = 0,
                                  speedy = FALSE)
{
  ageDC <- ageDiag + y[, 1]
  if (add.rmap.cut$breakpoint == FALSE) {

    nvarStrata <- (nvar + k * nstrata + 1)
    FD <- 1
    no <- length(add.rmap)
    if (is.null(strats)) {
      strats <- rep(1, nrow(x))
    }
    if (!is.null(add.rmap)) {
      nalpha <- nlevels(add.rmap)
    } else{
      nalpha = 0
    }

    f <- function(theta0) {
      ll <- FDbeta <- SDbeta <- 0
      if (nvar > 0) {
        rr <- exp(rowSums(t(theta0[1:nvar] * t(x))))
      } else{
        rr <- 1
      }
      if (nalpha) {
        levels(add.rmap) <- 1:nlevels(add.rmap)
        if (length(levels(add.rmap)) < 2) {
          alpha0 <- exp(theta0[nvarStrata])
          alpha <- exp(theta0[nvarStrata])
        } else{
          idxtheta0 <- (nvarStrata):(nvar + k * nstrata + nalpha)
          alpha0 <- exp(theta0[idxtheta0])
          Madd.rmap <- model.matrix( ~ add.rmap - 1)


          alpha <- (Madd.rmap %*% alpha0)
        }
      }else{
        alpha <- 1
      }

      for (h in 1:nstrata) {
        indextau <- (nvar + (h - 1) * k + 1):(nvar + h * k)
        tau <- t(exp(theta0[indextau]) * t(indic * (strats == h)))
        tauInt <- t(exp(theta0[indextau]) * t(integr1 * (strats == h)))
        chazard <- colSums(as.matrix(t(tau))) * rr
        ohazard <- chazard + alpha * ehazard
        eHazard <- alpha * ehazardInt
        tll <- colSums(as.matrix(-rr * rowSums(tauInt)
                                 - eHazard +
                                   event * log(ohazard)))
        ll <- tll + ll
      }
      return(ll)
    }


    gradient <- function(theta0) {
      FDbeta <- 0
      FDtau <- FDalpha <- NULL
      if (nvar > 0) {
        rr <- exp(rowSums(t(theta0[1:nvar] * t(x))))
      } else{
        rr <- 1
      }
      if (nalpha) {
        levels(add.rmap) <- 1:nlevels(add.rmap)
        if (length(levels(add.rmap)) < 2) {
          alpha0 <- exp(theta0[nvarStrata])
          alpha <- exp(theta0[nvarStrata])
        } else{
          idxtheta0 <- (nvarStrata):(nvar + k * nstrata + nalpha)
          alpha0 <- exp(theta0[idxtheta0])

          Madd.rmap <- model.matrix( ~ add.rmap - 1)
          alpha <- (Madd.rmap %*% alpha0)
        }
      } else{
        alpha <- 1
      }



      for (h in 1:nstrata) {
        indextau <- (nvar + (h - 1) * k + 1):(nvar + h * k)
        tau <- t(exp(theta0[indextau]) * t(indic * (strats == h)))
        tauInt <- t(exp(theta0[indextau]) * t(integr1 * (strats == h)))
        chazard <- colSums(as.matrix(t(tau))) * rr
        ohazard <- chazard + alpha * ehazard
        tFDtau <- colSums(
          as.matrix(-rr * tauInt + (event * tau * rr) / c(ohazard)))
        if (nstrata == 1) {
          FDtau <- tFDtau
        }
        else{
          FDtau <- c(FDtau, tFDtau)
        }

        if (nalpha) {
          ohazard <- chazard + alpha * ehazard

          if (length(levels(add.rmap)) < 2) {
            eHazard <- alpha * ehazardInt

            tFDalpha <- colSums(as.matrix(-eHazard +
                                            event * alpha * c(ehazard) / c(ohazard)))
            FDalpha <- tFDalpha
          } else{
            Malpha <- t(alpha0 * t(Madd.rmap))
            tFDalpha <- colSums(as.matrix(-Malpha * c(ehazardInt) +
                                            event * Malpha * c(ehazard) / c(ohazard)))
            FDalpha <- tFDalpha
          }
        }


        if (nvar != 0) {
          ohazard <- chazard + alpha * ehazard
          tFDbeta <- -x * rr * rowSums(tauInt) + event * x * chazard / c(ohazard)
          if (nstrata == 1) {
            FDbeta <- tFDbeta
          }
          else{
            FDbeta <- FDbeta + tFDbeta
          }
        }
        else{
          if (nalpha) {
            FD <- c(FDtau, FDalpha)
          }else{
            FD <- c(FDtau)
          }
        }
      }
      if (nvar != 0) {
        if (nalpha) {
          FD <- c(colSums(as.matrix(FDbeta)), FDtau, FDalpha)
        } else{
          FD <- c(colSums(as.matrix(FDbeta)), FDtau)
        }
      }
      FD
    }

    if (is.null(trace)) {
      trace <- 0
    }


    if (nalpha) {
      nalpha <- 0
      theta0[1:(nvar + k * nstrata)] <- optim(par = theta0[1:(nvar + k * nstrata)],
                                              fn = f,
                                              gr = gradient,
                                              method = "L-BFGS-B",
                                              control = list(REPORT = 1,
                                                             maxit = 500,
                                                             fnscale = -1,
                                                             trace = trace))$par
      if (length(levels(add.rmap)) < 2) {
        nalpha <- 1
      } else{
        nalpha <- nlevels(add.rmap)
      }
    }


    if(speedy) {
      max_cores <- detectCores()
      used_cores <- min(max_cores, (max_cores - 2))
      if (used_cores < 2) {
        stop("We didn't detect enough cores for speedy")
      }
      cl <- makeCluster(used_cores)
      setDefaultCluster(cl = cl)
      res <- optimParallel(
        par = theta0,
        fn = f,
        gr = gradient,
        hessian = TRUE,
        method = "L-BFGS-B",
        control = list(REPORT = 1,
                       maxit = 1000,
                       fnscale = -1,
                       trace = trace),
        parallel = list(loginfo = TRUE)
      )
      setDefaultCluster(cl = NULL)
      stopCluster(cl)


    }else{
      res <- optim(par = theta0,
                   fn = f,
                   gr = gradient,
                   method = "L-BFGS-B",
                   hessian = TRUE,
                   control = list(REPORT = 1,
                                  maxit = 1000,
                                  fnscale = -1,
                                  trace = trace))
    }

    ll <- res$value
    theta0 <- res$par
    FD <- -gradient(theta0)
    SD <- -res$hessian
    iter <- res$counts[2]
    message <- res$message
    convergence <- res$convergence

    return(
      list(
        theta0 = theta0,
        ll = ll,
        FD = FD,
        SD = SD,
        iter = iter,
        convergence = convergence,
        message = message
      )
    )


}
else{
  nvarStrata <- (nvar + k * nstrata + 1)
  FD <- 1
  no <- length(add.rmap)
  if (is.null(strats)) {
    strats <- rep(1, nrow(x))
  }
  if (!is.null(add.rmap)) {
    nalpha <- (length(add.rmap.cut$cut) + 1)*nlevels(add.rmap)
  } else{
    nalpha <- 0
  }

  f <- function(theta0) {
    ll <- FDbeta <- SDbeta <- 0
    if (nvar > 0) {
      rr <- exp(rowSums(t(theta0[1:nvar] * t(x))))
    } else{
      rr <- 1
    }
    if (nalpha) {
      levels(add.rmap) <- 1:nlevels(add.rmap)
      if (length(levels(add.rmap)) < 2) {
        alpha0 <- exp(theta0[nvarStrata])
        alpha <- exp(theta0[nvarStrata])
      } else{
        idxtheta0 <- (nvarStrata):(nvar + k * nstrata + (length(add.rmap.cut$cut) + 1)*nlevels(add.rmap))
        alpha0 <- exp(theta0[idxtheta0])

        interval <- sort(c( add.rmap.cut$cut, max(ageDC) + 1, min(ageDC) - 1))
        nbcut <- length(add.rmap.cut$cut) + 1

colnames_Madd.rmap <- c(unlist(
  lapply(1:nlevels(add.rmap),
         function(i)
           paste(
             paste0("add.rmap",1:nlevels(add.rmap))[i],
             1:nbcut, sep = "_"))))

ageDCgroup <- cut(ageDC, breaks = c(0, add.rmap.cut$cut, c(max(ageDC) + 1)))


Madd.rmap <-
  do.call("cbind", lapply(1:ncol(model.matrix( ~ add.rmap + 0)), function(i)
    c((model.matrix( ~ add.rmap + 0))[, i]) * (model.matrix( ~ ageDCgroup - 1))))



colnames(Madd.rmap) <- c(colnames_Madd.rmap)



  alpha <- (Madd.rmap %*% alpha0)


      }
    }else{
      alpha <- 1
    }

    for (h in 1:nstrata) {
      indextau <- (nvar + (h - 1) * k + 1):(nvar + h * k)
      tau <- t(exp(theta0[indextau]) * t(indic * (strats == h)))
      tauInt <- t(exp(theta0[indextau]) * t(integr1 * (strats == h)))
      chazard <- colSums(as.matrix(t(tau))) * rr
      ohazard <- chazard + alpha * ehazard
      eHazard <- alpha * ehazardInt
      tll <- colSums(as.matrix(-rr * rowSums(tauInt)
                               - eHazard +
                                 event * log(ohazard)))
      ll <- tll + ll
    }
    return(ll)
  }


  gradient <- function(theta0) {

    FDbeta <- 0
    FDtau <- FDalpha <- NULL
    if (nvar > 0) {
      rr <- exp(rowSums(t(theta0[1:nvar] * t(x))))
    } else{
      rr <- 1
    }
    if (nalpha) {
      levels(add.rmap) <- 1:nlevels(add.rmap)
      if (length(levels(add.rmap)) < 2) {
        alpha0 <- exp(theta0[nvarStrata])
        alpha <- exp(theta0[nvarStrata])
      } else{
        idxtheta0 <- (nvarStrata):(nvar + k * nstrata + (length(add.rmap.cut$cut) + 1)*nlevels(add.rmap))

        alpha0 <- exp(theta0[idxtheta0])
        interval <- sort(c(add.rmap.cut$cut, max(ageDC) + 1, min(ageDC) - 1))
        nbcut <- length(add.rmap.cut$cut) + 1
        colnames_Madd.rmap <- c(unlist(
          lapply(1:nlevels(add.rmap),
                 function(i)
                   paste(
                     paste0("add.rmap",1:nlevels(add.rmap))[i],
                     1:nbcut, sep = "_"))))

        ageDCgroup <- cut(ageDC, breaks = c(0, add.rmap.cut$cut, c(max(ageDC) + 1)))


        Madd.rmap <-
          do.call("cbind", lapply(1:ncol(model.matrix( ~ add.rmap + 0)), function(i)
            c((model.matrix( ~ add.rmap + 0))[, i]) * (model.matrix( ~ ageDCgroup - 1))))



        colnames(Madd.rmap) <- c(colnames_Madd.rmap)

        alpha <- (Madd.rmap %*% alpha0)

      }
    } else{
      alpha <- 1
    }



    for (h in 1:nstrata) {
      indextau <- (nvar + (h - 1) * k + 1):(nvar + h * k)
      tau <- t(exp(theta0[indextau]) * t(indic * (strats == h)))
      tauInt <- t(exp(theta0[indextau]) * t(integr1 * (strats == h)))
      chazard <- colSums(as.matrix(t(tau))) * rr
      ohazard <- chazard + alpha * ehazard
      tFDtau <- colSums(
        as.matrix(-rr * tauInt + (event * tau * rr) / c(ohazard)))
      if (nstrata == 1) {
        FDtau <- tFDtau
      }
      else{
        FDtau <- c(FDtau, tFDtau)
      }

      if (nalpha) {
        ohazard <- chazard + alpha * ehazard

        if (length(levels(add.rmap)) < 2) {
          eHazard <- alpha * ehazardInt

          tFDalpha <- colSums(as.matrix(-eHazard +
                                event * alpha * c(ehazard) / c(ohazard)))
          FDalpha <- tFDalpha
        } else{
          Malpha <- t(alpha0 * t(Madd.rmap))
          tFDalpha <- colSums(as.matrix(-Malpha * c(ehazardInt) +
                                event * Malpha * c(ehazard) / c(ohazard)))
          FDalpha <- tFDalpha
        }
      }


      if (nvar != 0) {
        ohazard <- chazard + alpha * ehazard
        tFDbeta <- -x * rr * rowSums(tauInt) + event * x * chazard / c(ohazard)
        if (nstrata == 1) {
          FDbeta <- tFDbeta
        }
        else{
          FDbeta <- FDbeta + tFDbeta
        }
      }
      else{
        if (nalpha) {
          FD <- c(FDtau, FDalpha)
        }else{
          FD <- c(FDtau)
        }
      }
    }
    if (nvar != 0) {
      if (nalpha) {
        FD <- c(colSums(as.matrix(FDbeta)), FDtau, FDalpha)
      } else{
        FD <- c(colSums(as.matrix(FDbeta)), FDtau)
      }
    }
    FD
  }

  if (is.null(trace)) {
    trace <- 0
  }


  if (nalpha) {
    # nalpha <- 0
    # theta0 [1:(nvar + k * nstrata)] <- optim(par = theta0[1:(nvar + k * nstrata)],
    #                                         fn = f,
    #                                         gr = gradient,
    #                                         method = "L-BFGS-B",
    #                                         control = list(REPORT = 1,
    #                                                        maxit = 500,
    #                                                        fnscale = -1,
    #                                                        trace = trace))$par

    theta0  <- optim(
      par = theta0,
      fn = f,
      gr = gradient,
      method = "L-BFGS-B",
      control = list(
        REPORT = 1,
        maxit = 500,
        fnscale = -1,
        trace = trace
      )
    )$par

    if (length(levels(add.rmap)) < 2) {
      nalpha <- 1
    } else{
      nalpha <- nlevels(add.rmap)
    }
  }


  if(speedy) {
    max_cores <- detectCores()
    used_cores <- min(max_cores, (max_cores - 2))
    if (used_cores < 2) {
      stop("We didn't detect enough cores for speedy")
      }
  cl <- makeCluster(used_cores)
  setDefaultCluster(cl = cl)
  res <- optimParallel(
    par = theta0,
    fn = f,
    gr = gradient,
    hessian = TRUE,
    method = "L-BFGS-B",
    control = list(REPORT = 1,
                   maxit = 1000,
                   fnscale = -1,
                   trace = trace),
    parallel = list(loginfo = TRUE)
  )
  setDefaultCluster(cl = NULL)
  stopCluster(cl)


}else{
  res <- optim(par = theta0,
               fn = f,
               gr = gradient,
               method = "L-BFGS-B",
               hessian = TRUE,
               control = list(REPORT = 1,
                              maxit = 1000,
                              fnscale = -1,
                              trace = trace))
}

  ll <- res$value
  theta0 <- res$par
  FD <- -gradient(theta0)
  SD <- -res$hessian
  iter <- res$counts[2]
  message <- res$message
  convergence <- res$convergence

  return(
    list(
      theta0 = theta0,
      ll = ll,
      FD = FD,
      SD = SD,
      iter = iter,
      convergence = convergence,
      message = message
    )
  )
}
}
