

esteve.ph.maxim <- function(x, y, theta0, nvar, k, indic, event, integr1, ehazard,
                            ehazardInt, control.iter.max, control.eps, Terms,
                            strats, nstrata, add.rmap, add.rmap.cut, ageDiag, trace) {
  ageDC <- ageDiag + y[, 1]

  if (add.rmap.cut$breakpoint == FALSE) {

    iter <- 0
    FD <- 1
    dum <- 1
    if (!is.null(add.rmap)) {
      nalpha <- nlevels(add.rmap)
    } else{
      nalpha = 0
    }

    if (is.null(strats)) {
      strats <- rep(1, nrow(x))
    }
    while (max(abs(FD)) > control.eps) {
      ll <-  0
      SDP <- FDtau <- SDtau <- FDalpha <- SDalpha <- NULL
      FDbeta <- SDbeta <- SDTauAlpha <- SDPBetaAlpha <- SDPBetaTau <- NULL

      if (iter > control.iter.max)
        stop(paste("Ran out of iterations",
                   control.iter.max,
                   "and did not converge."))
      iter <- iter + 1

      if (nvar > 0) {
        rr <- exp(rowSums(t(theta0[1:nvar] * t(x))))
      } else{
        rr <- 1
      }

      if (nalpha) {
        Ind.alpha <- outer(add.rmap, levels(add.rmap), '==')
        nvarstrata <- (nvar + k * nstrata + 1):(nvar + k * nstrata + nalpha)
        alpha <- t(t(Ind.alpha) * as.vector(exp(theta0[nvarstrata])))
      } else{
        alpha <- 1
      }

      for (h in 1:nstrata) {
        tauid <- (nvar + (h - 1) * k + 1):(nvar + h * k)
        tau <- t(exp(theta0[tauid]) * t(indic * (strats == h)))
        tauInt <- t(exp(theta0[tauid]) * t(integr1 * (strats == h)))
        chazard <- rowSums(tau * rr)
        if (nalpha) {
          ohazard <- chazard + rowSums(as.matrix(alpha * ehazard))
          tll <- sum(-rr * rowSums(tauInt) -
                       rowSums(as.matrix(alpha * ehazardInt)) +
                       event * log(ohazard),
                     na.rm = TRUE)
        } else{
          ohazard <- chazard + ehazard

          tll <- sum(-rr * rowSums(tauInt) -
                       rowSums(as.matrix(ehazardInt)) +
                       event * log(ohazard),na.rm = TRUE)
        }

        ll <- tll + ll
        ohazard <- chazard + rowSums(as.matrix(alpha * ehazard))
        tFDtau <- colSums(-tauInt * rr +
                            (event * tau * rr) / c(ohazard),
                          na.rm = T)
        tSDtau <- diag(colSums((
          -tauInt * rr +
            event * tau * rr * (c(chazard + rowSums(
              as.matrix(alpha * ehazard)
            )) - tau * rr) / c(chazard + rowSums(as.matrix(alpha * ehazard))) ^ 2
        ), na.rm = T),
        k,
        k)

        FDtau <- c(FDtau, tFDtau)
        SDtau <- cbind(SDtau, tSDtau)

        if (nalpha) {
          tFDalpha <- colSums(-alpha * ehazardInt + event * (alpha * ehazard) /
                                c(chazard + rowSums(as.matrix(alpha * ehazard))),
                              na.rm = T)

          TSDalpha <- diag(colSums(-alpha * ehazardInt +
                                     event * alpha * ehazard *
                                     (chazard + rowSums(
                                       as.matrix(alpha * ehazard))
                                      - alpha * ehazard) /
                                     (c(chazard + rowSums(
                                       as.matrix(alpha * ehazard))) ^ 2),
                                   na.rm = T),
                           nalpha,
                           nalpha)

          SDPTauAlpha <- -t(event * tau * rr) %*%
            ((alpha * ehazard) /
               c(chazard + rowSums(as.matrix(alpha * ehazard))) ^ 2)

          SDPBetaAlpha <- -t(x) %*% (event * alpha * (ehazard * chazard) /
                                       (c(chazard +
                                            rowSums(as.matrix(alpha * ehazard))
                                       ) ^ 2))

          FDalpha <- c(FDalpha, tFDalpha)
          SDalpha <- cbind(SDalpha, TSDalpha)
        }

        if (nvar != 0) {
          tFDbeta <- -t(x) %*% (rr * rowSums(tauInt)) +
            t(event * x) %*% (chazard /
                                c(chazard + rowSums(as.matrix(alpha * ehazard))))
          tSDbeta <- -t(x) %*% (x * rr * rowSums(tauInt)) +
            t(event * x) %*% (x * chazard *
                                rowSums(as.matrix(alpha * ehazard)) /
                                c(chazard + rowSums(as.matrix(alpha * ehazard))) ^ 2)

          tSDPBetaTau <- -t(x) %*% (tauInt * rr) +
            t(x) %*% (tau * event * rr *
                        (rowSums(as.matrix(alpha * ehazard))) /
                        c(chazard + rowSums(as.matrix(alpha * ehazard))) ^ 2)

          FDbeta <- c(FDbeta, tFDbeta)
          SDbeta <- cbind(SDbeta, tSDbeta)
          SDPBetaTau <- cbind(SDPBetaTau, tSDPBetaTau)
        }
      }
      if (nvar != 0) {
        if (nalpha) {
          FD <- c(FDbeta, FDtau, FDalpha)
          SD <- -rbind( cbind(SDbeta, SDPBetaTau, SDPBetaAlpha),
                        cbind(t(SDPBetaTau), SDtau, SDPTauAlpha),
                        cbind(t(SDPBetaAlpha), t(SDPTauAlpha), SDalpha))
          colnames(SD) <- names(theta0)
          rownames(SD) <- names(theta0)
        } else{
          FD <- c(FDbeta, FDtau)
          SD <- -rbind(cbind(SDbeta, SDPBetaTau), cbind(t(SDPBetaTau), SDtau))
          colnames(SD) <- names(theta0)
          rownames(SD) <- names(theta0)
        }
      } else{
        if (nalpha) {
          FD <- c(FDtau, FDalpha)
          SD <- rbind(cbind(SDtau, SDPTauAlpha), cbind(t(SDPTauAlpha), SDalpha))
          colnames(SD) <- names(theta0)
          rownames(SD) <- names(theta0)
        } else{
          FD <- FDtau
          SD <- -SDtau
          colnames(SD) <- names(theta0)
          rownames(SD) <- names(theta0)
        }
      }

      dum <- try(solve(qr(SD), FD))
      if (!is.numeric(dum)) {
        stop("Matrix not definite positive. Check for colinearity in the data set.")
        return(list(
          theta0 = theta0,
          ll = ll,
          FD = FD,
          SD = SD,
          iter = iter
        ))
      }

      theta0 <- theta0 + dum
      #diff.dum <- theta0- theta

      if (trace == TRUE) {
        cat("#### Iteration number:", iter, "####", "\n")
        cat("\n")
        cat("## Inverse of Hessian matrix (-H) at",
            paste('theta', iter, sep = '_'),
            "====",
            "\n")
        print(solve(SD))
        cat("\n")
        cat("## Gradiant at", paste('theta', iter, sep = '_'), "====", "\n")
        print(FD)
        cat("\n")
        cat("## Estimate values of",
            paste('theta', iter, sep = '_'),
            "====",
            "\n")
        print(theta0)
        cat("\n")
      }

      # If 'covtest' are required.
      # Used for the score test which needs only 1 iteration.
      if (control.iter.max == 1) {
        return(list(
          theta0 = theta0,
          ll = ll,
          FD = FD,
          SD = SD,
          iter = iter
        ))
      }
      NULL
    }
    list(theta0 = theta0, ll = ll, FD = FD, SD = SD, iter = iter)

  }else {
      stop("Please use optim method for breakpoint model")
    }


}
