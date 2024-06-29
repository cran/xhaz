
breakpoint_with_cut <- function(formula = formula,
                                data = data,
                                ratetable = ratetable,
                                rmap = rmap,
                                baseline = baseline,
                                pophaz = pophaz,
                                only_ehazard = only_ehazard,
                                add.rmap = add.rmap,
                                add.rmap.cut = add.rmap.cut,
                                interval = interval,
                                ratedata = ratedata,
                                subset = subset,
                                na.action = na.action,
                                init = init,
                                control = control,
                                optim = optim,
                                scale = scale,
                                trace = trace,
                                speedy = speedy,
                                nghq = nghq, m_int = m_int, rcall, method,...) {
  time_elapsed0 <- as.numeric(base::proc.time()[3])


    tsplitting <- splitting <- TRUE
    tnewdata <- data
    newdata2 <- tosplit(formula = formula,
                        add.rmap.cut = add.rmap.cut,
                        data = data, rmap = rmap,
                        interval = interval, subset = subset)



    splitting <- FALSE
    data <- newdata2$tdata2


    fit <- xhaz2(formula = formula,
                        data = data,
                        ratetable = ratetable, rmap = rmap,
                        baseline  = baseline,
                        pophaz = pophaz,
                        only_ehazard = only_ehazard,
                        add.rmap = add.rmap,
                        add.rmap.cut = add.rmap.cut,
                        splitting  = splitting,
                        interval = interval,
                        ratedata = ratedata,
                        subset  = subset,
                        na.action = na.action,
                        init = init,
                        control = control,
                        optim = optim,
                        scale = scale,
                        trace = trace,
                        speedy = speedy,
                        nghq = nghq, m_int = m_int,
                        rcall = rcall, method = method, ...)


    fit$splitting <- tsplitting
    return(fit)


}
