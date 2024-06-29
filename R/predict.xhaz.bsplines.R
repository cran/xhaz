#' @title Predictions of excess hazard and net Survival from a \code{bsplines}
#'  object
#'
#' @description Function to predict excess hazard and net survival based on
#'  an object of class \code{bsplines}. The function allows the
#'  predictions at several time points but not exceeding the maximum time of
#'  follow-up from the baseline model.
#'
#'
#' @param object an object of class \code{bsplines}
#'
#' @param new.data new.data where is covariates
#'
#' @param times.pts time in year scale to calculate the excess hazard. The
#' default value is NULL. In this case, time variable must be provided in the
#' new.data
#'
#' @param baseline default is survival baseline; put \code{baseline = FALSE}
#'  to estimate the net survival with covariates
#'
#' @param ... additional arguments affecting the predictions of excess hazard
#' and net survival
#'
#' @keywords predict.bsplines
#'
#' @return An object of class predxhaz, which is a list of data.frame. Each
#' element of the list contains the estimates of hazard and survival at a fixed
#'  time point. The return of this function can be used to produce graphics of
#'  excess hazard or net survival, when times.pts argument is provided. This
#'  object contains:
#'
#' \item{times.pts}{the times value in year at which the excess hazard
#'   and or the net survival have been estimated}
#'
#' \item{hazard}{the excess hazard values based on the model of interest}
#'
#' \item{survival}{the net survival values based on the model of interest}
#'
#'
#'
#' @author Juste Goungounga, Robert Darlin Mba, Nathalie Graff\'eo and Roch Giorgi
#'
#' @references Goungounga JA, Touraine C, Graff\'eo N, Giorgi R;
#' CENSUR working survival group. Correcting for misclassification
#' and selection effects in estimating net survival in clinical trials.
#' BMC Med Res Methodol. 2019 May 16;19(1):104.
#' doi: 10.1186/s12874-019-0747-3. PMID: 31096911; PMCID: PMC6524224.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31096911/}{PubMed})
#'
#' Touraine C, Graff\'eo N, Giorgi R; CENSUR working survival group.
#' More accurate cancer-related excess mortality through correcting
#' background mortality for extra variables.
#' Stat Methods Med Res. 2020 Jan;29(1):122-136.
#' doi: 10.1177/0962280218823234. Epub 2019 Jan 23. PMID: 30674229.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/30674229/}{PubMed})
#'
#' Mba RD, Goungounga JA, Graff\'eo N, Giorgi R; CENSUR working survival group.
#' Correcting inaccurate background mortality in excess hazard models
#' through breakpoints. BMC Med Res Methodol. 2020 Oct 29;20(1):268.
#' doi: 10.1186/s12874-020-01139-z. PMID: 33121436; PMCID: PMC7596976.
#'   (\href{https://pubmed.ncbi.nlm.nih.gov/33121436/}{PubMed})
#'
#'
#'
#'
#' @seealso \code{\link{xhaz}}, \code{\link{print.bsplines}}, \code{\link{print.constant}}
#'
#' @examples
#'
#' \donttest{
#' library("survival")
#' library("numDeriv")
#' library("survexp.fr")
#' library("splines")
#' data("dataCancer", package = "xhaz")   # load the data set in the package
#'
#' fit.phBS <- xhaz(
#'         formula = Surv(obs_time_year, event) ~ ageCentre + immuno_trt,
#'         data = dataCancer, ratetable = survexp.fr,
#'         interval = c(0, NA, NA, max(dataCancer$obs_time_year)),
#'         rmap = list(age = 'age', sex = 'sexx', year = 'year_date'),
#'         baseline = "bsplines", pophaz = "classic")
#'
#'
#' print(fit.phBS)
#'
#'
#' predicted <- predict(object = fit.phBS,
#'                      new.data = dataCancer[1:10,],
#'                      times.pts = c(seq(0,10,1)),
#'                      baseline = TRUE)
#'
#'
#' #a list of predicted hazard and survival at different time points
#' print(predicted)
#'
#'
#' #predicted hazard and survival at time points 10 years
#' print(predicted[[10]])
#' }
#' @export

predict.bsplines <- function(object,
                           new.data = NULL,
                           times.pts = NULL,
                           baseline = TRUE,
                           ...) {
  Call <- match.call()
  if (any(object$bsplines) == TRUE)
    stop("Predict.bplines is not yet implemented for non-proportional hazards setting\n")

  int <- (object$interval)
  if (inherits(object, "bsplines")) {
    coeffBS <- object$coefficients[1:5]
    if (is.null(new.data)) {
      new.data <- object$data
      xx <- as.data.frame(
        model.matrix((object$terms),
                     (new.data)))[, -1, drop = FALSE]
      if (is.null(times.pts)) {
        times.pts_init <- times.pts
        m <- eval(object$terms, sys.parent())
        m[[1]] <- as.name("model.frame")
        times.pts <- object$data[,toString(as.name(m[[2]][[2]]))]
      }else{
        times.pts_init <- times.pts
        m <- eval(object$terms, sys.parent())
        m[[1]] <- as.name("model.frame")
        index_time <- which(colnames(object$data) %in% c(toString(as.name(m[[2]][[2]]))))
        time_name <- colnames(object$data)[index_time]
        time_data <- data.frame(times.pts)

        new.data <- lapply(1:nrow(time_data),
                           function(i){
                             my_colnames <- c(colnames(new.data), time_name)
                             my_new.data <- data.frame(cbind(new.data,
                                                             rep(time_data[i,],
                                                                 nrow(new.data))))
                             colnames(my_new.data) <- my_colnames
                             return(my_new.data)
                           })
      }
    }else{
      if (is.null(times.pts)) {
        times.pts_init <- times.pts
        m <- eval(object$terms, sys.parent())
        m[[1]] <- as.name("model.frame")
        times.pts <- try(new.data[,toString(as.name(m[[2]][[2]]))], TRUE)
        if (inherits(times.pts, "try-error"))
          stop("Need to provides time variable in the new.data or in the time.pts parameter")

        index_event <- which(colnames(object$data) %in% c(toString(as.name(m[[2]][[3]]))))
        index_time <- which(colnames(object$data) %in% c(toString(as.name(m[[2]][[2]]))))

        event_name <- colnames(object$data)[index_event]
        event_data <- data.frame(rep(0, nrow(new.data)))
        colnames(event_data) <- event_name
        time_name <- colnames(object$data)[index_time]
        time_data <- data.frame(times.pts)
        colnames(time_data) <- time_name
        new.data <- data.frame(cbind(new.data, event_data, time_data))
        xx <- as.data.frame(model.matrix((object$terms),(new.data)))[, -1, drop = FALSE]

      }else {
        times.pts_init <- times.pts
        m <- eval(object$terms, sys.parent())
        m[[1]] <- as.name("model.frame")
        index_event <- which(colnames(object$data) %in% c(toString(as.name(m[[2]][[3]]))))
        index_time <- which(colnames(object$data) %in% c(toString(as.name(m[[2]][[2]]))))

        event_name <- colnames(object$data)[index_event]
        event_data <- data.frame(rep(0, nrow(new.data)))
        colnames(event_data) <- event_name
        time_name <- colnames(object$data)[index_time]
        time_data <- data.frame(times.pts)
        colnames(time_data) <- time_name
        new.data <- lapply(1:nrow(time_data),
                           function(i){
                             my_colnames <- c(colnames(new.data),event_name, time_name)
                             my_new.data <- data.frame(cbind(new.data,
                                                             event_data,
                                                             rep(time_data[i,],
                                                                 nrow(new.data))))
                             colnames(my_new.data) <- my_colnames
                             return(my_new.data)
                             })
        xx <- as.data.frame(model.matrix((object$terms),
                                         (new.data[[1]])))[, -1, drop = FALSE]

      }
    }

    if (object$pophaz == "classic") {
      nalpha <- 0

    } else if (object$pophaz == "rescaled" |
               object$pophaz == "corrected") {
      indxAlpha <- which(stringr::str_detect(names(object$coefficients),
                                             pattern = "alpha"))
      nalpha <- length(indxAlpha)
    }



    nPH <- object$nPH
    nTD <- object$nTD

    coeffPred <- object$coefficients[(5 + 5*nTD + 1):(5 + 5*nTD + nPH )]
    k <- 3
    knot <- c(int[2], int[3])
    delta <- sort(c(rep(c(int[1], int[4]), k), knot))
    object$linear.predictors <- exp(as.matrix(xx) %*% coeffPred)
    rrBetaZ <- t(object$linear.predictors)
    CMUint <- list()

    for (i in 1:length(times.pts)) {
      CMUint[[i]] <- exp(-(
        integrate(function(times.pts, coeffBS)
          (exp(
            apply((coeffBS) * t(
              splines::spline.des(knots = delta, x = times.pts, ord =  k)$design
            ), 2, sum)
          )), 0, times.pts[i], coeffBS)$value
      ))

    }

    CHBSplines <-
      exp(apply((coeffBS) * t(
        splines::splineDesign(
          knots = delta,
          x = times.pts,
          ord = k,
          outer.ok = FALSE
        )
      ), 2, sum))
    if (is.null(times.pts_init)) {
      mypred <- suppressWarnings(round(data.frame(times.pts = times.pts,
                                                        hazard = CHBSplines,
                                                        survival = unlist(CMUint)), 4))


      class(mypred) <- c("data.frame","predxhaz")
    } else{
      mypred <- lapply(1:nrow(time_data),
                       function(i){
                         suppressWarnings(
                           round(data.frame(times.pts = rep(times.pts[[i]], nrow(new.data[[i]])),
                                            hazard = rep(CHBSplines[i], nrow(new.data[[i]])),
                                            survival = rep(CMUint[[i]], nrow(new.data[[i]]))), 4))
                       })

      class(mypred) <- c("list","predxhaz")
    }

    attributes(mypred)$call <- Call
    attributes(mypred)$baseline <- object$baseline
    attributes(mypred)$pophaz <- object$pophaz
    attributes(mypred)$coefficients <- object$coefficients
    attributes(mypred)$intervall <- object$interval

if (max(times.pts) > max(attr(mypred, "interval"))) {
  stop( "time must be inferior or equal to max value in interval specified to estimate the model parameter")
}

    if (baseline)  {
      return(mypred)
    } else{
      if (is.null(times.pts_init)) {
        mypred$hazard <- c(mypred$hazard * rrBetaZ)
        mypred$survival <- c(mypred$survival^rrBetaZ)
      }else{
        for (i in 1:length(times.pts)) {
          mypred[[i]]$hazard <- c(mypred[[i]]$hazard * rrBetaZ)
          mypred[[i]]$survival <- c(mypred[[i]]$survival^rrBetaZ)
        }

      }
      return(mypred)
    }
  }
  invisible()
}
