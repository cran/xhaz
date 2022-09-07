#' @import survival


survSplit2 <-
  function(formula, data, subset, na.action = na.pass, cut, start = "tstart",
           id, zero = 0, episode, end = "tstop",
           event = "event") {

  Call <- match.call()
  if (missing(formula) || is.data.frame(formula)) {
    if (missing(data)) {
      if (!missing(formula)) {
        names(Call)[[2]] <- "data"
        data <- formula
      }
      else stop("a data frame is required")
    }
    if (missing(end) || missing(event))
      stop("either a formula or the end and event arguments are required")
    if (!(is.character(event) && length(event) == 1 && event %in%
          names(data)))
      stop("'event' must be a variable name in the data set")
    if (!(is.character(end) && length(end) == 1 && end %in%
          names(data)))
      stop("'end' must be a variable name in the data set")
    if (!(is.character(start) && length(start) == 1))
      stop("'start' must be a variable name")
    if (start %in% names(data))
      temp <- paste(start, end, event, sep = ",")
    else temp <- paste(end, event, sep = ",")
    formula <- as.formula(paste("Surv(", temp, ")~ ."))
  }
  else if (missing(formula))
    stop("either a formula or the end and event arguments are required")
  indx <- match(c("data", "weights", "subset"),
                names(Call), nomatch = 0)
  temp <- Call[c(1L, indx)]
  temp$formula <- formula
  temp$na.action <- na.action
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)
  Y <- model.response(mf)
  states <- attr(Y, "states")
  if (!is.Surv(Y))
    stop("the model must have a Surv object as the response")
  if (!(attr(Y, "type") %in% c("right", "mright",
                               "counting", "mcounting")))
    stop(paste("not valid for", attr(Y, "type"),
               "censored survival data"))
  nY <- ncol(Y)
  ymiss <- is.na(Y)
  if (nY == 2) {
    if (any(Y[!ymiss, 1] <= zero))
      stop("'zero' parameter must be less than any observed times")
    Y <- cbind(zero, Y)
  }
  temp <- (Y[!ymiss, 1] >= Y[!ymiss, 2])
  if (any(temp))
    stop("start time must be < stop time")
  if (!is.numeric(cut) || any(!is.finite(cut)))
    stop("cut must be a vector of finite numbers")
  cut <- unique(sort(cut))
  ntimes <- length(cut)
  n <- nrow(data)
  if (!missing(id)) {
    if (!is.character(id))
      stop("id must be a variable name")
    if (id %in% names(mf))
      stop("the suggested id name is already present")
    id <- make.names(id)
    if (id %in% names(mf))
      stop("the suggested id name is already present")
    mf[[id]] <- 1:nrow(mf)
  }

  storage.mode(Y) <- "double"
  Csurvsplit <-  Csurvsplit2()
  index <- .Call(Csurvsplit, mf$tage, Y[, 2], as.double(cut))

  newdata <- mf[index$row, -1, drop = FALSE]
  row.names(newdata) <- NULL
  attr(newdata, "terms") <- NULL
  status <- Y[index$row, 3]
  status[index$censor] <- 0
  if (!is.null(states))
    status <- factor(status, labels = c("censor", states))
  if (inherits(formula[[2]], "call") && formula[[2]][[1]] ==
      as.name("Surv")) {
    temp <- match.call(Surv, formula[[2]])
    if (nY == 2) {
      if (missing(end) && !is.null(temp[["time"]]) &&
          is.name(temp[["time"]]))
        end <- as.character(temp[["time"]])
      if (missing(event) && !is.null(temp$time2) && is.name(temp$time2))
        event <- as.character(temp$time2)
      if (missing(event) && !is.null(temp$event) && is.name(temp$event))
        event <- as.character(temp$event)
    }
    else {
      if (missing(end) && !is.null(temp[["time"]]) &&
          is.name(temp["time"]))
        start <- as.character(temp[["time"]])
      if (missing(end) && !is.null(temp$time2) && is.name(temp$time2))
        end <- as.character(temp$time2)
      if (missing(event) && !is.null(temp$event) && is.name(temp$event))
        event <- as.character(temp$event)
      if (missing(start) && !is.null(temp$time) && is.name(temp$time))
        start <- as.character(temp$time)
    }
    newdata[[start]] <- index$start
    newdata[[end]] <- index$end
    newdata[[event]] <- status
  }
  else {
    if (inherits(formula[[2]], "name") == FALSE)
      stop("left hand side not recognized")
    temp <- as.character(formula[[2]])
    newdata[temp] <- Surv(index$start, index$end, status)
  }

  if (!missing(episode)) {
    if (!is.character(episode))
      stop("episode must be a character string")

    newdata[[make.names(episode)]] <- index$interval + 1
  }
  newdata
}
