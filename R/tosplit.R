#' @import survival

tosplit <- function(formula = formula,
                    add.rmap.cut = add.rmap.cut,
                    data = data, na.action, rmap, interval, subset) {

  Call <- match.call()
  m <- match.call(expand.dots = FALSE)

  # indx <- match(c("formula", "data", "subset", "na.action"),
  #               names(Call), nomatch = 0)

  indx <- match(c("formula", "data", "na.action"),
                names(Call), nomatch = 0)

  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1L]] <- as.name("model.frame")
  special <- c("strata")
  Terms <- if (missing(data)) {
    terms(formula, special)
  }
  else{
    terms(formula, special, data = data)
  }
  temp$formula <- Terms

  m <- eval(temp, sys.parent())

  if (missing(na.action)) {
    na.action <- NULL
  } else if (length(attr(m, "na.action"))) {
    temp$na.action <- na.pass
    m <- eval(temp, sys.parent())
  }


  if (missing(data)) {
    stop("Missing data data frame in which to interpret
           the variables named in the formula.")
  } else{
    if (is.na(match(rmap$age, names(data))))
      stop("Must have informations for age on the data set.")
    if (is.na(match(rmap$sex, names(data))))
      stop("Must have informations for sex on the data set.")
    if (is.na(match(rmap$year, names(data))))
      stop("Must have informations for date on the data set.")
  }





  myvarnames <- colnames(model.matrix(Terms, m)[,-1, drop = FALSE])

  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object.")

  attr(Terms, "intercept") <- 1

  type <- attr(Y, "type")

  if (ncol(Y) == 2) {
    time <- Y[, 1]
    event <- Y[, 2]
  } else{
    time <- Y[, 2] - Y[, 1]
    event <- Y[, 3]
  }
  event[time > max(interval, na.rm = TRUE)] <- 0
  time[time > max(interval, na.rm = TRUE)] <- max(interval, na.rm = TRUE)

  data$tage2 <- data$tage <- data$ageDiag <- ageDiag <- data[, rmap$age]
  data$tageDC <- data$ageDC <- ageDC <- ageDiag + time
  data2 <- data

  data2$id <- 1:nrow(data2)

  data2$time_2 <- data2$time_old <- data2$time <- time

  data2$tageDC <- data2$ageDC <- data2$age + data2$time



  tdata2 <- survSplit2(Surv(tage2 + time_2, event == 1) ~ ., data2,
                       cut = add.rmap.cut$cut,
                       episode = "break_interval")

  tdata2$time <- with(tdata2, c(tstop - tstart))
  colnames(tdata2)[which(colnames(tdata2) == "time")] <- toString(Terms[[2]][[2]])
  colnames(tdata2)[which(colnames(tdata2) == "event")] <- toString(Terms[[2]][[3]])

  return(list(tdata2 = tdata2, Call = Call))
}
