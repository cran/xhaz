#' @title duplicate function
#'
#' @description  Duplicate data for survival analysis in the context of
#' competing risks, where an individual can experience only one of alternative
#' events, using the Lunn & McNeil (Biometrics, 1995) approaches.
#' Duplication of data proceeds as follows: Suppose that we study \code{J}
#' distinct types of events. Each observation concerning a given subject is
#' duplicated \code{J} times, with one row for each type of event. In addition,
#' \code{(J-1)} dummy variables are created, each indicating the type of event
#' in relation with that observation (\code{delta.j=1} if the event of type j
#' is the observed one and \code{0} otherwise).
#' Since, for a given subject, only the first occurring event is considered,
#' the status indicator equals \code{1} for that event and \code{0} for all the
#' others. In the case of a censored observation (dropout or administrative
#' censoring), the same principle applies also: duplication of each subject's
#' data is made \code{J} times with \code{(J-1)} dummy variables and a status
#' indicator equal to \code{0} for all observations.
#'
#' @param status the censoring status indicator (numeric vector),
#'  0=alive, 1=dead.
#'
#'
#' @param event the indicator of the event type (numeric vector).
#'  By default, the event==0 acts as the censoring indicator.
#'
#'
#' @param data a data frame containing the data to duplicate.
#'
#'
#' @keywords duplicate
#'
#' @return A data.frame containing the duplicated data with the new dummy
#' variables, named \code{delta.number_of_the_event}, indicating the type of
#' event.
#'
#' @author Roch Giorgi
#'
#' @references Lunn M and McNeil D. Applying Cox regression to competing risks.
#' Biometrics 1995;51:524-532
#' (\href{https://pubmed.ncbi.nlm.nih.gov/7662841/}{PubMed})
#'
#' @examples
#'
#' ## Create the simplest test data set
#' data1 <- data.frame(futime     = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
#'                     fustat     = c(0, 1, 1, 1, 0, 0, 1, 0, 1, 1),
#'                     firstevent = c(0, 2, 1, 2, 0, 0, 1, 0, 2, 2),
#'                     sex        = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0))
#'
#' ## Duplicate data1 with firstevent == 0 as the censoring indicator.
#' library(xhaz)
#' dupli.data <- duplicate(status=fustat, event=firstevent, data=data1)
#'
#'
#' data2 <- data.frame(futime = c(10, 2, 7, 3, 4, 9, 13, 2, 5, 9),
#'                     fustat = c(0, 1, 1, 1, 0, 0, 1, 0, 1, 1),
#'                     firstevent = c(3, 2, 1, 2, 3, 3, 1, 3, 2, 2),
#'                     sex = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0))
#'
#'
#'## Duplicate data1 with firstevent == 3 as the censoring indicator.
#'
#' dupli.data <- duplicate(status = fustat,
#'                         event = firstevent == 3,
#'                         data = data2)
#'
#'
#' # Joint modeling
#' coxph(Surv(futime, fustat) ~ delta.2 + sex + delta.2:(sex), data = dupli.data)
#'
#' coxph(Surv(futime, fustat) ~ delta.1 + sex + delta.1:(sex), data = dupli.data)
#'
#' # exemple using ccr.mevents data
#'
#'
#' ccr.mevents$loc.rec <- as.numeric(ccr.mevents$event == 1)
#' ccr.mevents$dist.rec <- as.numeric(ccr.mevents$event == 2)
#' ccr.mevents$death <- as.numeric(ccr.mevents$event == 3)
#' # Age centered to mean and scaled
#' ccr.mevents$agecr <-  scale(ccr.mevents$age, TRUE, TRUE)
#'
#' ## Duplication of the data with local recurrence as the reference
#' dupli.ccr.mevents <- duplicate(status = status,
#'                                event = event, data = ccr.mevents)
#' head(dupli.ccr.mevents)
#' # joint model including overall mortality modelling

#' fit <- coxph(Surv(time, status) ~ agecr + sexe + stage + delta.2 + delta.3,
#'              data = dupli.ccr.mevents)
#'
#' fit
#'
#' # add expected mortality from french life table to the data
#'
#'
#' library(survexp.fr)
#' fit.haz <- exphaz(formula = Surv(time, death) ~ 1,
#'                   data = dupli.ccr.mevents,
#'                   ratetable = survexp.fr, only_ehazard = TRUE,
#'                   rmap = list(age = 'age', sex = 'sexe', year = 'date_diag'))
#'
#'  dupli.ccr.mevents$mua <- fit.haz$ehazard * dupli.ccr.mevents$delta.3
#'
#' # joint model including excess hazard modelling
#' library(mexhaz)
#' fit.mort <- mexhaz(
#'     Surv(time, status) ~ delta.2 + delta.3,
#'     data = dupli.ccr.mevents, base = "exp.bs", degree = 3, knots = c(1),
#'     expected = "mua")
#'
#' fit.mort
#'
#'
#' @export

duplicate <- function(status, event, data){
  call <- match.call()
  status <- as.character(call[[2]])
  event <- as.character(call[[3]])
  ref.e <- 0
  if (length(call[[3]]) != 1) {
    event <- as.character(eval(expression(call[[3]][[2]])))
    ref.e <- call[[3]][[3]]
    }
  if (!is.numeric(data[, status]))
    stop("Status variable is not numeric.")
  if (!is.numeric(data[, event]))
    stop("Event variable is not numeric.")
  n.d <- length(unique(data[, event]))
  data <- cbind(data.frame(tid = c(1:nrow(data))), data)
  d.data <- data.frame(mapply(rep, data, n.d, SIMPLIFY = FALSE))
  d.data <- d.data[do.call(order, d.data), ]
  dimnames(d.data)[[1]] <- 1:nrow(d.data)
  temp.e <- sort(unique(data[, event]))
  temp.e <- c(ref.e, sort(temp.e[!temp.e == ref.e]))
  d.data[, event] <- rep(temp.e, nrow(d.data)/n.d)
  d.data <- cbind(d.data,
                  mapply(
                    rep,
                    data.frame(delta = contr.treatment(as.factor(temp.e))),
                    nrow(data)))
  d.data[, status] <- rep(0, nrow(d.data))
  d.data[match(paste(data$tid, data[, event]),
               paste(d.data$tid, d.data[, event])), status] <- data[, status]
  d.data <- d.data[d.data[,event] != ref.e, -1]
  dimnames(d.data)[[1]] <- 1:nrow(d.data)
  class(d.data) <- c("data.frame","duplicate")
  return(d.data)
}
