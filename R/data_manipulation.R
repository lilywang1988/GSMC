### Data manipulation ###
#' Trim the data
#'
#' Trim the data according to event number or time
#'
#' \code{data.trim} is to trim the data upto \code{t}, \code{data.trim.d} is to trim the data upto the cound \code{d}.
#'
#'
#' @param data There are two possible structures allowed for this input data. The first type needs to have \code{trimmed=F} and include variables: a \code{treatment} variable with "experimental" denoting treatment group, \code{cnsr} variable with value 1 denoting censoring, \code{ct} variable denoting event time from the origin of the study, which equals the sum of entering time \code{enterT} and the survival time (time to event or censoring). A dataset simulated from from R package \href{https://github.com/keaven/nphsim}{nphsim} should fit the first type well enough (see the example1). The second type can be any data.frame or data.table output from a \code{data.trim} function, including variables: \code{ct} denoting event time from the origin of the study or the sum of entering time and the survival time, \code{survival} denoting the survival time or time to event/censoring, \code{delta} as an event indicator, \code{enterT} as entering time (example 2). For the second type, we set \code{trimmed=T} to avoid extra computations, but should be fine if \code{trimmed=F}.
#' @param t Time of interest to pause/stop the study, which could be an interim stage or the final stage.
#' @param d Event counts to pause/stop the study.
#' @param trimmed Whether this data has been trimmed by \code{data.trim} or \code{data.trim.d} before.
#' @return Note that \code{data.trim} only outputs a data.table odered by \code{ct}, the event/censoring time since the start of the study (calendar scale), including variables in the input data.table/frame \code{data}, and additional/updated variables of event indicator \code{delta}, \code{ct}, follow-up time \code{survival} since the enrollment.
#' \code{data.trim.d} outpus a list of two components. The first component is the data censored with \code{d} events have been observed, ordered by \code{ct}, the event/censoring time since the start of the study (calendar scale). The second component is the time of the stopping point when \code{d} events have been observed.
#'
#' @seealso \code{\link{FH.frac.cal}}
#' @author Lili Wang
#'
data.trim <- function(t, data, trimmed = F){
  data.ord<- data[order(ct)]
  if(trimmed == F) 
    data.ord[, delta := 1-cnsr][, trt := ifelse(treatment == "experimental", 1, 0)]
  data.ord.t <- data.ord[enterT <= t] # only select the at-risk subjects
  data.ord.t[, delta := ifelse(ct > t, 0, delta)][,ct := pmin(ct, t)][,survival := ct - enterT]
  #print(paste0("t_time: ", interim_time, " and t: ", t))
  return(copy(data.ord.t))
}

# trim the data when there are enough deaths
#' @rdname data.trim
data.trim.d <- function(d, data, trimmed=F){
  data.ord <- data[order(ct)]
  if(trimmed == F) data.ord[, delta := 1 - cnsr][, trt := ifelse(treatment == "experimental", 1, 0)]

  row <- which.max(cumsum(data.ord$delta) == d)
  t<- as.numeric(data.ord[row, 'ct'])

  data.ord.t <- data.ord[enterT <= t] # only select the at-risk subjects
  data.ord.t[, delta := ifelse(ct > t, 0, delta)][,ct := pmin(ct, t)][, survival := ct - enterT]
  #print(paste0("t_time: ", interim_time, " and t: ", t))
  return(list(copy(data.ord.t), t))
}



