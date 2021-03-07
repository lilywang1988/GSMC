#' Basic function for Fleming-Harrington family weighted log-rank tests
#'
#' Basic function to build the table for the calculation of the Fleming-Harrington family of weighted log-rank tests.
#'
#'
#'
#' @param survival Time to event or censoring.
#' @param delta Event indicator: 1 for observed cases, 0 for censored cases.
#' @param trt Treatment assignment indicator: 1 for treatment group, 0 for control group.
#' @param rho,gamma Parameters for Fleming-Harrington family with \eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.
#'
#'
#' @return Build a table for Fleming-Harrington log-rank test which ouputs \emph{ordered} \code{survival} as follow-up times, \code{Surv} as predictable survival functions \eqn{S(t^-)}, \code{Surv.exact} as exact survival functions \eqn{S(t)}, \code{delta} as event indicators, \code{trt} as treatment assignement (treated=1, control=0), \code{weight} as weight functon calcualted from the predictable survival functions \code{Surv}.
#'
#' In addition, the output also include \code{O1} as the observed events from the treatment arm, \code{E1} as the expected events from the treatment arm, \code{Cov} as the estimated variance without considering the weights.
#' @author Lili Wang
#'
FH.table <- function(survival, delta, trt, rho, gamma){
  ord <- order(survival)
  survival <- survival[ord]
  delta <- delta[ord]
  trt <- trt[ord]
  n <- length(delta)
  if(n != length(survival)) stop("Unequal lengths of survival and delta")
  #delete the last delta=1 observation to avoid the situation of S=0
  if(delta[n] == 1) {
    survival = survival[ - n]
    delta = delta[ - n]
    trt = trt[ - n]
    n = n - 1
  }
  # Note the Surv is prepared for the weight, which is predictable, corresponding to survKM_minus.
  out <- data.table(
    survival = survival,
    Surv = c(1, cumprod(1 - delta / (n:1))[ - n]),
    Surv.exact = cumprod(1 - delta / (n:1)), 
    delta = delta, 
    trt = trt
    )
  Y <- n:1
  P1 <- rev(cumsum(rev(trt))) / Y
  P0 <- 1 - P1
  out[, weight := Surv^rho * (1 - Surv)^gamma][,O1 := trt*delta][,E1 := P1 * delta][,Cov := P1 * P0 * delta]
  return(out)
}

#' Basic function for standard log-rank test
#'
#' Build the table for log-rank test calculation.
#'
#'
#'
#' @inheritParams FH.table
#'
#' @return Build a table for log-rank test which ouputs \emph{ordered} \code{survival} as follow-up times, \code{delta} as event indicators,\code{trt} as treatment assignement (treated=1, control=0), \code{Y} as the at-risk numbers, \code{P1} as the proportion of treated set, \code{P0} as the proportion of the control set.
#'
#' In addition, the output also include \code{O1} as the observed events from the treatment arm, \code{E1} as the expected events from the treatment arm, \code{Cov} as the estimated variance.
#' @seealso \code{\link{FH.test}}, \code{\link{I_t.2}}, \code{\link{WLR.test.cov}},\code{\link{WLR.test.cor}}
#'
#'
#' @author Lili Wang
#'
#'
logrank.table <- function(survival, delta, trt){
  ord <- order(survival)
  survival <- survival[ord]
  delta <- delta[ord]
  trt <- trt[ord]
  n <- length(delta)
  if(n != length(survival)) stop("Unequal lengths of survival and delta")

  #delete the last delta=1 observation to avoid the situation of S=0
  if(delta[n] == 1) {
    survival = survival[ - n]
    delta = delta[ - n]
    trt = trt[ - n]
    n = n - 1
  }
  Y <- n:1
  P1 <- rev(cumsum(rev(trt))) / Y
  P0 <- 1 - P1
  O1 <- trt * delta
  E1 <- P1 * delta
  Cov <- P1 * P0 * delta
  out = data.table(O1, E1, Cov, P1, P0, Y, survival, delta, trt)
  return(out)
}

#' Calculate the survival functions
#'
#' Calculate the survival function, either the predictable one \eqn{S(t^-)} using \code{survKM_minus} or \eqn{S(t)} using \code{survKM_exact}.
#'
#' @param v Time vector to give the corresponding survival functions.
#' @param survival Input follow-up times.
#' @param delta Input event indicators.
#'
#' @return \code{survKM_minus} returns the predictable one \eqn{S(t^-)}, and \code{survKM_exact} returns \eqn{S(t)}.
#'
#' @author Lili Wang
survKM_minus <- function(v, survival, delta){
  v = as.vector(v)
  ord <- order(survival)
  survival <- survival[ord]
  delta <- delta[ord]
  n = length(delta)
  if(delta[n] == 1) {
    survival = survival[ - n]
    delta = delta[ - n]
    n = n - 1
  }
  if(is.unsorted(survival)) {
    ord = order(survival)
    survival = survival[ord]
    delta = delta[ord]
    }
  if(n != length(survival)) stop("Unequal lengths of survival and delta")
  c(1, cumprod(1 - delta / (n:1)))[
    sapply(v, function(tt) which.min(survival < tt))
    ] # corresponding to t^-, that is to make it predictable
}
#' @rdname survKM_minus
survKM_exact <- function(v, survival, delta){
  v <- as.vector(v)
  ord <- order(survival)
  survival <- survival[ord]
  delta <- delta[ord]
  n <- length(delta)
  if(delta[n] == 1) {
    survival = survival[ - n]
    delta = delta[ - n]
    n = n - 1
  }
  if(is.unsorted(survival)) {
    ord = order(survival)
    survival = survival[ord]
    delta = delta[ord];}
  if(n != length(survival)) stop("Unequal lengths of survival and delta")
  cumprod(1 - delta / (n:1))[sapply(v, function(tt) which.min(survival < tt))]
}

#'
#' Fleming-Harrington weighted log-rank tests
#'
#' Calculating the Fleming-Harrington weighted log-rank tests
#'
#'
#'
#' @param survival Time to event or censoring.
#' @param delta Event indicators.
#' @param trt Treatment assignment indicator with 1 denoting the treated group, and 0 denoting the placebo group.
#' @param rho First power parameter for the Fleming-Harrington weight which weighs on the early departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#' @param gamma Second power parameter for the Fleming-Harrington weight which weighs on the late departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#' @return A list 3 different components
#'   \item{O1}{Observed number of weighted events (with a multiplication of corresponding weights) in the treatment arm.}
#'   \item{E1 }{Expected number of weighted events (with a multiplication of corresponding weights) in the treatment arm.}
#'   \item{Z}{Weighted log-rank test statistic.}
#'
#' @author Lili Wang
#' @seealso \code{\link{WLR.test}}
#' @examples
#' \dontrun{
#' # Examples for FH.test and WLR.test
#' set.seed(12345)
#' data_temp<- nphsim(nsim=1,lambdaC=log(2)/6, lambdaE = c(log(2)/6,log(2)/6*0.7), ssC=250, intervals = c(2),ssE=250, gamma=500/14, R=14, eta=1e-5, fixEnrollTime = TRUE)$simd
#' data_final<-data.trim.d(100,data_temp)[[1]]
#' rho=1
#' gamma=0
#' # compare the 3 different ways below:
#' #library(survival)
#' sqrt(survdiff(Surv(survival,delta)~trt, data =data_final,rho=rho)$chisq)
#' FH.test(survival=data_final$survival,delta=data_final$delta,trt=data_final$trt,rho=rho,gamma=gamma)
#' WLR.test(survival=data_final$survival,delta=data_final$delta,trt=data_final$trt,w=function(...){survKM_minus(...)^rho*(1-survKM_minus(...))^gamma})
#' }
#'
FH.test <- function(survival, delta, trt, rho, gamma){
  FH_table <- FH.table(survival, delta, trt, rho, gamma)
  O1 <- sum(FH_table$weight * FH_table$O1)
  E1 <- sum(FH_table$weight * FH_table$E1)
  V <- sum(FH_table$weight^2 * FH_table$Cov)
  Z<-(O1 - E1) / sqrt(V)
  return(list(O1 = O1,E1 = E1,Z = Z))
}

#' Weighted log-rank tests with any input weight
#'
#' Weighted log-rank test for any input weight function.
#' @param survival Time to event or censoring.
#' @param delta Event indicator: 1 for observed cases, 0 for censored cases.
#' @param trt Treatment assignment indicator: 1 for treatment group, 0 for control
#' @param w Weight function, with default to be 1, which is similar to the use of input arbitray weight in \code{\link{WLR.test.cov}}. Please also refer to the examples as well.
#' @author Lili Wang
#' @seealso \code{\link{FH.test}}
#' @examples
#' \dontrun{
#' # Examples for FH.test and WLR.test
#' set.seed(12345)
#' data_temp<- nphsim(nsim=1,lambdaC=log(2)/6, lambdaE = c(log(2)/6,log(2)/6*0.7), ssC=250, intervals = c(2),ssE=250, gamma=500/14, R=14, eta=1e-5, fixEnrollTime = TRUE)$simd
#' data_final<-data.trim.d(100,data_temp)[[1]]
#' rho=1
#' gamma=0
#' # compare the 3 different ways below:
#' #library(survival)
#' sqrt(survdiff(Surv(survival,delta)~trt, data =data_final,rho=rho)$chisq)
#' FH.test(survival=data_final$survival,delta=data_final$delta,trt=data_final$trt,rho=rho,gamma=gamma)
#' WLR.test(survival=data_final$survival,delta=data_final$delta,trt=data_final$trt,w=function(...){survKM_minus(...)^rho*(1-survKM_minus(...))^gamma})
#' }
#' 
#'
WLR.test <- function(survival, delta, trt, w = function(v,...){1}){
  WLR_table <- logrank.table(survival, delta, trt)
  wt_vec <- w(
    v = WLR_table$survival, 
    survival = WLR_table$survival, 
    delta = WLR_table$delta
    )
  O1 <- sum(wt_vec * WLR_table$O1)
  E1 <- sum(wt_vec * WLR_table$E1)
  V <- sum(wt_vec^2 * WLR_table$Cov)
  Z <- (O1 - E1) / sqrt(V)
  return(list(O1 = O1, E1 = E1, Z = Z))
}
