#library(data.table) #imports >= 1.12.2
#' Basic functions
#'
#' Some basic functions for information prediction.
#'
#' To prepare the values for the prediction of information values. The control arm is following an exponential with rate \code{lambda}, the treatment arm is piece-wise exponential with hazard ratio with respect to the control arm to be 1 before the changing point \code{eps}, and \code{theta} after the change point. Details can be found in the appendix of the reference paper. 
#'
#' @param  eps Change point.
#' @param theta Hazard ratio after the change point (before the change point HR should be 1).
#' @param lambda Event hazard for the control arm.
#' @param t.star Time point we pause the study to check the cumulative results.
#' @param R End of the accrual period.
#' @param k,k1,k2,m Parameters to control the exponential power of the survival functions (the control arm for the null hypothesis or the weighted sum of two arms for the alternative hypothesis).
#' @param p Treatment assignment probability.
#' @param e Some convenience parameter to control the change point, which is usually set to be eps or tau
#' @author Lili Wang
#'
#' @return \code{getc} returns the \eqn{\exp(-\lambda*\epsilon*(1-\theta))} which is a multiplier for the survival and hazard of the treatment arm after the change point \code{eps}.
#'
#'@references
#' Wang, L., Luo, X., & Zheng, C. (2021). A Simulation-free Group Sequential Design with Max-combo Tests in the Presence of Non-proportional Hazards. Journal of Pharmaceutical Statistics.
#' 
## Basic functions which are derived in the derivation document for piece-wice exponential distributed survival curves.
getc <- function(theta, lambda, eps){
  exp( - (1 - theta) * lambda * eps)
  }
# Note that if e=tau, it is identical to 1/k+1/(k^2)/lamb/R*exp(-k*lamb*tau)*(1-exp(k*lamb*R))
#' @rdname getc
uv <- function(e, k, lambda, R, t.star){
  ifelse(
    e > (t.star - R),
    u(e, k, lambda, R, t.star),
    v(e, k, lambda)
    )
}
#' @rdname getc
v <- function(e, k, lambda){
  1 / k * (1 - exp( - k * lambda * e))
}
#' @rdname getc
u <- function(e, k, lambda, R, t.star){
  1 / k * (1 - exp( - k * lambda * max(t.star - R, 0))) + 
    1 / k * min(1, t.star / R) * exp( - k * lambda * max(t.star - R, 0))-
    1 / k * max(t.star - e, 0) / R * exp( - k * lambda * min(e, t.star))+
    1 / k^2 / R / lambda * (exp( - k * lambda * min(e, t.star)) - 
                              exp( - k * lambda * max(t.star - R, 0))
                            )
}
#' @rdname getc
h1 <- function(k1, k2, lambda, theta, eps, R, t.star){
  uv(eps, k1 + k2 + 1, lambda, R, t.star) + 
    getc(theta, lambda, eps)^(k1 + 1) * theta *(
      u(t.star, theta * (k1 + 1) + k2,lambda,R,t.star) - 
        uv(eps,theta*(k1+1)+k2,lambda,R,t.star)
      )
}
#' @rdname getc
h0 <- function(k1, k2, lambda, theta, eps, R, t.star){
  uv(eps,k1 + k2 + 1, lambda, R, t.star) + 
    getc(theta, lambda, eps)^k1 * (
      u(t.star,theta * k1 + k2 + 1, lambda, R, t.star) - 
        uv(eps, k1 * theta + k2 + 1, lambda, R, t.star)
      )
}

# under H1
# h.tilde is introduced for the convenient calculation of the information under H1
#' @rdname getc
h.tilde <- function(m, lambda, theta, eps, R, p, t.star){
  out <- 0
  for(i in 0:m){
    out = out + factorial(m) / factorial(i) / factorial(m-  i) * 
      (
        p^(i + 1) * (1 - p)^(m - i) * 
          h1(i, m - i, lambda, theta, eps, R, t.star) + 
        p^i * (1 - p)^(m - i + 1) * h0(i, m - i, lambda, theta, eps, R, t.star)
       )
    #print(factorial(m)/factorial(i)/factorial(m-i))
  }
  return(out)
}
#' Precdict information/covariance under null hypothesis
#'
#' Calulcation of the information/covariance based on a presumed survival function under the null.
#'
#' This function is prepared to calculate the predicted information/covariance purely based on the assumed survival function under the null hypothesis: an exponential distribution with hazard \code{lambda}.
#'
#' @param rho First power parameter for the Fleming-Harrington weight which weighs on the early departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#' @param gamma Second power parameter for the Fleming-Harrington weight which weighs on the late departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#' @param rho1,rho2 First power parameters for the two Fleming-Harrington weights, defined for covariance calculation.
#' @param gamma1,gamma2 Second power parameters for the two Fleming-Harrington weights, defined for covariance calculation.
#' @param lambda Event hazard for the control arm.
#' @param R End of the accrual period.
#' @param p Treatment assignment probability.
#' @param t.star Time point we pause the study to check the cumulative information under the null.
#' @author Lili Wang.
#' @seealso \code{\link{I.1}}
#' @references
#' Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.
#'
I.0 <- function(rho, gamma, lambda, R, p, t.star){
  if(rho == 0 & gamma == 0){
    u(t.star, 1, lambda, R, t.star) * p * (1 - p)
  }else if (rho == 1 & gamma == 0){
    u(t.star, 3, lambda, R, t.star) * p * (1 - p)
  }else if(rho == 0 & gamma == 1){
    (u(t.star, 1, lambda, R, t.star) + 
       u(t.star, 3, lambda, R, t.star) - 
       2 * u(t.star, 2, lambda, R, t.star)) * p * (1 - p)
  }else if(rho == 1 & gamma == 1){
    (u(t.star, 3, lambda, R, t.star) + 
       u(t.star, 5, lambda, R, t.star) - 
       2 * u(t.star, 4, lambda, R, t.star)) * p * (1 - p)
  }else{
    stop("improper input of rho or gamma, they must be either 0 or 1. ")
  }
}

#' @rdname I.0
I.0.cov <- function(rho1, gamma1, rho2, gamma2, lambda, R, p, t.star){
  if(rho1 + rho2 == 0 & gamma1 + gamma2 == 0){
    u(t.star, 1, lambda, R, t.star) * p * (1 - p)
  }else if (rho1 + rho2 == 1 & gamma1 + gamma2 == 0){
    u(t.star, 2, lambda, R, t.star) * p * (1 - p)
  }else if (rho1 + rho2 == 2 & gamma1 + gamma2 == 0){
    u(t.star, 3, lambda, R, t.star) * p * (1 - p)
  }else if (rho1 + rho2 == 0 & gamma1 + gamma2 == 1){
    (u(t.star, 1, lambda, R, t.star) - 
       u(t.star, 2, lambda, R, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 0 & gamma1 + gamma2 == 2){
    (u(t.star, 1, lambda,R, t.star) + 
       u(t.star, 3, lambda, R, t.star) - 
       2 * u(t.star, 2, lambda, R, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 1 & gamma1 + gamma2 == 1){
    (u(t.star, 2, lambda, R, t.star) - 
       u(t.star, 3, lambda, R, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 2 & gamma1 + gamma2 == 1){
    (u(t.star, 3, lambda, R, t.star) - 
       u(t.star, 4, lambda, R, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 1 & gamma1 + gamma2 == 2){
    (u(t.star, 2, lambda, R, t.star) + 
       u(t.star, 4, lambda, R, t.star) - 
       2 * u(t.star, 3, lambda, R, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 2 & gamma1 + gamma2 == 2){
    (u(t.star, 3, lambda, R, t.star) + 
       u(t.star, 5, lambda, R, t.star) - 
       2 * u(t.star, 4, lambda, R, t.star)) * p * (1 - p)
  }else{
    stop("improper input of rho or gamma, they must be either 0 or 1. ")
  }
}



#' Predicted information/covariance under the alternative hypothesis
#'
#' Calulcation of the information/covariance based on a presumed survival function under the alternative hypothesis.
#'
#' This function is prepared to calculate the predicted information/covariance purely based on the assumed survival function under the alternaitve hypothesis: the control group is following an exponential distribution with hazard \code{lambda}, while the treatment group is following a piece-wise exponential distribution with same hazard before \code{eps}, but a hazard equals \code{theta} times the \code{lambda} after \code{eps}.
#'
#' @inheritParams I.0
#' @param eps Change point.
#' @param theta Hazard ratio after the change point (before the change point HR should be 1).
#' @author Lili Wang.
#' @seealso \code{\link{I.0}}
#' @references
#' Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.
#'
I.1 <- function(rho, gamma, lambda, theta, eps, R, p, t.star){

  if(rho == 0 & gamma == 0){
    h.tilde(0, lambda, theta, eps, R, p, t.star) * p * (1 - p)
  }else if(rho == 1 & gamma == 0){
    h.tilde(2, lambda, theta, eps, R, p, t.star) * p * (1 - p)
  }else if(rho == 0 & gamma == 1){
    (h.tilde(0, lambda, theta, eps, R, p, t.star) + 
       h.tilde(2, lambda, theta, eps, R, p, t.star) - 
       2 * h.tilde(1, lambda, theta, eps, R, p, t.star)) * p * (1 - p)
  }else if(rho == 1 & gamma == 1){
    (h.tilde(2, lambda, theta, eps, R, p, t.star) + 
       h.tilde(4, lambda, theta, eps, R, p, t.star) - 
       2 * h.tilde(3, lambda, theta, eps, R, p, t.star)) * p * (1 - p)
  }else{
    stop("improper input of rho or gamma, they must be either 0 or 1. ")
  }
}

#' @rdname I.1
I.1.cov <- function(
  rho1, 
  gamma1, 
  rho2, 
  gamma2, 
  lambda, 
  theta, 
  eps, 
  R, 
  p, 
  t.star
  ){
  if(rho1 + rho2 == 0 & gamma1 + gamma2 == 0){
    h.tilde(0, lambda, theta, eps, R, p, t.star) * p * (1 - p)
  }else if (rho1 + rho2 == 1 & gamma1 + gamma2 == 0){
    h.tilde(1, lambda, theta, eps, R, p, t.star) * p * (1 - p)
  }else if (rho1 + rho2 == 2 & gamma1 + gamma2 == 0){
    h.tilde(2, lambda, theta, eps, R, p, t.star) * p * (1 - p)
  }else if (rho1 + rho2 == 0 & gamma1 + gamma2 == 1){
    (h.tilde(0, lambda, theta, eps, R, p, t.star) - 
       h.tilde(1, lambda, theta, eps, R, p, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 0 & gamma1 + gamma2 == 2){
    (h.tilde(0, lambda, theta, eps, R, p, t.star) + 
       h.tilde(2, lambda, theta, eps, R, p, t.star) - 
       2 * h.tilde(1, lambda, theta, eps, R, p, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 1 & gamma1 + gamma2 == 1){
    (h.tilde(1, lambda, theta, eps, R, p, t.star) - 
       h.tilde(2, lambda, theta, eps, R, p, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 2 & gamma1 + gamma2 == 1){
    (h.tilde(2, lambda, theta, eps, R, p, t.star) - 
       h.tilde(3, lambda, theta, eps, R, p, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 1 & gamma1 + gamma2 == 2){
    (h.tilde(1, lambda, theta, eps, R, p, t.star) + 
       h.tilde(3, lambda, theta, eps, R, p, t.star) - 
       2 * h.tilde(2, lambda, theta, eps, R, p, t.star)) * p * (1 - p)
  }else if (rho1 + rho2 == 2 & gamma1 + gamma2 == 2){
    (h.tilde(2, lambda, theta, eps, R, p, t.star) + 
       h.tilde(4, lambda, theta, eps, R, p, t.star) - 
       2*h.tilde(3 , lambda , theta , eps , R , p , t.star)) * p * (1 - p)
  }else{
    stop("improper input of rho or gamma, they must be either 0 or 1. ")
  }
}


#' Predicted cross-test correlation
#'
#' These two functions are to predict the correlation between two weighted log-rank tests at certain time \code{t.star} under either the null hypothesis (using \code{cor.0}) or the alternative hypothesis (using \code{cor.1}).
#'
#' These two functions are designed to calculate the predicted correlation between the two weighted log-rank tests at time \code{t.star} under the two hypotheses. The null hypothesis is an exponential distribution for both the treatment and control arms with hazard \code{lambda}, while the alternative hypothesis has the control group following an exponential distribution with hazard \code{lambda}, and the treatment group following a piece-wise exponential distribution with hazard \code{lambda} before \code{eps}, but a hazard \code{theta} times \code{lambda} after \code{eps}.
#'
#' @inheritParams I.1
#' @author Lili Wang.
#' @references
#' Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.
#'
cor.0 <- function(rho1, gamma1, rho2, gamma2, lambda, R, p,t.star){
  I.0.cov(
    rho1 = rho1, 
    gamma1 = gamma1, 
    rho2 = rho2,
    gamma2 = gamma2,
    lambda = lambda,
    R = R, 
    p = p,
    t.star = t.star)/
    sqrt(
      I.0(rho = rho1, 
          gamma = gamma1, 
          lambda = lambda, 
          R = R, p = p, 
          t.star = t.star)
      )/
    sqrt(
      I.0(rho = rho2,
          gamma = gamma2,
          lambda = lambda,
          R = R,
          p = p,
          t.star = t.star)
      )
}
#' @rdname cor.0
cor.1 <- function(
  rho1,
  gamma1,
  rho2,
  gamma2,
  lambda,
  theta,
  eps,
  R,
  p,
  t.star
  ){
  I.1.cov(
    rho1 = rho1,
    gamma1 = gamma1,
    rho2 = rho2,
    gamma2 = gamma2,
    lambda = lambda,
    theta = theta,
    eps = eps,
    R = R,
    p = p,
    t.star = t.star)/
    sqrt(
      I.1(rho = rho1,
          gamma = gamma1,
          lambda = lambda,
          theta = theta,
          eps = eps,
          R = R,
          p = p,
          t.star = t.star))/
    sqrt(
      I.1(rho = rho2,
          gamma = gamma2,
          lambda = lambda,
          theta = theta,
          eps = eps,
          R = R,
          p = p,
          t.star = t.star)
      )
}


#' Estimated information based on the data
#'
#' Estimate the information based on the data, which is the numerator of the information fraction.
#'
#' The \code{I_t} function estimates the information up to the maximum follow-up time in the data of \code{data_check}, which is identical to the numerator of the information fraction proposed by Hasegawa (2016):\eqn{ \hat{P}_1(t)\hat{P}_0(t)\int_0^t W(t,s)^2N(t,ds)}. Note that the datasets \code{data_check} and \code{data_ref} input here are output data from \code{data.trim} functions, or any datasets including \code{survival} as time to event or censoring, \code{delta} as event indicators, and \code{trt} denotes treatment assignment (1 is treatment, 0 is control). Note that \code{I_t.2} is another option which is slightly different from the one proposed in Hasegawa(2016), but is identical to the estimate of variance of the weighted log-rank test, which considers the total at-risk set \eqn{R(t)} and treatment arm \eqn{R_1(t)}: \eqn{\int_0^t\frac{R_1(s)R_0(s)}{R(s)^2}W(t,s)^2N(t,ds)}.
#'
#' @param rho First power parameter for the Fleming-Harrington weight which weighs on the early departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#'
#' @param gamma Second power parameter for the Fleming-Harrington weight which weighs on the late departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#'
#' @param data_ref Input reference dataset which provides the survival curves for the estimation. It could be some dataset entirely external to \code{data_check}. This dataset should include at lease the 3 variables: \code{survival} for the time to event or censoring, \code{delta} as the event indicator, and \code{trt} for the treatment assignment indicator. It will perfectly fit the output dataset from the \code{data.trim} functions.
#'
#' @param data_check Input dataset to check the estimated information. It should follow the sample format as \code{data_ref}, which includes three variables: \code{survival}, \code{delta} and \code{trt}.
#'
#' @return The returned value is the calculated information estimated from the input dataset \code{data_check} using the survival function estimated from \code{data_ref}.
#'
#' @examples
#' \dontrun{
#' # install.packages("devtools")
#' # library(devtools)
#' # install_github("keaven/nphsim")
#' library(nphsim)
#' eps<-2 # delayed effect
#' p<-0.5 #treatment assignment
#' b<-30 # an intrinsic parameter to decide the number of intervals per time unit
#' tau<- 18 # end of the study
#' R<-14 # accrual period [0,R]
#' omega<- tau-R
#' lambda<-log(2)/6 # control group risk hazard
#' theta<-0.7 # hazard ratio
#' lambda.trt<- lambda*theta #hazard after the change point for the treatment arm
#' rho<- 0 # parameter for the weights
#' gamma<-1 #parameter for the weights
#' alpha<-0.025 #type 1 error
#' beta<-0.1 #type 2 error
#' # First we decide the sample size:
#' size_FH <- sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)
#' n_FH <-size_FH$n
#' n_event_FH<-size_FH$n_event
#' accrual.rt<-n_FH/R # the needed arrual rate
#' #Generate data accordingly, use eta=1e-5 to inhibit censoring
#' data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p), gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd
#'
#' #Obtain the full information at the final stage based on the generated data
#' #Trim the data up to the final stage when n_event_FH events have been observed
#' data_temp1 <-data.trim.d(n_event_FH,data_temp)[[1]]
#' I_t(data_temp1,data_temp1,rho,gamma) # the estimated information at the final stage
#' #Trim the data up to certain event numbers at the interim stage when 60% of the events have been observed. Have been trimmed once to get data_temp1, no need to add additional variables, thus set the third argument to be F.
#' I_t.2(data_temp1,data_temp1,rho,gamma) # If we consider the change of the at-risk set, which is not necessary to be a fixed probability.
#' data_temp2 <- data.trim.d(ceiling(0.6*n_event_FH),data_temp1,F)[[1]]
#' I_t(data_temp1,data_temp2,rho,gamma) # Use the full dataset data_temp to provide the survival function, and check the estimated information for the trimmed data set data_temp2 with only 60% of the planned events have been observed.
#' I_t.2(data_temp1,data_temp2,rho,gamma) # If we consider the change of the at-risk set, which is not necessary to be a fixed probability.
#' }
#'
#' @seealso \code{\link{data.trim}}
#' @author Lili Wang
#' @references
#' Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.
#'
I_t <- function(data_ref, data_check, rho, gamma){
  p_hat <- mean(data_ref$trt) #obtain the estimated treatment assignment probability
  sum_KM_t <- data.table(
    surv = survKM_minus(
      v = data_check$survival, 
      survival = data_ref$survival,
      delta = data_ref$delta
      ),
    n.event = data_check$delta
   ) #obtain the survival functions for each subjects in data_check
  {if(rho == 0 && gamma == 0) 
    p_hat * (1 - p_hat) * sum(sum_KM_t$n.event) # For the G^(0,0) class
    else if(rho == 0 && gamma == 1) 
      p_hat * (1 - p_hat) * sum(sum_KM_t$n.event * (1 - sum_KM_t$surv)^2)
    else if(rho == 1 && gamma == 0 ) 
      p_hat * (1 - p_hat) * sum(sum_KM_t$n.event * (sum_KM_t$surv)^2)
    else if(rho == 1 && gamma == 1 ) 
      p_hat * (1 - p_hat) * sum(sum_KM_t$n.event * 
                                  ((sum_KM_t$surv) * (1 - sum_KM_t$surv))^2)
    else stop("Invalid input of rho or gamma: (0,0), (0,1), (1,0), or (1,1)")
  }
}
#' @rdname I_t
I_t.2 <- function(data_ref, data_check, rho, gamma){
  LRT.table <- logrank.table(
    survival = data_check$survival,
    delta = data_check$delta,
    trt = data_check$trt
    )
  surv <- survKM_minus(
    v = LRT.table$survival,
    survival = data_ref$survival,
    delta = data_ref$delta
    )
  {if(rho == 0 && gamma == 0) sum(LRT.table$Cov) # For the G^(0,0) class
    else if(rho == 0 && gamma == 1)  sum(LRT.table$Cov * (1 - surv)^2)
    else if(rho == 1 && gamma == 0 ) sum(LRT.table$Cov * (surv)^2)
    else if(rho == 1 && gamma == 1 ) sum(LRT.table$Cov * (surv * (1 - surv))^2)
    else stop("Invalid input of rho or gamma: (0,0), (0,1), (1,0), or (1,1)")
    }
}

#' Information fraction for Fleming-Harrington weighted log-rank test
#'
#' Monitor the raction for Fleming-Harrington weighted log-rank test for a vector of time points
#'
#' Calculation the information fraction for Fleming-Harrington family weighted log-rank tests using the monitored estimated information for numerator, and the predicted information \eqn{I_{max}} as denominator.
#'
#' @param data There are two possible structures allowed for this input data. The first type needs to have \code{trimmed=F} and include variables: a \code{treatment} variable with "experimental" denoting treatment group, \code{cnsr} variable with value 1 denoting censoring, \code{ct} variable denoting event time from the origin of the study, which equals the sum of entering time \code{enterT} and the survival time (time to event or censoring). A dataset simulated from from R package \href{https://github.com/keaven/nphsim}{nphsim} should fit the first type well enough (see the example1). The second type can be any data.frame or data.table output from a \code{data.trim} function, including variables: \code{ct} denoting event time from the origin of the study or the sum of entering time and the survival time, \code{survival} denoting the survival time or time to event/censoring, \code{delta} as an event indicator, \code{enterT} as entering time (example 2). For the second type, we set \code{trimmed=T} to avoid extra computations, but should be fine if \code{trimmed=F}.
#'
#' @param t_vec Follow-up time since the origin of the study (not that it's not following the survival time scale, but following the calendar time scale ), which could be a vector, to measure the information fraction for these time points.
#' @param I_max The evaluated \eqn{I_{max}}, which returned from function \code{I.1} or \code{I.0} by setting the \code{t.start=tau}, where \code{tau} is the end of the study.
#' @param rho First power parameter for the Fleming-Harrington weight which weighs on the early departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#' @param gamma Second power parameter for the Fleming-Harrington weight which weighs on the late departures: \eqn{S(t^-)^\rho(1-S(t^-))^\gamma}.
#' @param trimmed Logical indicator to show whether the \code{data} input has been "trimmed" by \code{data.trim} and \code{data.trim.d} before: adding variables like \code{delta} indicating events (=1), and \code{trt} distringuishing the treatment group (=1) from the control group (=0)
#' @return This function returns a vector of information fractions corresponding to the input time vector \code{t_vec}.
#' @seealso \code{\link{data.trim}}
#' @author Lili Wang
#' @references
#' Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.
#' @examples
#' \dontrun{
#' # install.packages("devtools")
#' # library(devtools)
#' # install_github("keaven/nphsim")
#' # library(nphsim)
#' eps<-2 # delayed effect
#' p<-0.5 #treatment assignment
#' b<-30 # an intrinsic parameter to decide the number of intervals per time unit
#' tau<- 18 # end of the study
#' R<-14 # accrual period [0,R]
#' omega<- tau-R
#' lambda<-log(2)/6 # control group risk hazard
#' theta<-0.7 # hazard ratio
#' lambda.trt<- lambda*theta #hazard after the change point for the treatment arm
#' rho<- 0 # parameter for the weights
#' gamma<-1 #parameter for the weights
#' alpha<-0.025 #type 1 error
#' beta<-0.1 #type 2 error
#' # First we decide the sample size:
#' size_FH <- sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)
#' n_FH <-size_FH$n
#' n_event_FH<-size_FH$n_event
#' accrual.rt<-n_FH/R # the needed arrual rate
#' #Generate data accordingly, use eta=1e-5 to inhibit censoring
#' data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p), gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd
#'
#' # Example 1 for FH.frac.cal: Set trimmed=F and work on the crude dataset from nphsim()
#' inf_frac_vec1<-FH.frac.cal(data_temp,c(10,15,18),I_denom,rho,gamma,trimmed=F)
#' inf_frac_vec1
#'
#' # Example 2 for FH.frac.cal: First trim the data before inputting into FH.frac.cal() setting trimmed=T, and obtain the whole spectrum.
#' I_denom<-I.1(rho, gamma,lambda,theta,eps,R,p,t.star=tau)*n_FH
#' tau.star=21 #in case the ratio=1 when t>tau
#' #Trim the data
#' data_temp2 <-data.trim(tau.star,data_temp)
#' t_seq <- seq(0.1,tau.star,0.1) # the time series to check the information fraction
#' inf_frac_vec2<-FH.frac.cal(data_temp2,t_seq,I_denom,rho,gamma,trimmed=T)
#' # WLRT at the interim
#' interim_index<- which.min(abs(inf_frac_vec2-0.6))
#' interim_time<-t_seq[interim_index]
#' interim_frac<-inf_frac_vec2[interim_index]
#'  # WLRT at the final
#'  final_index<- which.min(abs(inf_frac_vec2-1))
#'  final_time<-t_seq[final_index]
#'  final_frac<-inf_frac_vec2[final_index]
#' }

FH.frac.cal <- function(data, t_vec, I_max, rho, gamma, trimmed){
  t_vec <- as.vector(t_vec) # force t_vec to become a vector
  sapply(as.vector(unlist(t_vec)), function(tt){
    dtemp <- data.trim(tt, data, trimmed)
    #p_hat<-mean(dtemp$trt)
    # The problem of the commented codes below is that it is using exact survival S(t) not S(t^-)
    #KM_t<-tryCatch(survfit(Surv(survival,delta)~0,data=dtemp),
     #              error=function(e) NA)
   # sum_KM_t<-summary(KM_t)[c("time","surv","n.event")]
    #I_numer_t<-sum(p_hat*(1-p_hat)*(sum_KM_t$surv^rho*(1-sum_KM_t$surv)^gamma)^2*sum_KM_t$n.event)
    #I_numer_t/I_max
    I_t(dtemp, dtemp, rho,gamma) / I_max # unlike using survfit and summary(survfit object), which should use the survival function one row above because of the Fleming-Harrington weights are predictable.
  })
}

#' Estimate the covariance and correlation between two arbitrary weights
#'
#' These two functions estimate the covariance and correlations between the two arbitrary weight functions, which are not necessary to be Fleming-Harrington family.
#'
#' Any two weight functions can be assigned to arguments \code{w1} and \code{w2}. Two examples, one is Fleming-Harrington family and the other is not, are demonstrated in the examples section.

#' @param survival The time to event or censoring, not that, it's the follow-up time after entoring, you may also consider as the total at-risk time.
#' @param delta The event indicator, with 1 indicating observed events, and 0 indicating censoring.
#' @param trt The treatment assignment indicator, with 1 indicating treatment group, and 0 as control group.
#' @param w1 It has the default function which will return standard log-rank test with weight 1 and thus the function will be reduced to a variance for log-rank tests, and correlation always equals 1. If the two weights are identical, \code{WLR.test.cov} is equivalent to the estimated variance, and \code{WLR.test.cor} is always equal to 1. The function can be any non-negative functions with a basic argument v as the input time vector, which are corresponding to the follow-up times. Optionally, there are two additional variables,follow-up time \code{survival} and event indicator \code{delta} to make the weights dependent on the survival functions (like the Fleming-Harrington family). It would be better if the function itself has \code{...} as the last argument, so that it can be robust to any misspecification of the variable names, and thus, it will just ignore the misspecified ones. Please refer to the examples to figure out how to define the Fleming-Harrington and any other weight functions.
#' @param w2 Same requirements as the other argument \code{w1}. Just not that if they are identical, \code{WLR.test.cov} returns the variance like \code{I.t}, and  \code{WLR.test.cor} always returns 1.
#' @return The two functions, \code{WLR.test.cov} returns the covariance, \code{WLR.test.cor} returns the correlation coefficient estimate solely based on the input data.
#' @seealso \code{\link{cor.0}},\code{\link{cor.1}}, \code{\link{I_t}}.
#' @author Lili Wang
#' @references
#' Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.
#' @examples
#' \dontrun{
#' # install.packages("devtools")
#' # library(devtools)
#' # install_github("keaven/nphsim")
#' # library(nphsim)
#'
#' eps<-2 # delayed effect
#' p<-0.5 #treatment assignment
#' b<-30 # an intrinsic parameter to decide the number of intervals per time unit
#' tau<- 18 # end of the study
#' R<-14 # accrual period [0,R]
#' omega<- tau-R
#' lambda<-log(2)/6 # control group risk hazard
#' theta<-0.7 # hazard ratio
#' lambda.trt<- lambda*theta #hazard after the change point for the treatment arm
#' rho<- 0 # parameter for the weights
#' gamma<-1 #parameter for the weights
#' alpha<-0.025 #type 1 error
#' beta<-0.1 #type 2 error
#' # First we decide the sample size:
#' size_FH <- sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)
#' n_FH <-size_FH$n
#' n_event_FH<-size_FH$n_event
#' d_fixed<-ceiling(-0.6*n_event_FH)
#' accrual.rt<-n_FH/R # the needed arrual rate
#' #Generate data accordingly, use eta=1e-5 to inhibit censoring
#' data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt),
#'                   ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p),
#'                             gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd
#' # Example 1 for WLR.test.cov and WLR.test.cor: Fleming-Harrington family Weights
#' # I will let w1 be the default 1
#' # define a WLRT for w2 accodring to the rho and gamma defined above.
#' w2<-function(...){survKM_minus(...)^rho*(1-survKM_minus(...))^gamma}
#' data_interim<-data.trim.d(d_fixed,data_temp)[[1]] #data trimmed at the interim stage, the second enry on the list is the interim time, refer to function data.trim.d for details.
#' data_final<-data.trim.d(n_event_FH,data_temp)[[1]] #data trimmed at the final stage
#' WLR.test.cov(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2)
#' WLR.test.cor(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2)
#'
#' # The variance should be identical to the output from I_t and correlation is 1 if two weights are identical.
#' WLR.test.cov(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w1=w2,w2=w2)
#' I_t.2(data_interim,data_interim,rho,gamma)
#' WLR.test.cor(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w1=w2,w2=w2)
#'
#' #Example 2 for WLR.test.cov and WLR.test.cor: any Weights
#' w2_2<-function(v,...){1-exp(-v*0.25)}
#' WLR.test.cov(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2_2)
#' WLR.test.cor(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2_2)
#' }
#'
WLR.test.cov <- function(
  survival, 
  delta, 
  trt,
  w1 = function(v,...){1},
  w2 = function(v,...){1}){
  WLR_table <- logrank.table(survival, delta, trt)
  wt_vec1 <- w1(
    v = WLR_table$survival, 
    survival = WLR_table$survival,
    delta = WLR_table$delta
    )
  wt_vec2 <- w2(
    v = WLR_table$survival, 
    survival = WLR_table$survival,
    delta = WLR_table$delta
    )
  sum(wt_vec1 * wt_vec2 * WLR_table$Cov)
}
#' @rdname WLR.test.cov
WLR.test.cor <- function(
  survival, 
  delta, 
  trt,w1 = function(v,...){1}, 
  w2 = function(v,...){1}
  ){
  WLR_table <- logrank.table(survival, delta, trt)
  wt_vec1 <- w1(
    v = WLR_table$survival,
    survival = WLR_table$survival,
    delta = WLR_table$delta
    )
  wt_vec2 <- w2(
    v = WLR_table$survival,
    survival = WLR_table$survival,
    delta = WLR_table$delta
    )
  cov <- sum(wt_vec1 * wt_vec2 * WLR_table$Cov)
  v1 <- sum(wt_vec1^2 * WLR_table$Cov)
  v2 <- sum(wt_vec2^2 * WLR_table$Cov)
  (cor_val <- cov / sqrt(v1 * v2))
}


#' Approximate information for an arbitrary survival function
#'
#' An approximation alternative to the regular prediction of the information/covariance based on the assumed survival functions.
#'
#' @param t.star The ending time of the cumulative informaiton or covariance prediciton.
#' @param p Treatment assignment probability.
#' @param S1 Survival function for the treatment group.
#' @param S0 Survival function for the control gorup.
#' @param func The integrand function.
#' @param n.length The number of intervals spitted to obtain the approximate integration.
#' @author Lili Wang
#' @references
#' Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.
#' @examples
#' \dontrun{
#' # Examples for approx.I
#' eps<-2 # delayed effect
#' p<-0.5 #treatment assignment
#' tau<-18 # end of the study
#' R<-14 # accrual period [0,R]
#' lambda<-log(2)/6 # control group risk hazard
#' theta<-0.7 # hazard ratio
#' lambda.trt<- lambda*theta
#' rho<- 0 # parameter for the weights
#' gamma<-1 #parameter for the weights
#' S1<-function(x){
#'   ifelse(x>eps,exp(-theta*lambda*x)*getc(theta,lambda,eps),exp(-lambda*x))
#'   }
#'   S0<-function(x){
#'     exp(-lambda*x)
#'     }
#'     S_pool<-function(x){
#'       p*S1(x)+(1-p)*S0(x)
#'       }
#'       func<-function(x){
#'         min((tau-x)/R,1)*(S_pool(x)^rho*(1-S_pool(x))^gamma)^2
#'         }
#'  approx.I(t.star=tau,p,S1=S1,S0=S0,fun=func,n.length=1e6)
#'  I.1(rho,gamma,lambda,theta,eps,R,p,tau)
#'  # Change the cumulative information up to 10 instead of taus
#'  func2<-function(x){
#'    min((10-x)/R,1)*(S_pool(x)^rho*(1-S_pool(x))^gamma)^2
#'    }
#'    approx.I(t.star=10,p,S1=S1,S0=S0,fun=func2,n.length=1e6)
#'    I.1(rho,gamma,lambda,theta,eps,R,p,t.star=10)
#'    # Covariance approximation for two weights: 1 and G(0,1)
#'    rho1=rho2=0
#'    gamma1=0
#'    gamma2=1
#'    func3<-function(x){
#'   min((10-x)/R,1)*(S_pool(x)^rho1*(1-S_pool(x))^gamma1)*(S_pool(x)^rho2*(1-S_pool(x))^gamma2)
#'   }
#'   approx.I(t.star=10,p,S1=S1,S0=S0,fun=func3,n.length=1e6)
#'   I.1.cov(rho1,gamma1,rho2,gamma2,lambda,theta,eps,R,p,t.star=10)
#' }
#' 
approx.I <- function(
  t.star,
  p,
  S1 = function(x){1},
  S0 = function(x){1},
  func = function(x){1},
  n.length = 1e6
  ){
  out = 0
  x_seq <- seq(0, t.star, length.out = n.length)
  df <- function(x){
    p * S1(x) + (1 - p) * S0(x)
  }
  for(i in 2:n.length){
    x = x_seq[i]
    out = out - func(x) * (df(x_seq[i]) - df(x_seq[i - 1]))
  }
  out * p * (1 - p)
}




