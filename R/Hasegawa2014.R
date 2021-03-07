#' Sample size calculation for Fleming-Harrington weighted log-rank tests
#' This sample size calculation method was proposed by Hasegawa (2014).
#' This function is to calculate the sample size for Fleming-Harrington weighted log-rank tests with piece-wise exponential distributed survival curves in described in Hasegawa(2014).
#' 
#' @param eps The change point, before which, the hazard ratio is 1, and after which, the hazard ratio is theta
#' @param p Treatment assignment probability.
#' @param b The number of subintervals per time unit.
#' @param tau The end of the follow-up time  in the study. Note that this is identical to \eqn{T+\tau} in the paper from Hasegawa (2014).
#' @param omega The minimum follow-up time for all the patients.  Note that Hasegawa(2014) assumes that the accrual is uniform between time 0 and R, and there does not exist any censoring except for the administrative censoring at the ending time \eqn{\tau}. Thus this value omega is equivalent to \code{tau-R}. Through our simulation tests, we found that this function is quite robust to violations of these assumptions: dropouts, different cenosring rates for two  arms, and changing accrual rates. 
#' @param lambda The hazard for the control group.
#' @param lambda.trt The hazard for the treatment group after time eps.
#' @param rho The first parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.
#' @param gamma The second parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.
#' @param alpha Type I error.
#' @param beta Type II error. 
#' @return 
#' \item{n}{The needed sample size.}
#' \item{n_event}{The needed  event numbers for both arms together.}
#' \item{E.star}{The unit mean, correspoinding to \eqn{E^*} in Hasegawa(2014)}
#' \item{sum_D}{The cumulative D, and ceiling(n*D) is quivalent to n_vent. }
#' @references 
#' Ye, T., & Yu, M. (2018). A robust approach to sample size calculation in cancer immunotherapy trials with delayed treatment effect. Biometrics, 74(4), 1292-1300.
#' Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.
#' @author 
#' Lili Wang, Ting Ye
#' @note 
#' This function is based on a R function from Dr. Ting Ye's paper
#' : Ye, T., & Yu, M. (2018). A robust approach to sample size calculation in cancer immunotherapy trials with delayed treatment effect. Biometrics, 74(4), 1292-1300.
#' @examples
#' \dontrun{ 
#' # Example 1 from Hasegawa (2014)
#' p<-2/3
#' tau<-66
#' omega<-18
#' eps<-6
#' m1=21.7  #median survival time for placebo group
#' m2=25.8  # median survival time for treatment group
#' lambda<-log(2)/m1
#' lambda.trt<-log(2)*(m1-eps)/(m2-eps)/m1
#' theta=lambda.trt/lambda
#' alpha<-0.025
#' beta<-0.1
#' rho=0
#' gamma=1
#' b=30
#' sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)$n
#' #1974, identical to the paper's report
#' }
sample.size_FH <- function(
  eps, 
  p, 
  b, 
  tau, 
  omega,
  lambda,
  lambda.trt,
  rho, 
  gamma,
  alpha,
  beta
  ){
  # every sequence starts from 0
  # b<-30
  z.alpha <- qnorm(p = alpha, lower.tail = FALSE)
  z.beta <- qnorm(p = beta, lower.tail = FALSE)
  n_sub <- floor(b * tau)
  t <- c(0, seq(1, n_sub) / b)
  h_1 <- rep(lambda, (n_sub + 1)) #control
  h_2 <- c(
    rep(lambda, round(eps * b)),
    rep(lambda.trt, n_sub - round(eps * b) + 1)
    ) #treatment
  N_1 <- rep((1 - p),(n_sub + 1))
  N_2 <- rep(p, (n_sub + 1))
  for(i in 1:(n_sub-1)){
    N_1[i + 1] <- N_1[i] * (1 - h_1[i] / b - (t[i] > omega) / 
                            b / (tau - t[i]))
    N_2[i + 1] <- N_2[i] * (1 - h_2[i] / b-(t[i] > omega) / 
                              b / (tau - t[i]))
  }
  N_1[n_sub + 1] <- N_2[n_sub + 1] <- 0

  f_S_1 <- function(x) exp( - lambda * x)
  f_S_2 <- function(x) (x < eps) * exp( - lambda*x) + 
    (x >= eps) * exp( - (lambda * eps + lambda.trt * (x - eps)))
  #f_S_2_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-eps*lambda.trt*(1/theta-1))*exp(-lambda.trt*x)
  S_1 <- f_S_1(t)
  S_2 <- f_S_2(t)
  S <- (1 - p) * S_1 + p * S_2
  D <- (h_1 * N_1 + h_2 * N_2) / b * min(tau / R, 1) 
  theta_seq <- h_2 / h_1
  phi <- N_2 / N_1
  r <- S^rho * (1 - S)^gamma
  num_vec <- D * r * (phi * theta_seq / (1 + phi * theta_seq) - 
                        phi / (1 + phi))
  den_vec <- D * r^2 * phi / (1 + phi)^2
  E.star_num <- sum(num_vec[1:n_sub])
  E.star_den <- sqrt(sum(den_vec[1:n_sub]))
  E.star <- E.star_num / E.star_den
  n <- (z.alpha + z.beta)^2 / E.star^2
  n_event <- sum(D) * n
  return(list(
    n = ceiling(ceiling(n * p) / p), 
    n_event = ceiling(ceiling(n_event*p)/p),
    E.star = E.star,
    sum_D = sum(D[1:n_sub]),
    D = D[1:n_sub],
    den_vec = den_vec[1:n_sub],
    num_vec = num_vec[1:n_sub],
    time_vec = seq(1, n_sub) / b))
}

#' Calcualte the average hazard ratios
#' 
#' Calculate the average hazard ratios according to Kalbfleisch and Prentice (1981) or in the paper Hasegawa (2014) for piece-wise exponential survival functions (only one change point \code{eps}).
#' 
#' @param theta hazard ratio after eps between the treatment and the control group, assuming that the hazard rato is 1 before eps.
#' @param eps The change point, before which, the hazard ratio is 1, and after which, the hazard ratio is theta.
#' @param lambda The constant hazard for the control arm.
#' @param p Treatment assignment probability.
#' @references 
#' Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.
#' @author 
#' Lili Wang
#' @examples 
#' \dontrun{
#' # test
#' lambda=log(2)/6
#' theta=0.7
#' eps=2
#' avg.haz(theta,eps,lambda)
#' }
avg.haz <- function(theta, eps, lambda, p = 1 / 2){
  term1 <- 1 / 2 / p * (1 - exp(-2 * p * lambda * eps)) + 
    theta / p /(1 + theta) * exp( 
      - lambda * eps * (1 + p - theta * (1 - p))
      )
  term2 <- 1 / 2 / p * (1 - exp( - 2 * p * lambda * eps)) + 
    1 / p / (1 + theta) * exp(
      -lambda * eps * (1 + p - theta * (1 - p))
      )
  term1 / term2
}

