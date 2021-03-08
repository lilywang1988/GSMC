#' Boundary calculation for GSMC
#' 
#' Boundary calculation for interim analysis with max-combo tests based on correlation matrix and the alpha spending function.
#'
#' Suppose there are 2 stages (1 interim, 1 final), and two tests for a max-combo in each stage, then we have totally 4 test statistics. Let the alpha spending function to be \code{c(alpha1,alpha)}, and the first two (\eqn{Z_{11},Z_{12}}) share one cutoff value z1, the latter two share another two (\eqn{Z_{21},Z_{22}}) share another cutoff value z2. Controlling the type I error is equivalent to ensuring that \eqn{P(Z_{11}<z_1,Z_{12}<z_1)=\alpha_1} and \eqn{P(Z_{11}<z_1,Z_{12}<z_1,Z_{21}<z_2,Z_{22}<z_2)=\alpha} are both satisfied. Note that the vector \eqn{[Z_{11},Z_{12},Z_{21},Z_{22}]^T\sim MVN(0,\Sigma_0)}. \code{Sigma0} corresponds to \eqn{\Sigma_0}, \code{index} records the ordered stages of each test statistics, whose order should be identical to the order of rows or columns in \code{Sigma0}.  Specifically, in this example, \code{index} should be \code{c(1,1,2,2)}. \code{alpha_sp} is the alpha spending function, which records how much type I error you would like to spend up to every stage.
#'
#' @param Sigma0 correlation matrix for all the test statistics.
#' @param index vector of non-decreasing integer starting from 1 to indicate which stage each column or row of the correlation matrix \code{Sigma0} corresponds to.
#' @param alpha_sp vector of cumulative errors to spend up to each stage.
#' @param n.rep number of repeats to take the median for output since the called likelihood generator of a multivariate normal distribution \code{\link[mvtnorm]{pmvnorm}} is not determinant. The default \code{n.rep} value is 5. 
#' @return
#' \item{z_alpha}{boundary values for all the stages.}
#' \item{z_alpha_vec}{boundary values for all the test statistics following the \code{index}. }
#' @author Lili Wang
#' @references 
#' Wang, L., Luo, X., & Zheng, C. (2021). A Simulation-free Group Sequential Design with Max-combo Tests in the Presence of Non-proportional Hazards. Journal of Pharmaceutical Statistics.
#' @examples
#'  \dontrun{
#'   #install.packages("gsDesign")
#'   library(gsDesign)
#'   alpha=0.025
#'   beta=0.1
#'   # If there are two stages (K=2), with on interim stage and a final stage
#'   # First we obtain the errors spent at each stage to be identical to the
#'    ones from regular interim analysis assuming that the interim stage
#'     happened at 60% of events have been observed. The error spending
#'      function used below is O\'Brien-Fleming.
#'   x <- gsDesign::gsDesign(
#'   k = 2, 
#'   test.type = 1, 
#'   timing = 0.6, 
#'   sfu = "OF", 
#'   alpha = alpha, 
#'   beta = beta,
#'   delta = -log(0.7))
#'   (z <- x$upper$bound)
#'   x
#'   Sigma0_v <- rep(0.5, 6)
#'   Sigma0 <- matrix(1, ncol = 4, nrow = 4)
#'   Sigma0[upper.tri(Sigma0)] <- Sigma0_v
#'   Sigma0[lower.tri(Sigma0)] <- t(Sigma0)[lower.tri(t(Sigma0))]
#'   Sigma0
#'   # The error you would like to spend at the interim stage:
#'   alpha_interim <- pnorm(z[1],lower.tail = F) 
#'   
#'   zz <- Maxcombo.bd(
#'   Sigma0 = Sigma0,
#'   index = c(1, 1, 2, 2),
#'   alpha_sp = c(alpha_interim, alpha))
#'   
#'  
#'   # boundary value for each stage
#'   zz$z_alpha 
#'   # boundary value for each test statistic correponding to index
#'   zz$z_alpha_vec 
#'   mvtnorm::pmvnorm(
#'   upper = rep(zz$z_alpha[1], 2),
#'   corr = Sigma0[1:2,1:2]
#'   )[[1]]
#'   
#'   1-alpha_interim
#'   1-mvtnorm::pmvnorm(
#'   upper = zz$z_alpha_vec,
#'   corr = Sigma0
#'   )[[1]]
#'   
#'   alpha
#'   # What if we do not consider interim stage but with only a final stage? 
#'   zz1 <- Maxcombo.bd(
#'   Sigma0 = Sigma0[3:4,3:4],
#'   index = c(1,1),
#'   alpha_sp = c(alpha)
#'   )
#'   mvtnorm::pmvnorm(
#'   upper = rep(zz1$z_alpha, 2),
#'   corr = Sigma0[1:2, 1:2]
#'   )[[1]]
#'   
#'   1-alpha
#'   # This function will also fit 2 or any number of interims (K>=3)
#'   # Let there are 3 stages, Let us try controlling the error spent 
#'   at each stage.
#'   stage_p <- c(0.5,0.7,0.8,0.9)
#'   x <- gsDesign::gsDesign(k=5, test.type=1, timing=stage_p, sfu="OF", 
#'   alpha=alpha, beta=beta,delta=-log(0.7))
#'   (z <- x$upper$bound)
#'   alpha_sp<- cumsum(x$upper$prob[,1]) # the theoretical cumulative
#'    errors spent at each stage
#' # 2 tests per stage
#' Sigma0_v<-rep(0.5,choose(10,2))
#' Sigma0<-matrix(1, ncol=10,nrow=10)
#' Sigma0[upper.tri(Sigma0)]<- Sigma0_v
#' Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(t(Sigma0))]
#' Sigma0
#'
#' zz<-Maxcombo.bd(Sigma0 = Sigma0,index=c(1,1,2,2,3,3,4,4,5,5),alpha_sp=alpha_sp)
#'
#' zz$z_alpha # boundary value for each stage
#' zz$z_alpha_vec # boundary value for each test statistic correponding to index
#' # interim 1
#' mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1],2),corr=Sigma0[1:2,1:2])[[1]] # expected error spent at this stage
#' 1-alpha_sp[1] #compare with the expected error spent at this stage
#' # above two rows are very close to each other, same for the following pairs.
#' # interim 2
#' mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:2],each=2),corr=Sigma0[1:4,1:4])[[1]]
#' 1-alpha_sp[2]
#' # interim 3
#' mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:3],each=2),corr=Sigma0[1:6,1:6])[[1]]
#' 1-alpha_sp[3]
#' # interim 4
#' mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:4],each=2),corr=Sigma0[1:8,1:8])[[1]]
#' 1-alpha_sp[4]
#' # final stage
#' mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:5],each=2),corr=Sigma0[1:10,1:10])[[1]]
#' 1-alpha_sp[5]
#'  }
Maxcombo.bd <- function(
  Sigma0, 
  index, 
  alpha_sp, 
  n.rep = 5
  ){
  K  <- length(alpha_sp)
  if(length(unique(index)) != K){
    stop("Error: The index should be consistent with the length of alpha_sp!")
    }
  if(any(diff(alpha_sp) <= 0)) 
    stop("Error: alpha_sp should be monotone (strictly) increasing!")
  if(any(diff(index) < 0)) 
    stop("Error: index should be monotone increasing!")
  z_alpha <- NULL
  z_alpha_vec <- NULL # a repeat vector to record the test statistic for each intex
  for(k in 1:K){
    use_index <- which(index <= k)
    Sigma0_use <- Sigma0[use_index, use_index]
    if(k == 1){
      z_alpha[1] <- median(
        replicate(
          n.rep,
          mvtnorm::qmvnorm(1 - alpha_sp[k],
                  tail = "lower.tail",
                  corr = Sigma0_use)$quantile
          )
        )
    }else{
      z2_search <- function(z){
        median(
          replicate(
            n.rep,
            mvtnorm::pmvnorm(
              upper = c(
                z_alpha_vec, 
                rep(z, sum(index == k))
                ), 
              corr = Sigma0_use
              ) - (1-alpha_sp[k])
            )
          )
      }
      z_alpha[k] <- uniroot(
        f = z2_search, 
        interval = c(1, z_alpha[k - 1]))$root
    }
    z_alpha_vec = c(z_alpha_vec,rep(z_alpha[k],sum(index==k)))
  }
  return(list(z_alpha = z_alpha, z_alpha_vec = z_alpha_vec))
}

#' Sample size calculation
#'
#' Sample size calculation to control the type II error or the power of an interim analysis with Maxcombo tests.
#'
#' Assume that there are 2 stages (1 interm, 1 final), and two tests for a max-combo in each stage, then we have 4 test statistics, and the two cutoff values for the two stages have been determined by \code{Maxcombo.bd} in advance using their correlation matrix and the error spending function \eqn{\alpha_1, \alpha}. The goal of this function is to control the sample size n (number of patients for both arms) or d (observed events) to achieve the ideal type II error \eqn{\beta} or the power (\eqn{1-\beta}), i.e. \eqn{\P(Z_{11}<z_1,Z_{12}<z_1,Z_{21}<z_2,Z_{22}<z_2)=\beta}.
#'
#'
#' @param Sigma1 the correlation matrix under the alternative hypothesis.
#' @param mu1 the unit mu under the alternative hypothesis (the mean of the expectation of each subject scaled weighted log-rank test statistic, which can be approximated using the formula for \eqn{E^*} in Hasegawa 2014 paper. ).
#' @param z_alpha_vec same as the one exported from Maxcombo.bd, which is the boundaries for ordered test statistics, its order should be consistent to the rows and columns in \code{Sigma1}.
#' @param beta type II error.
#' @param interim_vec the vector of the interims in each stages, not that it should be a repeat vector with same iterim values for all the test statitics at same stages.
#' @param R end of the enrollment time, which is identical to \code{R} defined in other functions like \code{\link{I.1}}.
#' @param n_range the range ot the expected patient numbers.
#' @param sum_D same as the exported value from \code{\link{sample.size_FH}}, the summed \eqn{D^*} in Hasegawa (2014).
#' @param n.rep number of repeats to take the median for output
#' @return
#' \item{n}{the number of patients needed for the trial to achieve the predefined power.}
#' \item{d}{the number of events needed for the trial to achieve the predefined power.}
#' \item{sum_D}{the input \code{sum_D} value. }
#' @author Lili Wang
#' @references Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.
#' @seealso \code{\link{Maxcombo.beta.n}}
#' @examples
#' \dontrun{
#' # install.packages("mvtnorm")
#' library(mvtnorm)
#' # install.packages("gsDesign")
#' library(gsDesign)
#' alpha <- 0.025
#' beta <- 0.1
#' # If there are two stages (K = 2), with on interim stage and a final stage
#' # First we obtain the errors spent at each stage to be identical to the ones
#'  from regular interim analysis assuming that the interim stage happened at
#'   60% of events have been observed. The error spending function used below
#'    is O\'Brien-Fleming.
#' x <- gsDesign::gsDesign(
#' k = 2, 
#' test.type = 1, 
#' timing = 0.6, 
#' sfu = "OF", 
#' alpha = alpha, 
#' beta = beta,
#' delta = -log(0.7)
#' )
#' (z <- x$upper$bound)
#' x
#' Sigma0_v <- rep(0.5,6)
#' Sigma0 <- matrix(1, ncol = 4, nrow = 4)
#' Sigma0[upper.tri(Sigma0)]<- Sigma0_v
#' Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(t(Sigma0))]
#' Sigma0
#' alpha_interim <- pnorm(z[1],lower.tail = F) # The error you would like to spend at the interim stage
#' zz <- Maxcombo.bd(
#' Sigma0 = Sigma0,
#' index = c(1, 1, 2, 2),
#' alpha_sp = c(alpha_interim,alpha)
#' )
#' zz$z_alpha # boundary value for each stage
#' zz$z_alpha_vec # boundary value for each test statistic correponding to index
#' # Correlation matrix under the alternative hypothesis
#' Sigma1_v<-rep(0.5,6)
#' Sigma1<-matrix(1, ncol=4,nrow=4)
#' Sigma1[upper.tri(Sigma1)]<- Sigma1_v
#' Sigma1[lower.tri(Sigma1)]<- t(Sigma1)[lower.tri(t(Sigma1))]
#' Sigma1
#' # Define mu1
#' mu1=c(0.1,0.1,0.2,0.2)
#' # Obtain the sample size
#' Maxcombo.sz(
#' Sigma1 = Sigma1,
#' mu1 = mu1,
#' z_alpha_vec = zz$z_alpha_vec,
#' beta = 0.1,
#' interim_vec=c(10,10,18,18),
#' R = 14,
#' n_range = c(100,1000),
#' sum_D = 0.6)
#' # need 232 patients, 140 deaths
#' }
#'
#'
#'
Maxcombo.sz <- function(
  Sigma1, 
  mu1, 
  z_alpha_vec, 
  beta, 
  interim_vec, 
  R, 
  n_range, 
  sum_D, 
  n.rep = 5){
search_n <- function(n){
  median(
    replicate(
      n.rep,
      mvtnorm::pmvnorm(
        cor = Sigma1, 
        upper = z_alpha_vec - sqrt(n * pmin(interim_vec / R, 1)) * mu1) - (beta)
      )
    )
}
n = ceiling(uniroot(f = search_n, interval = n_range)$root)
d = ceiling(n * sum_D)
return(list(n = n, d = d, sum_D = sum_D))
}

#' The power for a vector of sample sizes
#'
#' To obtain a spectrum of power for a vector of numbers of subjects (n) using \code{Maxcombo.beta.n} or events (d) using \code{Maxcombo.beta.d}. 
#'
#' @inheritParams Maxcombo.sz
#' @param n_seq the sequence of number of patients.
#' @param d_seq the sequence of number of expected events.
#' @param n.rep number of repeats to take the median for output
#' @author Lili Wang
#' @seealso \code{\link{Maxcombo.sz}}
#' @examples
#' Sigma0_v<-rep(0.5,6)
#' Sigma0<-matrix(1, ncol=4,nrow=4)
#' Sigma0[upper.tri(Sigma0)]<- Sigma0_v
#' Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(Sigma0)]
#' Sigma0
#' alpha_stage <- c(0.01,0.025) # The error you would like to spend at the interim stage
#' zz <- Maxcombo.bd(Sigma0 = Sigma0,index=c(1,1,2,2),alpha_sp=alpha_stage)
#' zz$z_alpha # boundary value for each stage
#' zz$z_alpha_vec # boundary value for each test statistic corresponding 
#' to the index
Maxcombo.beta.n <- function(
  Sigma1, 
  mu1, 
  z_alpha_vec, 
  interim_vec, 
  R,
  n_seq, 
  n.rep = 5
  ){
  sapply(
    n_seq, 
    function(n){
      median(
        replicate(
          n.rep,
          mvtnorm::pmvnorm(
            cor = Sigma1, 
            upper = z_alpha_vec - sqrt(n * pmin(interim_vec / R, 1)) * mu1
          )[[1]]))
      }
    )
}
#' @rdname Maxcombo.beta.n
Maxcombo.beta.d <- function(
  Sigma1,
  mu1,
  z_alpha_vec,
  interim_vec,
  R,
  d_seq,
  sum_D, 
  n.rep = 5
  ){
  sapply(d_seq, function(d){
    median(
      replicate(
        n.rep,
        mvtnorm::pmvnorm(
          cor = Sigma1,
          upper = z_alpha_vec - 
            sqrt(d / sum_D * pmin(interim_vec / R, 1)) * mu1)[[1]])
      )
    }
  )
}




#' A stochastic prediction results
#'
#' A stochastic-process way of prediction of the expected event ratio (\eqn{D}), mean difference (\eqn{\mu}), and the information(variance) using \code{stoch_pred} or the covariance using \code{stoch_pred.cov}.
#'
#' @param eps delayed treatment effect time.
#' @param p probability of treatment assignment.
#' @param b the number of sub-intervals at each time point, the larger the finer splitting for more accurate computation. Usually \eqn{b = 30} is sufficient.
#' @param omega the minimum follow-up time for all the patients. Note that Hasegawa(2014) assumes that the accrual is uniform between time 0 and R, and there does not exist any censoring except for the administrative censoring at the ending time \eqn{\tau}. Thus this value omega is equivalent to \code{tau-R}. 
#' @param lambda the hazard for the control group.
#' @param theta the hazard ratio after the delayed time \code{eps} for the treatment arm.
#' @param rho,rho1,rho2 the first parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.
#' @param gamma,gamma1,gamma2 the second parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.
#' @param R the accrual period. 
#' 
#' @return 
#' \item{sum_D}{the mean expected event ratio. Once being multiplied by \code{n}, it will become the stochastically predicted event size. }
#' \item{inf or covariance}{the information/variance or covariance (averaged for each subject), should be multiplied by \code{n}, which gives the stochastically predicted information. }
#' \item{E.star}{the unit mean, corresponding to \eqn{E^*} in Hasegawa(2014), or the \eqn{\tilde{\mu}} of formula (8) in Wang et al(2021).}
#' \item{trt_vs_ctrl_N}{the ratio of the samples sizes between the two arms, treatment vs control, corresponding to the time vector \code{t_vec}.}
#' \item{t_vec}{the time sequence corresponding to \code{trt_vs_ctrl_N}.}
#' @author Lili Wang
#' @references 
#' Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135. 
#' Wang, L., Luo, X., & Zheng, C. (2021). A Simulation-free Group Sequential Design with Max-combo Tests in the Presence of Non-proportional Hazards. Journal of Pharmaceutical Statistics.
#' 
stoch_pred <- function(eps, p, b, tau, omega, lambda, theta, rho, gamma, R){
  n_sub <- floor(b * tau)
  t <- c(0, seq(1, n_sub) / b)
  h_1 <- rep(lambda, (n_sub + 1)) #control
  h_2 <- c(
    rep(lambda, min(round(eps * b), n_sub + 1)),
    rep(lambda * theta, max(n_sub - round(eps * b) + 1, 0))
    ) #treatment
  N_1 <- rep((1 - p),(n_sub + 1))
  N_2 <- rep(p, (n_sub + 1))
  for(i in 1:max(n_sub - 1,1)){
    N_1[i + 1] <- N_1[i] * (1 - h_1[i] / b-(t[i] > omega) / b / (tau - t[i]))
    N_2[i + 1] <- N_2[i] * (1 - h_2[i] / b-(t[i] > omega) / b / (tau - t[i]))
  }
  N_1[n_sub + 1] <- N_2[n_sub + 1] <- 0
  
  f_S_1 <- function(x) exp( - lambda * x)
  f_S_2 <- function(x) (x < eps) * exp( - lambda * x) +
    (x >= eps) * exp( - (lambda * eps + lambda * theta * (x - eps)))
  #f_S_2_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-eps*lambda.trt*(1/theta-1))*exp(-lambda.trt*x)
  S_1 <- f_S_1(t)
  S_2 <- f_S_2(t)
  S <- (1 - p) * S_1 + p * S_2
  D <- (h_1 * N_1 + h_2 * N_2) / b * min(tau / R, 1) #predicted total events at each time
  theta_seq <- h_2 / h_1
  phi <- N_2 / N_1
  r <- S^rho * (1 - S)^gamma
  num_vec <- D * r * (phi * theta_seq/(1 + phi * theta_seq)-phi / (1 + phi))
  den_vec <- D * r^2 * phi / (1 + phi)^2
  E.star_num <- sum(num_vec[1:n_sub])
  E.star_den <- sqrt(sum(den_vec[1:n_sub]))
  E.star <- E.star_num / E.star_den
  #time_vec=seq(1,n_sub)/b)
  return(
    list(
      sum_D = sum(D[1:n_sub]),
      inf = sum(den_vec[1:n_sub]), 
      E.star = E.star, 
      trt_vs_ctrl_N = phi[1:n_sub], 
      t_vec = seq(1, n_sub) / b
      )
    )
}
#' 
#' 
#' @rdname stoch_pred
#' 
stoch_pred.cov<-function(
  eps, 
  p, 
  b, 
  tau, 
  omega, 
  lambda, 
  theta, 
  rho1, 
  gamma1, 
  rho2, 
  gamma2, 
  R
  ){
  n_sub <- floor(b * tau)
  t <- c(0, seq(1, n_sub) / b)
  h_1 <- rep(lambda, (n_sub + 1)) #control
  h_2 <- c(
    rep(lambda, min(eps * b,n_sub + 1)), 
    rep(lambda * theta, max(n_sub - eps * b+1, 0))
    ) #treatment
  N_1 <- rep((1 - p),(n_sub + 1))
  N_2 <- rep(p, (n_sub + 1))
  for(i in 1:(n_sub - 1)){
    N_1[i + 1] <- N_1[i] * (1 - h_1[i] / b - (t[i] > omega) / b / (tau - t[i]))
    N_2[i + 1] <- N_2[i] * (1 - h_2[i] / b - (t[i] > omega) / b / (tau - t[i]))
  }
  N_1[n_sub + 1] <- N_2[n_sub + 1] <- 0
  
  f_S_1 <- function(x) exp( - lambda * x)
  f_S_2 <- function(x) (x < eps) * exp( - lambda * x)+
                     (x >= eps) * exp( - (lambda * eps + lambda * theta * (x - eps)))
  #f_S_2_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-eps*lambda.trt*(1/theta-1))*exp(-lambda.trt*x)
  S_1 <- f_S_1(t)
  S_2 <- f_S_2(t)
  S <- (1 - p) * S_1 + p * S_2
  #predicted total events at each time
  D <- (h_1 * N_1 + h_2 * N_2) / b * min(tau / R, 1) 
  
  theta_seq <- h_2 / h_1
  phi <- N_2 / N_1
  r1 <- S^rho1 * (1 - S)^gamma1
  r2 <- S^rho2 * (1 - S)^gamma2
  den_vec <- D * r1 * r2 * phi / (1 + phi)^2
  #time_vec=seq(1,n_sub)/b)
  return(
    list(
      sum_D = sum(D[1:n_sub]),
      cov = sum(den_vec[1:n_sub]), 
      trt_vs_ctrl_N = phi[1:n_sub], 
      t_vec = seq(1, n_sub) / b
      )
    )
}

#' Predicted sample sizes and boundaries for GSMC design 
#' 
#' Compute predicted sample size and boundaries for a group sequential design of max-combo tests.
#' 
#' Predict the sample sizes and boundaries to achieve the targeted type I errors (error spent at each stage) and power. Prediction approaches include the exact prediction or the stochastic prediction approach following 2-piece-wise exponential distributions given in the appendix of the reference paper. 
#' 
#' @param eps the change point, before which, the hazard ratio is 1, and after which, the hazard ratio is theta
#' @param p treatment assignment probability.
#' @param b the number of subintervals per time unit.
#' @param tau the end of the follow-up time  in the study. Note that this is identical to \eqn{T+\tau} in the paper from Hasegawa (2014).
#' @param omega the minimum follow-up time for all the patients.  Note that Hasegawa(2014) assumes that the accrual is uniform between time 0 and R, and there does not exist any censoring except for the administrative censoring at the ending time \eqn{\tau}. Thus this value omega is equivalent to \code{tau-R}. Through our simulation tests, however, we found that this function is quite robust to violations of these assumptions: dropouts, different censoring rates for two  arms, and changing accrual rates. 
#' @param  FHweights a list of pairs of parameters (rho, gamma) for the Fleming Harrington weighted log-rank tests:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}. The first one will provide the information fraction we will use to sign the stopping time. Now we only accommodate the surrogate information fraction, i.e., rho = 0 and gamma = 0.
#' @param interim_ratio a vector of ratios between 0 and 1 that sign the stopping time following the first type of information fraction in FHweights, i.e. the event sizes. Note that the last value must be 1, otherwise, the function will automatically append 1 for the final stage.
#' @param error_spend cumulative errors spent at each stage (including the final one). Must be of the same length as interim_ratio or one element shorter if the last value of interim_ratio is less than 1. The last element of error_spend should be type I error \code{alpha}.
#' @param lambda the hazard for the control group.
#' @param lambda.trt the hazard for the treatment group after time eps.
#' @param beta type II error, 1-power. 
#' @param stoch stochastic prediction or not (exact prediction), the former is default. 
#' @param range_ext the width to extend the range of sample sizes. If the predicted sample size is not found in the range, try a wider \code{range_ext}. It will automatically search for the likely range, but still possibly miss the best one. Its default value is 200. 
#' @param time_ext_multiplier compute the time window for the possible stopping time points by multiplying it with the total expected follow-up time \code{tau}. The default is 1.5. In other words, when tau = 18, the longest time we consider would be \eqn{18 * 1.5 = 27} months.
#' @param time_increment time increments to compute the predicted stopping time points, the finer the more accurate. 
#' @param n.rep same as the one in \code{\link{Maxcombo.sz}}. The number of repeats to take the median for output since the called likelihood generator of a multivariate normal distribution \code{\link[mvtnorm]{pmvnorm}} is not determinant. The default \code{n.rep} value is 5. 
#' @return 
#' \item{z_alpha_pred}{predicted boundary values for all the stages, length is equivalent to the input \code{interim_ratio} or \code{error_spend}.}
#' \item{z_alpha_vec_pred}{predicted boundary values for all the test statistics following \code{index}.}
#' \item{d_fixed}{the required observed events at each stage. }
#' \item{n_FH}{the total required number of subjects according to the defined hazards.}
#' \item{n_event_FH}{the total required number of events according to the defined hazards.}
#' \item{index}{records the ordered stages for each tests, starting from 1 and ending at the length of \code{interim_ratio} or \code{error_spend}. It is actually \code{rep(1:length(interim_ratio), each = length(FHweights))}. }
#' \item{interim_pred0}{predicts stopping time points under the null hypothesis following the order of \code{index}.}
#' \item{interim_pred1}{predicts stopping time points under the alternative hypothesis following the order of \code{index}.}
#' \item{Sigma0}{predicts correlation matrix under the null hypothesis with each row and column following the test statistics corresponding to \code{index}. }
#' \item{Sigma1}{predicts correlation matrix under the alternative hypothesis with each row and column following the test statistics corresponding to \code{index}. }
#' \item{mu1}{the predicts unit mean under the alternative hypothesis, the \eqn{\tilde{\mu}} in formula (5) of the reference paper. The test statistics follow \code{index}.  It is also the mean of the expectation of each subject scaled weighted log-rank test statistic, which can be approximated using the formula for \eqn{E^*} in Hasegawa 2014 paper. Under null, the predicted mean is otherwise 0, implying no treatment effect.  }
#' \item{stoch}{input \code{stoch} boolean variable, \code{TRUE} if stochastic prediction is enabled, \code{FALSE} otherwise. The default is \code{TRUE}. }
#' \item{FHweights}{input \code{FHweights} list. }
#' \item{interim_ratio}{input \code{interim_ratio} vector.}
#' \item{error_spend}{input \code{error_spend} vector. }
#' 
#' @author Lili Wang
#' @references 
#' Wang, L., Luo, X., & Zheng, C. (2021). A Simulation-free Group Sequential Design with Max-combo Tests in the Presence of Non-proportional Hazards. Journal of Pharmaceutical Statistics.
#' @examples 
#' \dontrun{
#' ### Parameters
#' FHweights <- list(
#'   c(0,0),
#'   c(0,1),
#'   c(1,0)
#'   )
#'   n_FHweights <- length(FHweights)
#'   # stop when what proportion of the events have been observed.
#'   fixed_death_p <-  c(0.6, 0.7, 0.8)
#'   interim_ratio <- c(fixed_death_p,1)
#'   n_stage <- length(interim_ratio)
#'   # treatment assignment.
#'   p <- 1/2 
#'   # end of the study, chronological assuming it's on the alternative arm.
#'   tau <- 18
#'   # end of the accrual period, chronological.
#'   R <- 14 
#'   # minimum potential follow-up time, not considering dropouts or other subject-specific censoring.
#'   omega <- (tau-R) 
#'   # waiting time before the change point: can be the delayed effect or the crossing effect
#'   eps <- 2
#'   # event hazard of the control arm.
#'   lambda <- log(2)/6 
#'   # hazard ratio after the change point eps under the alternative hypothesis. 
#'   theta <- 0.6 
#'   # event hazard for the treatment arm under the alternative hypothesis and after the change point eps.
#'   lambda.trt <- lambda*theta  
#'   # type I error under the control.
#'   alpha <- 0.025 
#'   # type II error under the alternative. 
#'   beta <- 0.1 
#'   # Obtain the cumulative errors spent at each stage
#'   error_spend <- c(0.005, 0.01, 0.015, alpha)
#'   # number of subintervals for a time unit
#'   b <- 30
#'   res <- GSMC_design(
#'   FHweights,
#'   interim_ratio,
#'   error_spend,
#'   eps, 
#'   p, 
#'   b, 
#'   tau,]
#'   omega,
#'   lambda,
#'   lambda.trt,
#'   rho, 
#'   gamma,
#'   beta, 
#'   stoch = F
#'   )
#' }
#' 
GSMC_design <- function(
  FHweights,
  interim_ratio,
  error_spend,
  eps, 
  p, 
  b, 
  tau, 
  omega,
  lambda,
  lambda.trt,
  rho, 
  gamma,
  beta,
  stoch = TRUE,
  range_ext = 200,
  time_ext_multiplier = 1.5,
  time_increment = 0.01,
  n.rep = 5
){
  n_stage <- length(error_spend)
  n_stage_ratio <- length(interim_ratio)
  if(interim_ratio[n_stage_ratio] < 1){
    n_stage_ratio <- n_stage_ratio + 1
    interim_ratio <- c(interim_ratio, 1)
  }
  if(n_stage_ratio != n_stage){
    stop("If the last interim_ratio is not 1, it will be automatically appended by 1 to become c(interim_ratio, 1). The final length of interim_ratio must be identical to that of error_spend. ")
  }
  n_FHweights <- length(FHweights)
  if(FHweights[[1]][1]!=0 || FHweights[[1]][2]!=0){
    stop("The first weight must be the standard log-rank test: rho = gamma = 0. ")
  } 
  if(!all(unlist(FHweights)%in% c(0,1))){
    stop("Parameters in FHweights must be either 0 or 1. ")
  } 
  alpha <- error_spend[n_stage]
  size_FH_ls <- lapply(
    FHweights, 
    function(weight)
      sample.size_FH(
        eps, 
        p, 
        b, 
        tau, 
        omega, 
        lambda,
        lambda.trt,
        weight[1],
        weight[2],
        alpha,
        beta)
  )
  
  sum_D <- size_FH_ls[[1]]$sum_D # unit death, should be identical to size_FH_lower$sum_D
  size_pre_vec <- sapply(size_FH_ls, '[[','n')
  # size_seq is the range of the possible sizes using the old method calculation. 
  size_range <- pmax(1, range(size_pre_vec) + 
                       c(-range_ext, +range_ext))
  
  size_seq <- size_range[1]:size_range[2]
  
  t_ls <- seq(
    # ensure that it starts from at least with 2 subintervals in stoch_pred
    2/b, 
    ceiling(tau * time_ext_multiplier),
    time_increment
    )
  
  if(stoch){
    # Under the null, predict the interim and final stage time. 
    interim_pred0 <- lapply( 
      interim_ratio,
      function(prop){
        t_ls[
          which.min(
            abs(
              sum_D * prop - 
                sapply(t_ls,
                       function(t){
                         stoch_pred(
                           eps,
                           p,
                           b,
                           t,
                           omega = max(t - R, 0),
                           lambda,
                           1,
                           0,
                           0,
                           R)$sum_D})))]
      })
    
    
    # Under the alternative, predict the interim and final stage time.
    interim_pred1 <- lapply( 
      interim_ratio,
      function(prop){
        t_ls[
          which.min(
            abs(
              sum_D * prop - 
                sapply(t_ls, function(t){
                  stoch_pred(
                    eps,
                    p,
                    b,
                    t,
                    omega = max(t - R,0),
                    lambda, 
                    theta, 
                    0,
                    0,
                    R)$sum_D})))]
      })
    
    nvalues <- n_stage * n_FHweights
    Sigma1 <- Sigma0 <- matrix(NA,
                               nrow = nvalues,
                               ncol = nvalues)
    # Predict the correlation values under the null (indexed by 0) and the alternative (indexed by 1). 
    
    for(i in 1:nvalues){
      for(j in i:nvalues){
        i_test <- ifelse(i %% n_FHweights == 0, 
                         n_FHweights,
                         i %% n_FHweights )
        i_stage <- ceiling(i / n_FHweights)
        j_test <- ifelse(j %% n_FHweights == 0, 
                         n_FHweights, 
                         j %% n_FHweights )
        j_stage <- ceiling(j / n_FHweights)  
        if(i_test == j_test && i_stage != j_stage){
          Sigma0[i,j] <- sqrt(
            stoch_pred(
              eps,
              p,
              b,
              interim_pred0[[i_stage]],
              omega = max(interim_pred0[[i_stage]]-R, 0),
              lambda,
              1,
              FHweights[[i_test]][1],
              FHweights[[i_test]][2],
              R)$inf / 
              stoch_pred(
                eps,
                p,
                b,
                interim_pred0[[j_stage]],
                omega = max(interim_pred0[[j_stage]]-R, 0),
                lambda,
                1,
                FHweights[[j_test]][1],
                FHweights[[j_test]][2],
                R)$inf
          )
          
          Sigma1[i,j] <- sqrt(
            stoch_pred(
              eps,
              p,
              b,
              interim_pred1[[i_stage]],
              omega = max(interim_pred1[[i_stage]]-R, 0),
              lambda,
              theta,
              FHweights[[i_test]][1],
              FHweights[[i_test]][2],
              R)$inf / 
              stoch_pred(
                eps,
                p,
                b,
                interim_pred1[[j_stage]],
                omega = max(interim_pred1[[j_stage]]-R, 0),
                lambda,
                theta,
                FHweights[[j_test]][1],
                FHweights[[j_test]][2],
                R)$inf
          )
          
        }else if(i_test != j_test && i_stage == j_stage){
          Sigma0[i,j] <- stoch_pred.cov(
            eps,
            p,
            b,
            interim_pred0[[i_stage]],
            omega = max(interim_pred0[[i_stage]] - R, 0),
            lambda,
            1,
            FHweights[[i_test]][1],
            FHweights[[i_test]][2],
            FHweights[[j_test]][1],
            FHweights[[j_test]][2],
            R)$cov / sqrt(
              stoch_pred(
                eps,
                p,
                b,
                interim_pred0[[i_stage]],
                omega = max(interim_pred0[[i_stage]] - R, 0),
                lambda, 
                1,
                FHweights[[i_test]][1],
                FHweights[[i_test]][2],
                R)$inf * stoch_pred(
                  eps,
                  p,
                  b,
                  interim_pred0[[i_stage]],
                  omega = max(interim_pred0[[i_stage]] - R, 0),
                  lambda,
                  1,
                  FHweights[[j_test]][1],
                  FHweights[[j_test]][2],
                  R)$inf
            ) 
          
          Sigma1[i,j] <- stoch_pred.cov(
            eps,
            p,
            b,
            interim_pred1[[i_stage]],
            omega = max(interim_pred1[[i_stage]] - R, 0),
            lambda,
            theta,
            FHweights[[i_test]][1],
            FHweights[[i_test]][2],
            FHweights[[j_test]][1],
            FHweights[[j_test]][2],
            R)$cov / sqrt(
              stoch_pred(
                eps,
                p,
                b,
                interim_pred1[[i_stage]],
                omega = max(interim_pred1[[i_stage]] - R, 0),
                lambda, 
                theta,
                FHweights[[i_test]][1],
                FHweights[[i_test]][2],
                R)$inf * stoch_pred(
                  eps,
                  p,
                  b,
                  interim_pred1[[i_stage]],
                  omega = max(interim_pred1[[i_stage]] - R, 0),
                  lambda,
                  theta,
                  FHweights[[j_test]][1],
                  FHweights[[j_test]][2],
                  R)$inf
            ) 
          
        }else if (i_test != j_test && i_stage != j_stage){
          if(is.na(
            Sigma0[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ])){
            Sigma0[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ] <- stoch_pred.cov(
              eps,
              p,
              b,
              interim_pred0[[i_stage]],
              omega = max(interim_pred0[[i_stage]] - R, 0),
              lambda,
              1,
              FHweights[[i_test]][1],
              FHweights[[i_test]][2],
              FHweights[[j_test]][1],
              FHweights[[j_test]][2],
              R)$cov / sqrt(
                stoch_pred(
                  eps,
                  p,
                  b,
                  interim_pred0[[i_stage]],
                  omega = max(interim_pred0[[i_stage]] - R, 0),
                  lambda, 
                  1,
                  FHweights[[i_test]][1],
                  FHweights[[i_test]][2],
                  R)$inf * stoch_pred(
                    eps,
                    p,
                    b,
                    interim_pred0[[i_stage]],
                    omega = max(interim_pred0[[i_stage]] - R, 0),
                    lambda,
                    1,
                    FHweights[[j_test]][1],
                    FHweights[[j_test]][2],
                    R)$inf
              ) 
          } 
          if(is.na(
            Sigma1[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ])){
            Sigma1[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ] <- stoch_pred.cov(
              eps,
              p,
              b,
              interim_pred1[[i_stage]],
              omega = max(interim_pred1[[i_stage]] - R, 0),
              lambda,
              theta,
              FHweights[[i_test]][1],
              FHweights[[i_test]][2],
              FHweights[[j_test]][1],
              FHweights[[j_test]][2],
              R)$cov / sqrt(
                stoch_pred(
                  eps,
                  p,
                  b,
                  interim_pred1[[i_stage]],
                  omega = max(interim_pred1[[i_stage]] - R, 0),
                  lambda, 
                  theta,
                  FHweights[[i_test]][1],
                  FHweights[[i_test]][2],
                  R)$inf * stoch_pred(
                    eps,
                    p,
                    b,
                    interim_pred1[[i_stage]],
                    omega = max(interim_pred1[[i_stage]] - R, 0),
                    lambda,
                    theta,
                    FHweights[[j_test]][1],
                    FHweights[[j_test]][2],
                    R)$inf
              ) 
          } 
          if(is.na(
            Sigma0[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ])){
            Sigma0[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ] <- sqrt(
              stoch_pred(
                eps,
                p,
                b,
                interim_pred0[[i_stage]],
                omega = max(interim_pred0[[i_stage]]-R, 0),
                lambda,
                1,
                FHweights[[j_test]][1],
                FHweights[[j_test]][2],
                R)$inf / 
                stoch_pred(
                  eps,
                  p,
                  b,
                  interim_pred0[[j_stage]],
                  omega = max(interim_pred0[[j_stage]]-R, 0),
                  lambda,
                  1,
                  FHweights[[j_test]][1],
                  FHweights[[j_test]][2],
                  R)$inf
            )
            
          }
          
          if(is.na(
            Sigma1[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ])){
            Sigma1[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ] <- sqrt(
              stoch_pred(
                eps,
                p,
                b,
                interim_pred1[[i_stage]],
                omega = max(interim_pred1[[i_stage]]-R, 0),
                lambda,
                theta,
                FHweights[[j_test]][1],
                FHweights[[j_test]][2],
                R)$inf / 
                stoch_pred(
                  eps,
                  p,
                  b,
                  interim_pred1[[j_stage]],
                  omega = max(interim_pred1[[j_stage]]-R, 0),
                  lambda,
                  theta,
                  FHweights[[j_test]][1],
                  FHweights[[j_test]][2],
                  R)$inf
            )
            
          }
          
          
          Sigma0[i,j] <- Sigma0[
            i_test + (i_stage - 1) * n_FHweights,
            j_test + (i_stage - 1) * n_FHweights  
          ] * Sigma0[
            j_test + (i_stage-1) * n_FHweights,
            j_test + (j_stage-1) * n_FHweights 
          ]
          
          Sigma1[i,j] <- Sigma1[
            i_test + (i_stage - 1) * n_FHweights,
            j_test + (i_stage - 1) * n_FHweights  
          ] * Sigma1[
            j_test + (i_stage-1) * n_FHweights,
            j_test + (j_stage-1) * n_FHweights 
          ]
          
        }else{
          Sigma0[i,j] = 1
          Sigma1[i,j] = 1
        }
        
      }
    }
    #identical(Sigma0[upper.tri(Sigma0)],  Sigma0[upper.tri(t(Sigma0))])
    Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(t(Sigma0))]
    Sigma1[lower.tri(Sigma1)]<- t(Sigma1)[lower.tri(t(Sigma1))]  
    
  }else{
    # Exact prediction 
    # Event time
    interim_pred0 <- lapply( 
      interim_ratio,
      function(prop){
        t_ls[
          which.min(
            abs(
              sum_D * prop - 
                sapply(t_ls,
                       function(t) {
                         h.tilde(0, lambda, 1, eps, R, p, t)})))]
      })
    
    interim_pred1 <- lapply( 
      interim_ratio,
      function(prop){
        t_ls[
          which.min(
            abs(
              sum_D * prop - 
                sapply(t_ls, function(t) {
                  h.tilde(0, lambda, theta, eps, R, p, t)})))]
      })
    
    nvalues <- n_stage * n_FHweights
    Sigma1 <- Sigma0 <- matrix(NA, nrow = nvalues, ncol = nvalues)
    # Predict the correlation values under the null (indexed by 0) and the alternative (indexed by 1). 
    
    for(i in 1:nvalues){
      for(j in i:nvalues){
        i_test <- ifelse(i %% n_FHweights == 0, 
                         n_FHweights, 
                         i %% n_FHweights )
        i_stage <- ceiling(i / n_FHweights)
        j_test <- ifelse(j %% n_FHweights == 0, 
                         n_FHweights, 
                         j %% n_FHweights )
        j_stage <- ceiling(j / n_FHweights) 
        if(i_test == j_test && i_stage != j_stage){
          Sigma0[i,j] <-  sqrt(
            I.0(
              FHweights[[i_test]][1],
              FHweights[[i_test]][2],
              lambda = lambda,
              R = R,
              p = p,
              t.star =  interim_pred0[[i_stage]]) / I.0(
                FHweights[[i_test]][1],
                FHweights[[i_test]][2],
                lambda = lambda,
                R = R,
                p = p,
                t.star =  interim_pred0[[j_stage]]
              ))
          
          Sigma1[i,j] <-  sqrt(
            I.1(
              FHweights[[i_test]][1],
              FHweights[[i_test]][2],
              lambda = lambda,
              theta = theta,
              eps = eps,
              R = R,
              p = p,
              t.star = interim_pred1[[i_stage]]) / I.1(
                FHweights[[i_test]][1],
                FHweights[[i_test]][2],
                lambda = lambda, 
                theta = theta,
                eps = eps,
                R = R,
                p = p,
                t.star = interim_pred1[[j_stage]]
              ))
          
          
        }else if(i_test != j_test && i_stage == j_stage){
          Sigma0[i,j] <- cor.0(
            rho1 = FHweights[[i_test]][1],
            gamma1 = FHweights[[i_test]][2],
            rho2= FHweights[[j_test]][1],
            gamma2= FHweights[[j_test]][2],
            lambda = lambda, 
            R = R,
            p = p,
            t.star = interim_pred0[[i_stage]]
          )
          
          Sigma1[i,j] <- cor.1(
            rho1 = FHweights[[i_test]][1],
            gamma1 = FHweights[[i_test]][2],
            rho2 = FHweights[[j_test]][1],
            gamma2 = FHweights[[j_test]][2],
            lambda = lambda,
            theta = theta,
            eps = eps,
            R = R,
            p = p,
            t.star = interim_pred1[[i_stage]]
          )
          
        }else if (i_test != j_test && i_stage != j_stage){
          if(is.na(
            Sigma0[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ])){
            Sigma0[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ] <- cor.0(
              rho1 = FHweights[[i_test]][1],
              gamma1 = FHweights[[i_test]][2],
              rho2= FHweights[[j_test]][1],
              gamma2= FHweights[[j_test]][2],
              lambda = lambda, 
              R = R,
              p = p,
              t.star = interim_pred0[[i_stage]]
            )
          } 
          if(is.na(
            Sigma1[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ])){
            Sigma1[
              i_test + (i_stage - 1) * n_FHweights,
              j_test + (i_stage - 1) * n_FHweights     
            ] <- cor.1(
              rho1 = FHweights[[i_test]][1],
              gamma1 = FHweights[[i_test]][2],
              rho2 = FHweights[[j_test]][1],
              gamma2 = FHweights[[j_test]][2],
              lambda = lambda,
              theta = theta,
              eps = eps,
              R = R,
              p = p,
              t.star = interim_pred1[[i_stage]]
            )
          } 
          if(is.na(
            Sigma0[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ])){
            Sigma0[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ] <- sqrt(
              I.0(
                FHweights[[j_test]][1],
                FHweights[[j_test]][2],
                lambda = lambda,
                R = R,
                p = p,
                t.star =  interim_pred0[[i_stage]]) / I.0(
                  FHweights[[j_test]][1],
                  FHweights[[j_test]][2],
                  lambda = lambda,
                  R = R,
                  p = p,
                  t.star =  interim_pred0[[j_stage]]
                ))
          }
          
          if(is.na(
            Sigma1[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ])){
            Sigma1[
              j_test + (i_stage-1) * n_FHweights,
              j_test + (j_stage-1) * n_FHweights     
            ] <- sqrt(
              I.1(
                FHweights[[j_test]][1],
                FHweights[[j_test]][2],
                lambda = lambda,
                theta = theta,
                eps = eps,
                R = R,
                p = p,
                t.star = interim_pred1[[i_stage]]) / I.1(
                  FHweights[[j_test]][1],
                  FHweights[[j_test]][2],
                  lambda = lambda, 
                  theta = theta,
                  eps = eps,
                  R = R,
                  p = p,
                  t.star = interim_pred1[[j_stage]]
                ))
            
          }
          
          
          Sigma0[i,j] <- Sigma0[
            i_test + (i_stage - 1) * n_FHweights,
            j_test + (i_stage - 1) * n_FHweights  
          ] * Sigma0[
            j_test + (i_stage-1) * n_FHweights,
            j_test + (j_stage-1) * n_FHweights 
          ]
          
          Sigma1[i,j] <- Sigma1[
            i_test + (i_stage - 1) * n_FHweights,
            j_test + (i_stage - 1) * n_FHweights  
          ] * Sigma1[
            j_test + (i_stage-1) * n_FHweights,
            j_test + (j_stage-1) * n_FHweights 
          ]
          
        }else{
          Sigma0[i,j] = 1
          Sigma1[i,j] = 1
        }
        
      }
    }
    
    #identical(Sigma0[upper.tri(Sigma0)],  Sigma0[upper.tri(t(Sigma0))])
    Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(Sigma0)]
    Sigma1[lower.tri(Sigma1)]<- t(Sigma1)[lower.tri(Sigma1)]
    
  }
  
  ## Obtain he predicted boundaries
  ### short z_alpha
  index <- rep(1:n_stage, each = n_FHweights)
  z_alpha_pred <- Maxcombo.bd(
    Sigma0 = Sigma0,
    index = index,
    alpha_sp = error_spend)$z_alpha
  ### long z_alpha_vec
  z_alpha_vec_pred <- rep(z_alpha_pred, each = n_FHweights)
  ### For only one final stage
  z_final_alpha_pred <- Maxcombo.bd(
    Sigma0 = Sigma0[(nvalues-n_FHweights+1):nvalues, 
                    (nvalues-n_FHweights+1):nvalues],
    index = c(1, 1),
    alpha_sp = c(alpha),
    n.rep = n.rep
    )$z_alpha
  
  ## Obtain the predicted mean: 
  ### add the "-" before the calculated E* to be consistent with the later calculations using -WLRT/SLRT other than WLRT/SLRT, to make the test statistics tend to be positive under the alternative. Note need to do so if we are not changing the sign of the WLRT/SLRTs. Just make sure that the means and the test statistics are consistent. 
  mu1 <- NULL
  k <-  0
  for(i in 1:n_stage){
    for(j in 1:n_FHweights){
      k <- k + 1
      mu1[k] <- -sample.size_FH(
        eps,
        p,
        b,
        interim_pred1[[i]],
        omega = max(interim_pred1[[i]] - R, 0),
        lambda,
        lambda.trt,
        FHweights[[j]][1], 
        FHweights[[j]][2],
        alpha,
        beta
      )$E.star
    }
  }
  
  ## Obtain the predicted sample sizes
  n_FH_ls<-Maxcombo.sz(
    Sigma1 = Sigma1,
    mu1 = mu1,
    z_alpha_vec = z_alpha_vec_pred,
    beta = beta,
    interim_vec = rep(unlist(interim_pred1), each = n_FHweights),
    R = R,
    n_range = size_range,
    sum_D = sum_D,
    n.rep = n.rep)
  
  n_FH <- n_FH_ls$n
  n_event_FH <- n_FH_ls$d
  # the number of events needed to pause the study at each stage
  d_fixed <- ceiling(interim_ratio * n_event_FH) 
 
  # cat("power")
  # 1-Maxcombo.beta.n(Sigma1,
  #                   mu1,
  #                   z_alpha_vec_pred,
  #                   rep(unlist(interim_pred1), each = n_FHweights),
  #                   R,
  #                   n_FH)
  # cat("power at interim stage(s)")
  # sapply(1:n_stage, function(stage){
  #   1-Maxcombo.beta.n(
  #     Sigma1[1:(n_FHweights*(stage)),1:(n_FHweights*(stage))],
  #     mu1[1:(n_FHweights*(stage))],
  #     z_alpha_vec_pred[1:(n_FHweights*(stage))],
  #     rep(unlist(interim_pred1),each = n_FHweights)[1:(n_FHweights*(stage))], R,
  #     n_FH)
  # })
  out <- list(
    z_alpha_pred = z_alpha_pred,
    z_alpha_vec_pred = z_alpha_vec_pred, 
    z_final_alpha_pred = z_final_alpha_pred,
    d_fixed = d_fixed,
    n_FH = n_FH,
    n_event_FH = n_event_FH,
    index = index,
    interim_pred0 = interim_pred0,
    interim_pred1 = interim_pred1,
    Sigma0 = Sigma0,
    Sigma1 = Sigma1,
    mu1 = mu1,
    FHweights = FHweights,
    interim_ratio = interim_ratio,
    error_spend = error_spend
  )
  return(out)
}
