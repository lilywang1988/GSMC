R package GSMC: a simulation-free group sequential design with max-combo
tests in the presence of non-proportional hazards
================
[Lili Wang](mailto:lilywang@umich.edu)
2021-03-07

## Purpose

This R package is to implement the proposed [simulation-free method for
group sequential design with maxcombo
tests](https://arxiv.org/abs/1911.05684 "Hey, direct me to the paper!")(GS-MC).
The goal here is to improve the detection power of the weighted log-rank
tests (WLRTs) when the non-proportional hazards are present. Moreover,
the method is simulation-free, and allows multiple checks (interims)
before ending the study.

For details about the calculation of the information and covariance
using estimation (data-driven) and prediction methods please find
another R package [IAfrac](https://github.com/lilywang1988/IAfrac). The
recent version (\>= 0.2.1) has all the functions in package `IAfrac`
imported already so there is no need to attach it separately.

We hope this package will be improved gradually. Please let us know your
feedback.

## Preparation

To install and use this R package from Github, please uncomment the
codes.

``` r
# install.packages("devtools")
# library(devtools)
# 
# Optional: to obtain the expected errors spent at each stage
# install.packages("gsDesign")
library(gsDesign) 

# Information and covariance calculation; sample size computation using the Hasegawa proposal
# install_github("lilywang1988/GSMC")
library(GSMC) 
```

## Sample size and threshold prediction for three tests and 4 stages (3 interim and 1 final)

The detailed difference between the two prediction approaches can be
found in the [methodology
paper](https://arxiv.org/abs/1911.05684 "Hey, direct me to the methodology paper!").

To sum up, the differences include: 1. The exact method is faster. 2.
The exact approach is limited to piece-wise exponential distribution. 3.
When the number of subintervals at each time unit (or value b) is large
in stochastic approach, both methods produce very comparable results.

### Stochastic prediction approach

We consider three different Fleming Harrington type weighted log-rank
tests, and **the first one must be the simple log-rank test**. We will
use the event number ratios (surrogate information fractions in the
Hasegawa 2016 paper) to monitor the stopping time points.

``` r
### Parameters
# The first pair of FHweight parameters must be a simple log-rank test
FHweights <- list(
  c(0,0), # must be c(0,0)
  c(0,1),
  c(1,0)
)
n_FHweights <- length(FHweights)

# stop when what proportions of the events have been observed.
interim_ratio <- c(0.6, 0.7, 0.8, 1)
n_stage <- length(interim_ratio)

# treatment assignment
p <- 1/2 
# end of the study, chronological, assuming it's on the alternative arm
tau <- 18
# end of the accrual period, chronological
R <- 14 
# minimum potential follow-up time
omega <- (tau-R) 
# waiting time before the delayed effect starts
eps <- 2
# event hazard of the control arm.
lambda <- log(2)/6 
# hazard ratio after the change point eps under the alternative hypothesis. 
theta <- 0.6 
# event hazard for the treatment arm under the alternative hypothesis and after the change point eps.
lambda.trt <- lambda * theta  
# cumulative type I error under the control.
alpha <- 0.025 
# cumulative type II error under the alternative. 
beta <- 0.1 

# Obtain the cumulative errors spent at each stage following O\'Brien-Fleming
x <- gsDesign(
  k = n_stage, 
  test.type = 1, 
  timing = interim_ratio[-n_stage], 
  sfu = "OF", 
  alpha = alpha, 
  beta = beta, 
  delta = -log(theta)
)
# obtain the cumulative errors spent at each stage
error_spend <- cumsum(x$upper$prob[,1])
# correct the last entry
error_spend[length(error_spend)] <- alpha 


# number of subintervals for a time unit, the larger the finer 
# and more accurate; 30 should be good enough for most cases
b <- 30

# Obtain the settings from stochastic prediction
setting_stoch <- GSMC_design(
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
  stoch = T
)
# Boundaries at each stage
setting_stoch$z_alpha_pred
#> [1] 2.858002 2.698675 2.548058 2.316725
# Event numbers to pause/stop
setting_stoch$d_fixed
#> [1] 189 221 252 315

# Boundary if only consider the last stage (without interim stages)
setting_stoch$z_final_alpha_pred
#> [1] 2.136161

# Predicted type I error at each stages
sapply(1:n_stage, function(stage){
   1-Maxcombo.beta.n(
     setting_stoch$Sigma0[1:(n_FHweights * (stage)), 
                          1:(n_FHweights * (stage))],
     0, # mean is 0 under null
     setting_stoch$z_alpha_vec_pred[1:(n_FHweights * (stage))],
     rep(
       unlist(setting_stoch$interim_pred0), 
       each = n_FHweights)[1:(n_FHweights * (stage))],
     R,
     setting_stoch$n_FH)
   })
#> [1] 0.004038821 0.008164345 0.013367219 0.024983513
# vs. the designed errors spent at each stage
error_spend
#> [1] 0.004035615 0.008187950 0.013471082 0.025000000
# Predicted power at each stage
sapply(1:length(interim_ratio), function(stage){
   1-Maxcombo.beta.n(
     setting_stoch$Sigma1[1:(n_FHweights * (stage)), 
                          1:(n_FHweights * (stage))],
     setting_stoch$mu1[1:(n_FHweights * (stage))],
     setting_stoch$z_alpha_vec_pred[1:(n_FHweights * (stage))],
     rep(
       unlist(setting_stoch$interim_pred1), 
       each = n_FHweights)[1:(n_FHweights * (stage))],
     R,
     setting_stoch$n_FH)
   })
#> [1] 0.2927187 0.5020424 0.6794439 0.9001261
```

### Exact prediction approach

``` r
# Obtain the settings from exact prediction
setting_exact <- GSMC_design(
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
  stoch = F # disable stochastic approach
)
# Boundaries at each stage
setting_exact$z_alpha_pred
#> [1] 2.867160 2.694316 2.548900 2.316180
# Event numbers to pause/stop
setting_exact$d_fixed
#> [1] 189 220 252 314

# Boundary if only consider the last stage (without interim stages)
setting_exact$z_final_alpha_pred
#> [1] 2.135569

# Predicted type I error at each stage
sapply(1:length(interim_ratio), function(stage){
   1-Maxcombo.beta.n(
     setting_exact$Sigma0[1:(n_FHweights * (stage)), 
                          1:(n_FHweights * (stage))],
     0, # mean is 0 under null
     setting_exact$z_alpha_vec_pred[1:(n_FHweights * (stage))],
     rep(
       unlist(setting_exact$interim_pred0), 
       each = n_FHweights)[1:(n_FHweights * (stage))],
     R,
     setting_exact$n_FH)
   })
#> [1] 0.004010273 0.008090862 0.013413109 0.024922578
# vs. the designed errors spent at each stage
error_spend
#> [1] 0.004035615 0.008187950 0.013471082 0.025000000
# Predicted power at each stage
sapply(1:length(interim_ratio), function(stage){
   1-Maxcombo.beta.n(
     setting_exact$Sigma1[1:(n_FHweights * (stage)), 
                          1:(n_FHweights * (stage))],
     setting_exact$mu1[1:(n_FHweights * (stage))],
     setting_exact$z_alpha_vec_pred[1:(n_FHweights * (stage))],
     rep(
       unlist(setting_exact$interim_pred1), 
       each = n_FHweights)[1:(n_FHweights * (stage))],
     R,
     setting_exact$n_FH)
   })
#> [1] 0.2906341 0.5054211 0.6794938 0.9001279
```

## References

1.  Lakatos, E. (1988). Sample sizes based on the log-rank statistic in
    complex clinical trials. Biometrics, 229-241.

2.  Hasegawa, T. (2014). Sample size determination for the weighted
    log‐rank test with the Fleming–Harrington class of weights in
    cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.

3.  Hasegawa, T. (2016). Group sequential monitoring based on the
    weighted log‐rank test statistic with the Fleming–Harrington class
    of weights in cancer vaccine studies. Pharmaceutical statistics,
    15(5), 412-419.

4.  Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019).
    Design and monitoring of survival trials in complex scenarios.
    Statistics in medicine, 38(2), 192-209.

5.  Wang, L., Luo, X., & Zheng, C. (2021). A Simulation-free Group
    Sequential Design with Max-combo Tests in the Presence of
    Non-proportional Hazards. Journal of Pharmaceutical Statistics.

Other relevant R packages:
