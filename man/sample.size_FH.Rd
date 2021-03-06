% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hasegawa2014.R
\name{sample.size_FH}
\alias{sample.size_FH}
\title{Sample size calculation for Fleming-Harrington weighted log-rank tests
This sample size calculation method was proposed by Hasegawa (2014).
This function is to calculate the sample size for Fleming-Harrington weighted log-rank tests with piece-wise exponential distributed survival curves in described in Hasegawa(2014).}
\usage{
sample.size_FH(
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
)
}
\arguments{
\item{eps}{The change point, before which, the hazard ratio is 1, and after which, the hazard ratio is theta}

\item{p}{Treatment assignment probability.}

\item{b}{The number of subintervals per time unit.}

\item{tau}{The end of the follow-up time  in the study. Note that this is identical to \eqn{T+\tau} in the paper from Hasegawa (2014).}

\item{omega}{The minimum follow-up time for all the patients.  Note that Hasegawa(2014) assumes that the accrual is uniform between time 0 and R, and there does not exist any censoring except for the administrative censoring at the ending time \eqn{\tau}. Thus this value omega is equivalent to \code{tau-R}. Through our simulation tests, we found that this function is quite robust to violations of these assumptions: dropouts, different cenosring rates for two  arms, and changing accrual rates.}

\item{lambda}{The hazard for the control group.}

\item{lambda.trt}{The hazard for the treatment group after time eps.}

\item{rho}{The first parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.}

\item{gamma}{The second parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.}

\item{alpha}{Type I error.}

\item{beta}{Type II error.}
}
\value{
\item{n}{The needed sample size.}
\item{n_event}{The needed  event numbers for both arms together.}
\item{E.star}{The unit mean, correspoinding to \eqn{E^*} in Hasegawa(2014)}
\item{sum_D}{The cumulative D, and ceiling(n*D) is quivalent to n_vent. }
}
\description{
Sample size calculation for Fleming-Harrington weighted log-rank tests
This sample size calculation method was proposed by Hasegawa (2014).
This function is to calculate the sample size for Fleming-Harrington weighted log-rank tests with piece-wise exponential distributed survival curves in described in Hasegawa(2014).
}
\note{
This function is based on a R function from Dr. Ting Ye's paper
: Ye, T., & Yu, M. (2018). A robust approach to sample size calculation in cancer immunotherapy trials with delayed treatment effect. Biometrics, 74(4), 1292-1300.
}
\examples{
\dontrun{ 
# Example 1 from Hasegawa (2014)
p<-2/3
tau<-66
omega<-18
eps<-6
m1=21.7  #median survival time for placebo group
m2=25.8  # median survival time for treatment group
lambda<-log(2)/m1
lambda.trt<-log(2)*(m1-eps)/(m2-eps)/m1
theta=lambda.trt/lambda
alpha<-0.025
beta<-0.1
rho=0
gamma=1
b=30
sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)$n
#1974, identical to the paper's report
}
}
\references{
Ye, T., & Yu, M. (2018). A robust approach to sample size calculation in cancer immunotherapy trials with delayed treatment effect. Biometrics, 74(4), 1292-1300.
Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.
}
\author{
Lili Wang, Ting Ye
}
