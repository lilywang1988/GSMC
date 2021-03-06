% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Maxcombo_size.R
\name{Maxcombo.bd}
\alias{Maxcombo.bd}
\title{Boundary calculation for GSMC}
\usage{
Maxcombo.bd(Sigma0, index, alpha_sp, n.rep = 5)
}
\arguments{
\item{Sigma0}{correlation matrix for all the test statistics.}

\item{index}{vector of non-decreasing integer starting from 1 to indicate which stage each column or row of the correlation matrix \code{Sigma0} corresponds to.}

\item{alpha_sp}{vector of cumulative errors to spend up to each stage.}

\item{n.rep}{number of repeats to take the median for output since the called likelihood generator of a multivariate normal distribution \code{\link[mvtnorm]{pmvnorm}} is not determinant. The default \code{n.rep} value is 5.}
}
\value{
\item{z_alpha}{boundary values for all the stages.}
\item{z_alpha_vec}{boundary values for all the test statistics following the \code{index}. }
}
\description{
Boundary calculation for interim analysis with max-combo tests based on correlation matrix and the alpha spending function.
}
\details{
Suppose there are 2 stages (1 interim, 1 final), and two tests for a max-combo in each stage, then we have totally 4 test statistics. Let the alpha spending function to be \code{c(alpha1,alpha)}, and the first two (\eqn{Z_{11},Z_{12}}) share one cutoff value z1, the latter two share another two (\eqn{Z_{21},Z_{22}}) share another cutoff value z2. Controlling the type I error is equivalent to ensuring that \eqn{P(Z_{11}<z_1,Z_{12}<z_1)=\alpha_1} and \eqn{P(Z_{11}<z_1,Z_{12}<z_1,Z_{21}<z_2,Z_{22}<z_2)=\alpha} are both satisfied. Note that the vector \eqn{[Z_{11},Z_{12},Z_{21},Z_{22}]^T\sim MVN(0,\Sigma_0)}. \code{Sigma0} corresponds to \eqn{\Sigma_0}, \code{index} records the ordered stages of each test statistics, whose order should be identical to the order of rows or columns in \code{Sigma0}.  Specifically, in this example, \code{index} should be \code{c(1,1,2,2)}. \code{alpha_sp} is the alpha spending function, which records how much type I error you would like to spend up to every stage.
}
\examples{
 \dontrun{
  #install.packages("gsDesign")
  library(gsDesign)
  alpha=0.025
  beta=0.1
  # If there are two stages (K=2), with on interim stage and a final stage
  # First we obtain the errors spent at each stage to be identical to the
   ones from regular interim analysis assuming that the interim stage
    happened at 60\% of events have been observed. The error spending
     function used below is O\'Brien-Fleming.
  x <- gsDesign::gsDesign(
  k = 2, 
  test.type = 1, 
  timing = 0.6, 
  sfu = "OF", 
  alpha = alpha, 
  beta = beta,
  delta = -log(0.7))
  (z <- x$upper$bound)
  x
  Sigma0_v <- rep(0.5, 6)
  Sigma0 <- matrix(1, ncol = 4, nrow = 4)
  Sigma0[upper.tri(Sigma0)] <- Sigma0_v
  Sigma0[lower.tri(Sigma0)] <- t(Sigma0)[lower.tri(t(Sigma0))]
  Sigma0
  # The error you would like to spend at the interim stage:
  alpha_interim <- pnorm(z[1],lower.tail = F) 
  
  zz <- Maxcombo.bd(
  Sigma0 = Sigma0,
  index = c(1, 1, 2, 2),
  alpha_sp = c(alpha_interim, alpha))
  
 
  # boundary value for each stage
  zz$z_alpha 
  # boundary value for each test statistic correponding to index
  zz$z_alpha_vec 
  mvtnorm::pmvnorm(
  upper = rep(zz$z_alpha[1], 2),
  corr = Sigma0[1:2,1:2]
  )[[1]]
  
  1-alpha_interim
  1-mvtnorm::pmvnorm(
  upper = zz$z_alpha_vec,
  corr = Sigma0
  )[[1]]
  
  alpha
  # What if we do not consider interim stage but with only a final stage? 
  zz1 <- Maxcombo.bd(
  Sigma0 = Sigma0[3:4,3:4],
  index = c(1,1),
  alpha_sp = c(alpha)
  )
  mvtnorm::pmvnorm(
  upper = rep(zz1$z_alpha, 2),
  corr = Sigma0[1:2, 1:2]
  )[[1]]
  
  1-alpha
  # This function will also fit 2 or any number of interims (K>=3)
  # Let there are 3 stages, Let us try controlling the error spent 
  at each stage.
  stage_p <- c(0.5,0.7,0.8,0.9)
  x <- gsDesign::gsDesign(k=5, test.type=1, timing=stage_p, sfu="OF", 
  alpha=alpha, beta=beta,delta=-log(0.7))
  (z <- x$upper$bound)
  alpha_sp<- cumsum(x$upper$prob[,1]) # the theoretical cumulative
   errors spent at each stage
# 2 tests per stage
Sigma0_v<-rep(0.5,choose(10,2))
Sigma0<-matrix(1, ncol=10,nrow=10)
Sigma0[upper.tri(Sigma0)]<- Sigma0_v
Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(t(Sigma0))]
Sigma0

zz<-Maxcombo.bd(Sigma0 = Sigma0,index=c(1,1,2,2,3,3,4,4,5,5),alpha_sp=alpha_sp)

zz$z_alpha # boundary value for each stage
zz$z_alpha_vec # boundary value for each test statistic correponding to index
# interim 1
mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1],2),corr=Sigma0[1:2,1:2])[[1]] # expected error spent at this stage
1-alpha_sp[1] #compare with the expected error spent at this stage
# above two rows are very close to each other, same for the following pairs.
# interim 2
mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:2],each=2),corr=Sigma0[1:4,1:4])[[1]]
1-alpha_sp[2]
# interim 3
mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:3],each=2),corr=Sigma0[1:6,1:6])[[1]]
1-alpha_sp[3]
# interim 4
mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:4],each=2),corr=Sigma0[1:8,1:8])[[1]]
1-alpha_sp[4]
# final stage
mvtnorm::pmvnorm(upper=rep(zz$z_alpha[1:5],each=2),corr=Sigma0[1:10,1:10])[[1]]
1-alpha_sp[5]
 }
}
\references{
Wang, L., Luo, X., & Zheng, C. (2021). A Simulation-free Group Sequential Design with Max-combo Tests in the Presence of Non-proportional Hazards. Journal of Pharmaceutical Statistics.
}
\author{
Lili Wang
}
