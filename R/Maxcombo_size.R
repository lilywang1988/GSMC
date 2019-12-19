# Author: Lili Wang
# Date: 20190823
# Modified: 20190917
n.rep=5
#' Boundary calculation for GS-MC
#'
#' Boundary calculation for interim analysis with max-combo tests based on correlation matrix and the alpha spending function.
#'
#' For example, when there are 2 stages (1 interm, 1 final), and two tests for a max-combo in each stage, then we have 4 test statistics. Let the alpha spending function to be \code{c(alpha1,alpha)}, and the first two (\eqn{Z_{11},Z_{12}}) share one cutoff value z1, the latter two share another two (\eqn{Z_{21},Z_{22}}) share another cutoff value z2. The goals is to ensure that \eqn{P(Z_{11}<z_1,Z_{12}<z_1})=\alpha_1} and \eqn{P(Z_{11}<z_1,Z_{12}<z_1},Z_{21}<z_2,Z_{22}<z_2)=\alpha}.Note that the vector \eqn{[Z_{11},Z_{12},Z_{21},Z_{22}]^T\sim MVN(0,\Sigma_0)}. \code{Sigma0} corresponds to \equn{\Sigma_0}, \code{index} records the ordered stages of each test statistics, whose order should be indentical to the order of rows or coloumns in \code{Sigma0}, in this example \code{index} should be \code{c(1,1,2,2)}.\code{alpha_sp} is the alpha spending function, which records how much type I error you would like to spend up to each stage.
#'
#' @param Sigma0 Correlation matrix for all the test statistics.
#' @param index Vector of non-decreasing integer starting from 1 to indicate which stage each column or row of the correlation matrix \code{Sigma0} corresponds to.
#' @param alpha_sp Vector of errors to spend up to each stage.
#' @return
#' \item{z_alpha}{Boundary values for all the stages.}
#' \item{z_alpha_vec}{Boundary values for all the test statistics correponding to index. }
#' @author Lili Wang
#' @examples
#'   #install.packages("gsDesign")
#'   library(gsDesign)
#'   alpha=0.025
#'   beta=0.1
#'   # If there are two stages (K=2), with on interim stage and a final stage
#'   # First we obtain the errors spent at each stage to be identical to the ones from regular interim analysis assuming that the interim stage happened at 60% of events have been observed. The error spending function used below is O'Brien-Fleming.
#'   x <- gsDesign(k=2, test.type=1, timing=0.6, sfu="OF", alpha=alpha, beta=beta,delta=-log(0.7))
#'   (z <- x$upper$bound)
#'   x
#'   Sigma0_v<-rep(0.5,6)
#'   Sigma0<-matrix(1, ncol=4,nrow=4)
#'   Sigma0[upper.tri(Sigma0)]<- Sigma0_v
#'   Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(t(Sigma0))]
#'   Sigma0
#'   alpha_interim<-pnorm(z[1],lower.tail = F) # The error you would like to spend at the interim stage
#'   zz<-Maxcombo.bd(Sigma0 = Sigma0,index=c(1,1,2,2),alpha_sp=c(alpha_interim,alpha))
#'
#'   zz$z_alpha # boundary value for each stage
#'   zz$z_alpha_vec # boundary value for each test statistic correponding to index
#'   pmvnorm(upper=rep(zz$z_alpha[1],2),corr=Sigma0[1:2,1:2])[[1]]
#'   1-alpha_interim
#'   1-pmvnorm(upper =zz$z_alpha_vec,corr=Sigma0)[[1]]
#'   alpha
#'   # What if we do not consider interim stage but with only a final stage? (K=1)
#'   zz1<-Maxcombo.bd(Sigma0 = Sigma0[3:4,3:4],index=c(1,1),alpha_sp=c(alpha))
#'   pmvnorm(upper=rep(zz1$z_alpha,2),corr=Sigma0[1:2,1:2])[[1]]
#'   1-alpha
#'   # This function will also fit 2 or any number of interims (K>=3)
#'   # Let there are 3 stages, Let us try controlling the error spent at each stage.
#'   stage_p<-c(0.5,0.7,0.8,0.9)
#'   x <- gsDesign(k=5, test.type=1, timing=stage_p, sfu="OF", alpha=alpha, beta=beta,delta=-log(0.7))
#'   (z <- x$upper$bound)
#'   alpha_sp<- cumsum(x$upper$prob[,1]) # the theoretical cumulative errors spent at each stage
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
#' pmvnorm(upper=rep(zz$z_alpha[1],2),corr=Sigma0[1:2,1:2])[[1]] # expected error spent at this stage
#' 1-alpha_sp[1] #compare with the expected error spent at this stage
#' # above two rows are very close to each other, same for the following pairs.
#' # interim 2
#' pmvnorm(upper=rep(zz$z_alpha[1:2],each=2),corr=Sigma0[1:4,1:4])[[1]]
#' 1-alpha_sp[2]
#' # interim 3
#' pmvnorm(upper=rep(zz$z_alpha[1:3],each=2),corr=Sigma0[1:6,1:6])[[1]]
#' 1-alpha_sp[3]
#' # interim 4
#' pmvnorm(upper=rep(zz$z_alpha[1:4],each=2),corr=Sigma0[1:8,1:8])[[1]]
#' 1-alpha_sp[4]
#' # final stage
#' pmvnorm(upper=rep(zz$z_alpha[1:5],each=2),corr=Sigma0[1:10,1:10])[[1]]
#' 1-alpha_sp[5]


#library(mvtnorm)  # 1.0-10 version (dependent)
#old name is Maxcombo_bound
Maxcombo.bd<-function(Sigma0,index,alpha_sp){
  K=length(alpha_sp)
  if(length(unique(index))!=K){stop("Error: The index should be consistent with the length of alpha_sp!")}
  if(any(diff(alpha_sp)<=0)) stop("Error: alpha_sp should be monotone (strictly) increasing!")
  if(any(diff(index)<0)) stop("Error: index should be monotone increasing!")
  z_alpha<-NULL
  z_alpha_vec<-NULL # a repeat vector to record the test statistic for each intex
  for(k in 1:K){
    use_index=which(index<=k)
    Sigma0_use<-Sigma0[use_index,use_index]
    if(k==1){
      z_alpha[1]=median(replicate(n.rep,qmvnorm(1-alpha_sp[k],tail="lower.tail",corr=Sigma0_use)$quantile))
    }else{
      z2_search<-function(z){
        median(replicate(n.rep,pmvnorm(upper=c(z_alpha_vec,rep(z,sum(index==k))),corr=Sigma0_use)-(1-alpha_sp[k])))
      }
      z_alpha[k]=uniroot(f=z2_search,interval=c(1,z_alpha[k-1]))$root
    }
    z_alpha_vec=c(z_alpha_vec,rep(z_alpha[k],sum(index==k)))
  }
  return(list(z_alpha=z_alpha,z_alpha_vec=z_alpha_vec))
}

#' Sample size calculation
#'
#' Sample size calculation to control the type II error or the power of an interim analysis with Maxcombo tests.
#'
#' Assume that there are 2 stages (1 interm, 1 final), and two tests for a max-combo in each stage, then we have 4 test statistics, and the two cutoff values for the two stages have been determined by \code{Maxcombo.bd} in advance using their correlation matrix and the error spending function \eqn{\alpha_1, \alpha}. The goal of this function is to control the sample size n (number of patients for both arms) or d (observed events) to achieve the ideal type II error \eqn{\beta} or the power (\eqn{1-\beta}), i.e. \eqn{\P(Z_{11}<z_1,Z_{12}<z_1,Z_{21}<z_2,Z_{22}<z_2)=\beta}.
#'
#'
#' @param Sigma1 The correlation matrix under the alternative hypothesis.
#' @param mu1 The unit mu under the alternative hypothesis (the mean of the expectation of each subject scaled weighted log-rank test statistic, which can be approximated using the fomula for \equn{E^*} in Hasegawa 2014 paper. ).
#' @param z_alpha_vec Same as the one exported from Maxcombo.bd, which is the boundaries for ordered test statistics, its order should be consistent to the rows and columns in \code{Sigma1}.
#' @param beta Type II error.
#' @param interim_vec The vector of the interims in each stages, not that it should be a repeat vector with same iterim values for all the test statitics at same stages.
#' @param R End of the enrollment time, which is identical to \code{R} defined in other functions like \code{\link{I.1}}.
#' @param n_range The range ot the expected patient numbers.
#' @param sum_D Same as the exported value from \code{\link{sample.size_FH}}, the summed \eqn{D^*} in Hasegawa (2014).
#' @return
#' \item{n}{The number of patients needed for the trial to achieve the predefined power.}
#' \item{d}{The number of events needed for the trial to achieve the predefined power.}
#' \item{sum_D}{The input \code{sum_D} value. }
#' @author Lili Wang
#' @references Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.
#' @seealso \code{\link{Maxcombo.beta.n}}
#' @examples
#' #install.packages("mvtnorm")
#' library(mvtnorm)
#' #install.packages("gsDesign")
#' library(gsDesign)
#' alpha=0.025
#' beta=0.1
#' # If there are two stages (K=2), with on interim stage and a final stage
#' # First we obtain the errors spent at each stage to be identical to the ones from regular interim analysis assuming that the interim stage happened at 60% of events have been observed. The error spending function used below is O'Brien-Fleming.
#' x <- gsDesign(k=2, test.type=1, timing=0.6, sfu="OF", alpha=alpha, beta=beta,delta=-log(0.7))
#' (z <- x$upper$bound)
#' x
#' Sigma0_v<-rep(0.5,6)
#' Sigma0<-matrix(1, ncol=4,nrow=4)
#' Sigma0[upper.tri(Sigma0)]<- Sigma0_v
#' Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(t(Sigma0))]
#' Sigma0
#' alpha_interim<-pnorm(z[1],lower.tail = F) # The error you would like to spend at the interim stage
#' zz<-Maxcombo.bd(Sigma0 = Sigma0,index=c(1,1,2,2),alpha_sp=c(alpha_interim,alpha))
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
#' Maxcombo.sz(Sigma1=Sigma1,mu1=mu1,z_alpha_vec=zz$z_alpha_vec,beta=0.1,interim_vec=c(10,10,18,18),R=14,n_range=c(100,1000),sum_D=0.6)
#' # need 232 patients, 140 deaths
#'
#'
#'
Maxcombo.sz<-function(Sigma1,mu1,z_alpha_vec,beta,interim_vec,R,n_range,sum_D){
search_n<-function(n){
  median(replicate(n.rep,pmvnorm(cor=Sigma1,upper=z_alpha_vec-sqrt(n*pmin(interim_vec/R,1))*mu1)-(beta)))
}
n=ceiling(uniroot(f=search_n,interval = n_range)$root)
d=ceiling(n*sum_D)
return(list(n=n,d=d,sum_D=sum_D))
}

#' The type II errors/Powers for a range of sample sizes
#'
#' To obtain a spectrum of powers or type II errors for a range of sample sizes n or d
#'
#' @inheritParams Maxcombo.sz
#' @param n_seq The sequence of number of patients.
#' @param d_seq The sequence of number of expected events.
#' @author Lili Wang
#' @seealso \code{\link{Maxcombo.sz}}
#' @examples
#' #install.packages("mvtnorm")
#' #library(mvtnorm)
#' #install.packages("gsDesign")
#' #library(gsDesign)
#' alpha=0.025
#' beta=0.1
#' # If there are two stages (K=2), with on interim stage and a final stage
#' # First we obtain the errors spent at each stage to be identical to the ones from regular interim analysis assuming that the interim stage happened at 60% of events have been observed. The error spending function used below is O'Brien-Fleming.
#' x <- gsDesign(k=2, test.type=1, timing=0.6, sfu="OF", alpha=alpha, beta=beta,delta=-log(0.7))
#' (z <- x$upper$bound)
#' x
#' Sigma0_v<-rep(0.5,6)
#' Sigma0<-matrix(1, ncol=4,nrow=4)
#' Sigma0[upper.tri(Sigma0)]<- Sigma0_v
#' Sigma0[lower.tri(Sigma0)]<- t(Sigma0)[lower.tri(t(Sigma0))]
#' Sigma0
#' alpha_interim<-pnorm(z[1],lower.tail = F) # The error you would like to spend at the interim stage
#' zz<-Maxcombo.bd(Sigma0 = Sigma0,index=c(1,1,2,2),alpha_sp=c(alpha_interim,alpha))
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
#' Maxcombo.sz(Sigma1=Sigma1,mu1=mu1,z_alpha_vec=zz$z_alpha_vec,beta=0.1,interim_vec=c(10,10,18,18),R=14,n_range=c(100,1000),sum_D=0.6)
#' # need 232 patients, 140 deaths
#'
#' #Obatain the spectrum of powers or type II errors in the input range
#' power_n<-1-Maxcombo.beta.n(Sigma1=Sigma1,mu1=mu1,z_alpha_vec=zz$z_alpha_vec,interim_vec=c(10,10,18,18),R=14,n_seq=seq(100,1000,50))
#' plot(x=seq(100,1000,50),y=power_n,type="l",col=1,lwd=2,main=expression(paste(1-beta," vs n")),ylab = expression(paste(1-beta)),xlab="n"     )
#' power_d<-1-Maxcombo.beta.n(Sigma1=Sigma1,mu1=mu1,z_alpha_vec=zz$z_alpha_vec,interim_vec=c(10,10,18,18),R=14,n_seq=seq(60,600,30))
#' plot(x=seq(60,600,30),y=power_d,type="l",col=1,lwd=2,main=expression(paste(1-beta," vs d")),ylab = expression(paste(1-beta)),xlab="d"     )
#' 
#' 
#'
#'
Maxcombo.beta.n <- function(Sigma1,mu1,z_alpha_vec,interim_vec,R,n_seq){
  sapply(n_seq,function(n){median(replicate(n.rep,pmvnorm(cor=Sigma1,upper=z_alpha_vec-sqrt(n*pmin(interim_vec/R,1))*mu1)[[1]]))})

}
#' @rdname Maxcombo.beta.n
Maxcombo.beta.d <- function(Sigma1,mu1,z_alpha_vec,interim_vec,R,d_seq,sum_D){
  sapply(d_seq,function(d){median(replicate(n.rep,pmvnorm(cor=Sigma1,upper=z_alpha_vec-sqrt(d/sum_D*pmin(interim_vec/R,1))*mu1)[[1]]))})
}




#' A stochastic-process way of prediction
#'
#' A stochastic-process way of prediction of the expected event counts, mean difference, and the information(variance) or the covariance
#'
#' @param eps Delayed treatment effect time.
#' @param p Probability of treatment assignment.
#' @param b The number of subintervals at each time point.
#' @param omega The minimum follow-up time for all the patients.  Note that Hasegawa(2014) assumes that the accrual is uniform between time 0 and R, and there does not exist any censoring except for the administrative censoring at the ending time \eqn{\tau}. Thus this value omega is equivalent to \code{tau-R}. Through our simulation tests, we found that this function is quite robust to violations of these assumptions: dropouts, different cenosring rates for two  arms, and changing accrual rates.
#' @param lambda The hazard for the control group.
#' @param theta The hazard ratio after the delayed time \code{eps} for the treatment arm.
#' @param rho,rho1,rho2 The first parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}.
#' @param gamma,gamma1,gamma2 The second parameter for Fleming Harrington weighted log-rank test:\eqn{W(t)=S^\rho(t^-)(1-S(t^-))^\gamma}}.
#' @param R The accrual period. 
#' 
#' @return 
#' \item{sum_D}{The mean expected event ratio, multiplied by \code{n}, the sample size, it is equal to the stochastically predicted number of events. }
#' \item{inf or covariance}{The information/variance or covariance (averaged for each subject) , should multiplied by \code{n}, which is the sample size to obtain the stochastically predicted information. }
#' \item{E.star}{The unit mean, correspoinding to \eqn{E^*} in Hasegawa(2014)}
#' \item{trt_vs_ctrl_N}{The ratio of the samples sizes between the two arms, treatment vs control, corresonding to the time vector \code{t_vec}.}
#' \item{t_vec}{The time sequence corresponding to \code{trt_vs_ctrl_N}.}
#' @author Lili Wang
#' @references Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.
#'
stoch_pred<-function(eps,p,b,tau,omega,lambda,theta,rho,gamma,R){
  n_sub<-floor(b*tau)
  t<-c(0,seq(1,n_sub)/b)
  h_1<-rep(lambda,(n_sub+1)) #control
  h_2<-c(rep(lambda,min(eps*b,n_sub+1)),rep(lambda*theta,max(n_sub-eps*b+1,0))) #treatment
  N_1<-rep((1-p),(n_sub+1))
  N_2<-rep(p,(n_sub+1))
  for(i in 1:(n_sub-1)){
    N_1[i+1]<-N_1[i]*(1-h_1[i]/b-(t[i]>omega)/b/(tau-t[i]))
    N_2[i+1]<-N_2[i]*(1-h_2[i]/b-(t[i]>omega)/b/(tau-t[i]))
  }
  N_1[n_sub+1]<-N_2[n_sub+1]<-0
  
  f_S_1<-function(x) exp(-lambda*x)
  f_S_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-(lambda*eps+lambda*theta*(x-eps)))
  #f_S_2_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-eps*lambda.trt*(1/theta-1))*exp(-lambda.trt*x)
  S_1<-f_S_1(t)
  S_2<-f_S_2(t)
  S<-(1-p)*S_1+p*S_2
  D<-(h_1*N_1+h_2*N_2)/b*min(tau/R,1) #predicted total events at each time
  theta_seq<-h_2/h_1
  phi<-N_2/N_1
  r<-S^rho*(1-S)^gamma
  num_vec<-D*r*(phi*theta_seq/(1+phi*theta_seq)-phi/(1+phi))
  den_vec<-D*r^2*phi/(1+phi)^2
  E.star_num<-sum(num_vec[1:n_sub])
  E.star_den<-sqrt(sum(den_vec[1:n_sub]))
  E.star<-E.star_num/E.star_den
  #time_vec=seq(1,n_sub)/b)
  return(list(sum_D=sum(D[1:n_sub]),inf=sum(den_vec[1:n_sub]), E.star=E.star, trt_vs_ctrl_N=phi[1:n_sub], t_vec=seq(1,n_sub)/b))
}
#' 
#' 
#' @rdname stoch_pred
#' 
stoch_pred.cov<-function(eps,p,b,tau,omega,lambda,theta,rho1,gamma1,rho2,gamma2,R){
  n_sub<-floor(b*tau)
  t<-c(0,seq(1,n_sub)/b)
  h_1<-rep(lambda,(n_sub+1)) #control
  h_2<-c(rep(lambda,min(eps*b,n_sub+1)),rep(lambda*theta,max(n_sub-eps*b+1,0))) #treatment
  N_1<-rep((1-p),(n_sub+1))
  N_2<-rep(p,(n_sub+1))
  for(i in 1:(n_sub-1)){
    N_1[i+1]<-N_1[i]*(1-h_1[i]/b-(t[i]>omega)/b/(tau-t[i]))
    N_2[i+1]<-N_2[i]*(1-h_2[i]/b-(t[i]>omega)/b/(tau-t[i]))
  }
  N_1[n_sub+1]<-N_2[n_sub+1]<-0
  
  f_S_1<-function(x) exp(-lambda*x)
  f_S_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-(lambda*eps+lambda*theta*(x-eps)))
  #f_S_2_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-eps*lambda.trt*(1/theta-1))*exp(-lambda.trt*x)
  S_1<-f_S_1(t)
  S_2<-f_S_2(t)
  S<-(1-p)*S_1+p*S_2
  D<-(h_1*N_1+h_2*N_2)/b*min(tau/R,1) #predicted total events at each time
  theta_seq<-h_2/h_1
  phi<-N_2/N_1
  r1<-S^rho1*(1-S)^gamma1
  r2<-S^rho2*(1-S)^gamma2
  den_vec<-D*r1*r2*phi/(1+phi)^2
  #time_vec=seq(1,n_sub)/b)
  return(list(sum_D=sum(D[1:n_sub]),cov=sum(den_vec[1:n_sub]), trt_vs_ctrl_N=phi[1:n_sub], t_vec=seq(1,n_sub)/b))
}

