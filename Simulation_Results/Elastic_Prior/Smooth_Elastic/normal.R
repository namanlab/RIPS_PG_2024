#--------------------------------------------------------------------------------#
# This file includes R code for implementing the Elastic Prior                   #
# in two arms superiority trial with normal endpoints based on one analysis.     #
# Three functions are included:                                                            #
# 1) decide_para() for determining tuning parameter a and b in logistic elastic  #
#    function                                                                    #
# 2) sample_poster() for sampling posterior of mean value for control arm        #
# 3) normal_logistic_elastic() for generating operating characteristics of EP1,  #
#    including probability of claiming efficacy, prior effective sample size     #
#--------------------------------------------------------------------------------#

library(LaplacesDemon)
library(invgamma)

## settings
nc <- 25 # current control size
nt <- 50 # current treatment size
n0 <- 50 # historical control size
sigc <- 1 # control sd
sigt <- 1 # treatment sd
sig0 <- 1 # historical sd
uc <- 1 # true mean of control


#################
################# EP1 with smooth elastic function #######################
#################

#-----Function to decide tuning parameter a and b in logistic elastic function------------#
# Inputs:
# c: pre-specified tuning parameter controlling the shape of elastic function. Here, we fix it by 1
# x0: historical data
# n0: sample size for historical data
# nc: sample size for control arm
# gamma: clinically highly meaningful difference,  -- 
# q1: q1th percentile of congruence measure T in the homogeneous case, 0.95 --
# q2: q2th percentile of congruence measure T in the heterogeneous case, 0.02 --
# small: value of elastic function in the heterogeneous case, 0.01 --
# large: value of elastic function in the homogeneous case, 0.99 --
# R: the number of simulations
# 
# Outputs:
# a: tuning parameter in elastic function 
# b: tuning parameter in elastic function
decide_para <- function(c, x0, n0, nc, gamma, q1, q2, small, large, R){
  set.seed(1)
  u0 <- mean(x0)
  sig0 <- sd(x0)
  mc <- c(u0, u0 + gamma,  u0 - gamma)
  t <- matrix(NA, R, length(mc))
  for (i in 1:R) {
    for (j in 1:length(mc)) {
      xc <- rnorm(nc, mc[j], sig0)
      sp <- ((n0-1)*sig0^2 + (nc-1)*var(xc))/(n0 + nc - 2) # pooled variance
      t[i,j] <- max(n0, nc)^(-1/4)*abs(u0-mean(xc))/(sqrt(sp/n0 + sp/nc))
    }
  }
  quant1 <- quantile(t[,1], q1)
  quant2 <- quantile(t[,2], q2)
  quant3 <- quantile(t[,3], q2)
  KS_homo <- quant1
  KS_hete <- min(quant2, quant3)
  b <- log((1-large)*small/((1-small)*large))/((log(KS_homo))^c-(log(KS_hete))^c)
  a <- log((1-large)/large)-b*(log(KS_homo))^c
  return(list(a=a, b=b))
}

#-----Function to sample posterior of mean value for control arm------------#
# Inputs:
# x0: historical data
# n0: sample size for historical data
# xc: control data
# nc: sample size for control arm
# gt: value of elastic function at current congruence measure t between historical and control data
# sim: number of simulated trial
# nburn: number of simulated trial in burn-in process
# 
# Outputs:
# muc_post: posterior of mean value for control arm
sample_poster <- function(x0, n0, xc, nc, gt, sim=20000, nburn=10000){
  muc_post <- numeric(sim-nburn)
  muc_ini <- mean(xc)
  muc <- muc_ini
  for (s in 1:sim) {
    # sample control variance from posterior
    alpha <- (1 + nc)/2
    beta <- nc*(var(xc) + (mean(xc) - muc)^2)/2
    sig <- rinvgamma(1, alpha, beta)
    # sample mean from posterior
    D <- var(x0)/(n0*gt)
    mu <- (nc*mean(xc)*D + sig*mean(x0))/(nc*D + sig)
    var <- sig*D/(D*nc + sig)
    muc <- rnorm(1, mu, sqrt(var)) 
    if(s > nburn){
      muc_post[s-nburn] <- muc
    }
  }
  return(list(muc_post=muc_post))
}

#-----Function to obtain the operating characteristics of EP1-------------------#
# Inputs:
# a: calibrated parameter obtained from decide_para()
# b: calibrated parameter obtained from decide_para()
# c: pre-specified tuning parameter controlling the shape of elastic function. Here, we fix it by 1
# n0: sample size for historical data
# x0: historical data
# nc: sample size for control arm
# uc: true mean value for control arm
# sigc: true standard deviation for control arm
# nt: sample size for treatment arm
# ut: true mean value for treatment arm
# sigt: true standard deviation for treatment arm
# cutoff: the efficacy evaluation probability cutoff, calibrated through simulations
# ntrial: number of simulated trial
# 
# Outputs:
# prob.rej: probability of claiming efficacy
# EHSS: prior effective sample size
normal_logistic_elastic <- function(a, b, c, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff, ntrial){
  rej_null <- 0
  ess <- NULL
  for(trial in 1:ntrial){
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rnorm(nc, uc, sigc)
    xt <- rnorm(nt, ut, sigt)
    # calculate statistic and g(t)
    sp <- ((n0-1)*var(x0) + (nc-1)*var(xc))/(n0 + nc - 2) # pooled variance
    T <- max(n0, nc)^(-1/4)*abs(mean(x0)-mean(xc))/(sqrt(sp/n0 + sp/nc))
    gt <- 1/(1 + exp(a + b*(log(T))^c))
    if(gt == 0) gt <- 0.00001 # numerical stability
    # record prior effective sample size
    ess <- c(ess, gt*n0)
    # posterior for control arm
    temp <- sample_poster(x0, n0, xc, nc, gt)
    muc <- temp$muc_post
    # posterior for treatment arm
    t.par1 <- mean(xt)
    t.par2 <- var(xt)/nt
    mut <- rst(10000, t.par1, sqrt(t.par2), nt - 1)
    # probability that treatment is superior to control
    pp <- mean(mut > muc)
    if(pp > cutoff){
      rej_null <- rej_null + 1
    }
  }
  EHSS <- mean(ess)
  prob.rej <- rej_null/ntrial
  cat("probability of claiming efficacy is", prob.rej, "\n")
  cat("effective historical sample size is", formatC(EHSS, digits = 2, format = "f"), sep = " ", "\n")
}


#########----------------------------------------------------------


######### scenario: match u0 = uc ##############
## historical data
set.seed(8172)
u0 <- 1 # mean of historical
x0 <- rnorm(n0, u0, sig0)

## obtain calibrated a and b in smooth elastic function
para <- decide_para(c=1, x0, n0, nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 50000)
a <- para$a 
b <- para$b 

## type I error
ut <- 1 # true mean of treatment
timestart <- Sys.time()
normal_logistic_elastic(a, b, c=1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

## power
ut <- 1.5 # true mean of treatment
timestart <- Sys.time()
normal_logistic_elastic(a, b, c=1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

#########----------------------------------------------------------


######### scenario: match u0 != uc ##############
## historical data
set.seed(8172)
u0 <- 1.2 # mean of historical
x0 <- rnorm(n0, u0, sig0)

## obtain calibrated a and b in smooth elastic function
para <- decide_para(c=1, x0, n0, nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 50000)
a <- para$a 
b <- para$b 

## type I error
ut <- 1 # true mean of treatment
timestart <- Sys.time()
normal_logistic_elastic(a, b, c=1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

## power
ut <- 1.5 # true mean of treatment
timestart <- Sys.time()
normal_logistic_elastic(a, b, c=1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

#########----------------------------------------------------------

######### scenario: match u0 != uc ##############
## historical data
set.seed(8172)
u0 <- 1.5 # mean of historical
x0 <- rnorm(n0, u0, sig0)

## obtain calibrated a and b in smooth elastic function
para <- decide_para(c=1, x0, n0, nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 50000)
a <- para$a 
b <- para$b 

## type I error
ut <- 1 # true mean of treatment
timestart <- Sys.time()
normal_logistic_elastic(a, b, c=1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

## power
ut <- 1.5 # true mean of treatment
timestart <- Sys.time()
normal_logistic_elastic(a, b, c=1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)






