#--------------------------------------------------------------------------------#
# This file includes R code for implementing the Elastic Prior                   #
# in two arms superiority trial with normal endpoints based on one analysis.     #
# Three functions are included:                                                   #
# 1) T_distr() for generating the distribution of congruence measure T in the    #
#    homogeneous case                                                           #
# 2) sample_poster() for sampling posterior of mean value for control arm       #
# 3) normal_step_elastic() for generating operating characteristics of EP2,      #
#    including probability of claiming efficacy, prior effective sample size     #
#--------------------------------------------------------------------------------#

library(LaplacesDemon)
library(invgamma)



## settings
nc <- 25 
nt <- 50
n0 <- 50
sigc <- 1
sigt <- 1
sig0 <- 1
uc <- 1


#################
################# EP2 with step elastic function #######################
#################


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


### return distribution of congruence measure T in the homogeneous case ###
# Inputs:
# n0: sample size for historical data
# x0: historical data
# nc: sample size for control arm
# R: number of simulations
T_distr <- function(n0, x0, nc, R=50000){
  set.seed(1)
  u0 <- mean(x0)
  sig0 <- sd(x0)
  t <- numeric(R)
  for (i in 1:R) {
    xc <- rnorm(nc, u0, sig0)
    sp <- ((n0-1)*var(x0) + (nc-1)*var(xc))/(n0 + nc - 2) # pooled variance
    t[i] <- max(n0, nc)^(-1/4)*abs(u0-mean(xc))/(sqrt(sp/n0 + sp/nc))
  }
  return(t)
}


#-----Function to obtain the operating characteristics of EP2-------------------#
# Inputs:
# T1: q1th percentile of congruence measure T in the homogeneous case
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
normal_step_elastic <- function(T1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff, ntrial){
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
    if(T <= T1){
      gt <- 1
      temp <- sample_poster(x0, n0, xc, nc, gt)
      muc <- temp$muc_post
    }else{
      gt <- 0
      muc <- rst(10000, mean(xc), sqrt(var(xc)/nc), nc - 1)
    }
    # record prior effective sample size
    ess <- c(ess, gt*n0)
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




######### match u0 == uc ##############
## historical data
set.seed(8172)
u0 <- 1
x0 <- rnorm(n0, u0, sig0)

## obtain q1th percentile of congruence measure T
t <- T_distr(n0, x0, nc, R=50000)
T1 <- quantile(t, 0.987)

## type I error
ut <- 1
timestart <- Sys.time()
normal_step_elastic(T1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

## power
ut <- 1.5
timestart <- Sys.time()
normal_step_elastic(T1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

######### match u0 != uc ##############
## historical data
set.seed(8172)
u0 <- 1.2
x0 <- rnorm(n0, u0, sig0)

## obtain q1th percentile of congruence measure T
t <- T_distr(n0, x0, nc, R=50000)
T1 <- quantile(t, 0.987)

## type I error
ut <- 1
timestart <- Sys.time()
normal_step_elastic(T1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

## power
ut <- 1.5
timestart <- Sys.time()
normal_step_elastic(T1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

######### match u0 != uc ##############
## historical data
set.seed(8172)
u0 <- 1.5
x0 <- rnorm(n0, u0, sig0)

## obtain q1th percentile of congruence measure T
t <- T_distr(n0, x0, nc, R=50000)
T1 <- quantile(t, 0.987)

## type I error
ut <- 1
timestart <- Sys.time()
normal_step_elastic(T1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

## power
ut <- 1.5
timestart <- Sys.time()
normal_step_elastic(T1, n0, x0, nc, uc, sigc, nt, ut, sigt, cutoff=0.918, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)


