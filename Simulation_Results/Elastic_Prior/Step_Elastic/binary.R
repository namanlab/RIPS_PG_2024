#---------------------------------------------------------------------------------#
# This file includes R code for implementing the Elastic Prior                    #
# in two arms superiority trial with binary endpoint based on one analysis.       #
# Two functions are included:                                                    #
# 1) T_distr() for generating the distribution of congruence measure T in the     #
#    homogeneous case                                                             #
# 2) binary_step_elastic() for generating operating characteristics of EP2,       #
#    including probability of claiming efficacy, prior effective sample size      #
#---------------------------------------------------------------------------------#


library(LaplacesDemon)
library(invgamma)


## settings
nc <- 40 # current control size
nt <- 80 # current treatment size
n0 <- 100 # historical control size
pc <- 0.4 # true mean of control

### return KS distribution under homogeneous and heterogeneous ###
# Inputs:
# x0: historical data
# n0: sample size for historical data
# nc: sample size for control arm
# R: number of simulations
T_distr <- function(x0, n0, nc, R=50000){
  set.seed(2)
  p0 <- x0/n0
  K <- numeric(R)
  for (i in 1:R) {
    y <- rbinom(1, nc, p0)
    phat <- (x0 + y)/(n0 + nc)
    obs <- cbind(c(x0, y), c(n0-x0, nc-y))
    exc <- cbind(c(n0, nc)*phat, c(n0, nc)*(1-phat))
    Ka <- max(n0, nc)^(-1/4)*sum((obs-exc)^2/exc)
    K[i] <- Ka
  }
  return(K)
}


#-----Function to obtain the operating characteristics of EP2-------------------#
# Inputs:
# T1: q1th percentile of congruence measure T in the homogeneous case
# n0: sample size for historical data
# x0: historical data
# nc: sample size for control arm
# pc: response rate for control arm
# nt: sample size for treatment arm
# pt: response rate for treatment arm
# cutoff: the probability cutoff, calibrated through simulations
# t_alpha: parameter in Beta prior of response rate for treatment arm
# t_beta: parameter in Beta prior of response rate for treatment arm
# c_alpha: parameter in Beta prior of response rate for control arm
# c_beta: parameter in Beta prior of response rate for control arm
# ntrial: number of simulated trial
#
# Outputs:
# prob.rej: probability of claiming efficacy
# EHSS: prior effective sample size
binary_step_elastic <- function(T1, n0, x0, nc, pc, nt, pt, cutoff, t_alpha, t_beta, c_alpha, c_beta, ntrial){
  rej_null <- 0
  ess <- NULL
  for (trial in 1:ntrial) {
    set.seed(trial+100)
    # generate control data
    xc <- rbinom(1, nc, pc)
    # calculate g(t)
    phat <- (x0 + xc)/(n0 + nc)
    O <- cbind(c(x0, xc), c(n0-x0, nc-xc))
    E <- cbind(c(n0, nc) * phat, c(n0, nc) * (1 - phat))
    T <- max(n0, nc)^(-1/4)*sum((O - E)^2/E)
    # obtain posterior sample for control arm
    if(T <= T1){
      gt <- 1
      c_pa1 <- (c_alpha + x0)*gt + xc
      c_pa2 <- (c_beta + n0 - x0)*gt + nc - xc
      pc_post <- rbeta(10000, c_pa1, c_pa2)
    }else{
      gt <- 0
      pc_post <- rbeta(10000, c_alpha + xc, c_beta + nc - xc)
    }
    # record prior effective sample size
    ess <- c(ess, gt*n0)
    # get posterior sample for treatment arm
    xt <- rbinom(1, nt, pt)
    pt_post <- rbeta(10000, t_alpha + xt, t_beta + nt - xt)
    # probability that treatment is superior to control
    pp <- mean(pt_post > pc_post)
    if(pp > cutoff){
      rej_null <- rej_null + 1
    }
  }
  prob.rej <- rej_null/ntrial
  EHSS <- mean(ess)
  cat("probability of claiming efficacy is", prob.rej, "\n")
  cat("prior effective sample size is", formatC(EHSS, digits = 2, format = "f"), sep = " ", "\n")
}



#########----------------------------------------------------------

######### scenario: match u0 = uc ##############
## historical data
p0 <- 0.4 # mean of historical
x0 <- 40

## obtain q1th percentile of congruence measure T
K <- T_distr(x0, n0, nc, R=50000)
T1 <- quantile(K, 0.91)

### type I error
pt <- 0.4 # true mean of treatment
timestart <- Sys.time()
binary_step_elastic(T1, n0, x0, nc, pc, nt, pt, cutoff=0.945, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

### power
pt <- 0.6 # true mean of treatment
timestart <- Sys.time()
binary_step_elastic(T1, n0, x0, nc, pc, nt, pt, cutoff=0.945, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)


#########----------------------------------------------------------

######### scenario: match u0 != uc ##############
## historical data
p0 <- 0.5 # mean of historical
x0 <- 50

## obtain q1th percentile of congruence measure T
K <- T_distr(x0, n0, nc, R=50000)
T1 <- quantile(K, 0.91)

### type I error
pt <- 0.4 # true mean of treatment
timestart <- Sys.time()
binary_step_elastic(T1, n0, x0, nc, pc, nt, pt, cutoff=0.945, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

### power
pt <- 0.6 # true mean of treatment
timestart <- Sys.time()
binary_step_elastic(T1, n0, x0, nc, pc, nt, pt, cutoff=0.945, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)


#########----------------------------------------------------------

######### scenario: match u0 != uc ##############
## historical data
p0 <- 0.6 # mean of historical
x0 <- 60

## obtain q1th percentile of congruence measure T
K <- T_distr(x0, n0, nc, R=50000)
T1 <- quantile(K, 0.91)

### type I error
pt <- 0.4 # true mean of treatment
timestart <- Sys.time()
binary_step_elastic(T1, n0, x0, nc, pc, nt, pt, cutoff=0.945, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

### power
pt <- 0.6 # true mean of treatment
timestart <- Sys.time()
binary_step_elastic(T1, n0, x0, nc, pc, nt, pt, cutoff=0.945, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)














