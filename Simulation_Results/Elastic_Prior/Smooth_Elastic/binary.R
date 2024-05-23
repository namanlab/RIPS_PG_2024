#---------------------------------------------------------------------------------#
# This file includes R code for implementing the Elastic Prior                    #
# in two arms superiority trial with binary endpoint based on one analysis.       #
# Two functions are included:                                                    #
# 1) decide_para() for determining tuning parameter a and b in logistic elastic   #
#    function                                                                     #
# 2) binary_logistic_elastic() for generating operating characteristics of EP1,   #
#    including probability of claiming efficacy, prior effective sample size      #
#---------------------------------------------------------------------------------#


library(LaplacesDemon)
library(invgamma)


## settings
nc <- 40 # current control size
nt <- 80 # current treatment size
n0 <- 100 # historical control size
pc <- 0.4 # true mean of control

#-----Function to decide tuning parameter a and b in logistic elastic function------------#
# Inputs:
# x0: historical data
# n0: sample size for historical data
# nc: sample size for control arm
# gamma: clinically highly meaningful difference
# q1: q1th percentile of congruence measure T in the homogeneous case
# q2: q2th percentile of congruence measure T in the heterogeneous case
# c: pre-specified tuning parameter controlling the shape of elastic function. Here, we fix it by 1
# small: value of elastic function in the heterogeneous case
# large: value of elastic function in the homogeneous case
# R: the number of simulations
#
# Outputs:
# a: tuning parameter in elastic function
# b: tuning parameter in elastic function
decide_para <- function(x0, n0, nc, gamma, q1, q2, c, small, large, R=50000){
  set.seed(2)
  p0 <- x0/n0
  if(p0-gamma<0){
    p <- c(p0, p0 + gamma)
  }else if(p0 + gamma>1){
    p <- c(p0, p0 - gamma)
  }else{
    p <- c(p0, p0 - gamma, p0 + gamma)
  }
  K <- matrix(NA, R, length(p))
  for (i in 1:R) {
    for (j in 1:length(p)) {
      y <- rbinom(1, nc, p[j])
      phat <- (x0 + y)/(n0 + nc)
      obs <- cbind(c(x0, y), c(n0-x0, nc-y))
      exc <- cbind(c(n0, nc)*phat, c(n0, nc)*(1-phat))
      Ka <- max(n0, nc)^(-1/4)*sum((obs-exc)^2/exc)
      K[i,j] <- Ka
    }
  }
  if(length(p)==2){
    quant1 <- quantile(K[,1], probs = q1)
    quant2 <- quantile(K[,2], probs = q2)
    K_homo <- quant1
    K_hete <- quant2
  }else{
    quant1 <- quantile(K[,1], probs = q1)
    quant2 <- quantile(K[,2], probs = q2)
    quant3 <- quantile(K[,3], probs = q2)
    K_homo <- quant1
    K_hete <- min(quant2, quant3)
  }
  b <- log((1-large)*small/((1-small)*large))/((log(K_homo))^c-(log(K_hete))^c)
  a <- log((1-large)/large)-b*(log(K_homo))^c
  return(list(a=a, b=b))
}

#-----Function to obtain the operating characteristics of EP1-------------------#
# Inputs:
# a: calibrated parameter obtained from decide_para()
# b: calibrated parameter obtained from decide_para()
# c: pre-specified tuning parameter in elastic function, fixed by c=1
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
binary_logistic_elastic <- function(a, b, c, n0, x0, nc, pc, nt, pt, cutoff, t_alpha, t_beta, c_alpha, c_beta, ntrial){
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
    gt <- 1/(1+exp(a + b*(log(T))^c))
    if(gt < 0.01) gt <- 0.01
    # record prior effective sample size
    ess <- c(ess, gt*n0)
    # obtain posterior sample for control arm
    c_pa1 <- (c_alpha + x0)*gt + xc
    c_pa2 <- (c_beta + n0 - x0)*gt + nc - xc
    pc_post <- rbeta(10000, c_pa1, c_pa2)
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

## obtain calibrated a and b in smooth elastic function
para <- decide_para(x0, n0, nc, gamma=0.24, q1=0.9, q2=0.1, c=1, small=0.01, large=0.99, R=50000)
a <- para$a 
b <- para$b

### type I error
pt <- 0.4 # true mean of treatment
timestart <- Sys.time()
binary_logistic_elastic(a, b, c=1, n0, x0, nc, pc, nt, pt, cutoff=0.942, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

### power
pt <- 0.6 # true mean of treatment
timestart <- Sys.time()
binary_logistic_elastic(a, b, c=1, n0, x0, nc, pc, nt, pt, cutoff=0.942, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)


#########----------------------------------------------------------

######### scenario: match u0 != uc ##############
## historical data
p0 <- 0.5 # mean of historical
x0 <- 50

## obtain calibrated a and b in smooth elastic function
para <- decide_para(x0, n0, nc, gamma=0.24, q1=0.9, q2=0.1, c=1, small=0.01, large=0.99, R=50000)
a <- para$a 
b <- para$b

### type I error
pt <- 0.4 # true mean of treatment
timestart <- Sys.time()
binary_logistic_elastic(a, b, c=1, n0, x0, nc, pc, nt, pt, cutoff=0.942, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

### power
pt <- 0.6 # true mean of treatment
timestart <- Sys.time()
binary_logistic_elastic(a, b, c=1, n0, x0, nc, pc, nt, pt, cutoff=0.942, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)


#########----------------------------------------------------------

######### scenario: match u0 != uc ##############
## historical data
p0 <- 0.6 # mean of historical
x0 <- 60

## obtain calibrated a and b in smooth elastic function
para <- decide_para(x0, n0, nc, gamma=0.24, q1=0.9, q2=0.1, c=1, small=0.01, large=0.99, R=50000)
a <- para$a 
b <- para$b

### type I error
pt <- 0.4 # true mean of treatment
timestart <- Sys.time()
binary_logistic_elastic(a, b, c=1, n0, x0, nc, pc, nt, pt, cutoff=0.942, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)

### power
pt <- 0.6 # true mean of treatment
timestart <- Sys.time()
binary_logistic_elastic(a, b, c=1, n0, x0, nc, pc, nt, pt, cutoff=0.942, t_alpha=0.1, t_beta=0.1, c_alpha=0.1, c_beta=0.1, ntrial=100)
timeend <- Sys.time()
print(timeend - timestart)






