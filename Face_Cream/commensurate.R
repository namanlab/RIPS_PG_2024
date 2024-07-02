library(MCMCpack)
library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)
library(mvtnorm)

run_simulation <- function(sig, tau, uh, nh, uc, Xc, Zc, N, R, cutoff){
  
  set.seed(42)
  xh <- rnorm(nh, uh, sqrt(sig^2 + tau^2))
  power <- 0
  
  for(trial in 1:R){
    
    set.seed(100 + trial)
    nc <- ncol(Zc)
    overall_noise <- rnorm(2*nc, 0, sig^2)
    indiv_noise <- rnorm(nc, 0, tau^2)
    xc <- Xc %*% uc + Zc %*% indiv_noise + overall_noise
    pc <- ncol(Xc)
    
    # posterior for control arm -- TO CHANGE DEPENDING ON METHOD
    res_inter <- get_control_prost(Xc, Zc, xc, xh, nc, nh, N, pc)
    muc <- res_inter$prost_samples
    print(apply(muc, 2, mean))
    res <- 0
    for (i in 1:(pc - 1)){
      curv = mean(muc[,i] <= muc[,pc])
      res <- res + ifelse(curv >= cutoff, 1, 0)
    }
    power <- power + ifelse(res == pc - 1, 1, 0)
    print(power)
  }
  cat("Power:", power/R)
  
}


get_control_prost <- function(Xc, Zc, xc, xh, nc, nh, N, pc) {
  
  # Define the prior distribution for a0 depending on tau
  log_prior_a0 <- function(a0, tau) {
    g_tau <- function(tau) {
      return(max(tau, 1))
    }
    return((g_tau(tau) - 1) * log(a0))
  }
  
  # Define the prior distribution for tau
  log_prior_tau <- function(tau) {
    log(dcauchy(log(tau), 0, 30)) - log(tau)
  }
  
  # Logit and inverse logit functions
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  # Define the posterior distribution for control arm using commensurate power prior
  log_posterior <- function(theta) {
    mu <- theta[1:pc]
    sigma2 <- exp(theta[pc + 1])
    tau2 <- exp(theta[pc + 2])
    a0 <- inv_logit(theta[pc + 3])
    theta0 <- theta[pc + 4]
    tau <- exp(theta[pc + 5])
    
    # Likelihood
    log_lik_current <- dmvnorm(t(xc), mean = Xc %*% mu, sigma = tau2*(Zc %*% t(Zc)) + sigma2*diag(2*nc), log = T)
    log_lik_hist <- sum(dnorm(xh, theta0, sqrt(sigma2 + tau2), log = T))
    
    # Prior distributions
    log_prior_sigma <- -sigma2
    log_prior_tau <- -tau2
    log_prior_a0_val <- log_prior_a0(a0, tau)
    log_prior_comm <- log_prior_tau(tau)
    log_prior_mu <- -0.5 * tau * (mu[pc] - theta0)^2
    
    integ_log <- -a0/(2*sigma2)*sum(xh^2) + nh*a0*mean(xh)^2/(2*sigma2) - log(sqrt(sigma2)/sqrt(nh*a0))
    jacob_log <- log(sigma2) + log(tau2) + log(a0*(1 - a0)) + log(tau)
    
    # Posterior
    log_post <- log_lik_current + (a0 * log_lik_hist) + log_prior_sigma + 
      log_prior_a0_val + log_prior_tau + log_prior_comm + log_prior_mu + 
      jacob_log #- integ_log
    
    if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
      return(-Inf)  # Return -Inf for invalid log-posterior values
    }
    return(log_post)
  }
  
  init_values <- c(rep(0, pc), log(1), log(1), logit(0.8), mean(xh), log(1))  # Initial values for MCMC sampling
  # Perform MCMC sampling for posterior
  samples <<- MCMCmetrop1R(log_posterior, theta.init = init_values, tune = 1,
                           mcmc = N, burnin = 1000, thin = 1, force.samp = T,
                           optim.method = "Nelder-Mead")
  
  # Extract posterior samples for mu, sigma, a0, theta0, and tau
  mu_samples <- samples[, 1:pc]
  sigma2_samples <- exp(samples[, pc + 1])
  tau2_samples <- inv_logit(samples[, pc + 2])
  list(prost_samples = mu_samples)
}



#########--------------------------------------------------------------#########
#################################### MODEL #####################################
#########--------------------------------------------------------------#########

gen_Z <- function(n){
  bm <- diag(n)
  res <- matrix(rep(0, 2*n^2), nrow = 2*n)
  for (i in 1:(2*n)){
    print(bm[floor(i/2),])
    res[i,] <- bm[ceiling(i/2),]
  }
  return(res)
}

gen_X <- function(n, p){
  res <- diag(p)
  for (i in 1:ceiling(2*n/p)){
    res <- rbind(res, diag(p))
  }
  res <- res[1:(2*n), ]
  return(res)
}

## settings
nc <- 45
pc <- 3
uc <- c(2.41, 2.45, 2.59)
Zc <- gen_Z(nc)
Xc <- gen_X(nc, pc)
sig = 0.2
tau = 0.2
uh = 2.67
nh = 50
run_simulation(sig, tau, uh, nh, uc, Xc, Zc, 10000, 100, 0.95)

