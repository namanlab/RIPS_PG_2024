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
  ess <- NULL # effective sample size
  timestart <- Sys.time() # computation time
  quantile_interval_count <- 0 # number of times muc is in quantile interval
  width_quantile_interval <- NULL # width of credible interval
  point_est <- NULL # point estimator based on control prior
  
  # distrn:
  distr_bf <- NULL
  distr_at <- NULL
  
  pc <- ncol(Xc)
  
  for(trial in 1:R){
    
    set.seed(100 + trial)
    nc <- ncol(Zc)
    overall_noise <- rnorm(2*nc, 0, sig^2)
    indiv_noise <- rnorm(nc, 0, tau^2)
    xc <- Xc %*% uc + Zc %*% indiv_noise + overall_noise
    
    # posterior for control arm -- TO CHANGE DEPENDING ON METHOD
    res_inter <- get_control_prost(Xc, Zc, xc, xh, nc, nh, N, pc)
    muc <- res_inter$prost_samples
    
    res <- 0
    for (i in 1:(pc - 1)){
      curv = mean(muc[,i] <= muc[,pc])
      curv = max(curv, 1 - curv) # 2 sided
      res <- res + ifelse(curv >= cutoff, 1, 0)
    }
    power <- power + ifelse(res == pc - 1, 1, 0)
    
    # effective sample size:
    ess <- c(ess, res_inter$ess)
    
    # Count if the true value is within the quantile interval
    quantile_interval <- quantile(muc[,pc], probs = c(0.025, 0.975))
    quantile_interval_count <- quantile_interval_count + 
      ifelse((quantile_interval[1] <= uc[pc]) && (quantile_interval[2] >= uc[pc]), 1, 0)
    # Width of the quantile interval
    width_quantile_interval <- c(width_quantile_interval, 
                                 quantile_interval[2] - quantile_interval[1])
    # point estimator
    val_point_est <- mean(muc[,pc])
    point_est <- c(point_est, val_point_est)
    
    # distrns
    if (trial == R){
      distr_at <- muc
      distr_bf <- matrix(rep(0, N*pc), nrow = N)
      for (i in 1:pc){
        xt_cur = xc[which(Xc[,i] == 1)]
        t.par1 <- mean(xt_cur)
        t.par2 <- var(xt_cur)/length(xt_cur)
        distr_bf[,i] <- rst(N, t.par1, sqrt(t.par2), length(xt_cur) - 1)
      }
    }
  }
  timeend <- Sys.time()
  EHSS <- mean(ess)
  quantile_interval_count_mean <-  quantile_interval_count/R
  bias_point_est <- mean(point_est) - uc[pc]
  var_point_est <- var(point_est)
  mse_point_est <- bias_point_est^2 + var_point_est
  width_quantile_interval_mean <- mean(width_quantile_interval)
  cat("Power", power/R, "\n")
  cat("effective historical sample size is", formatC(EHSS, digits = 2, format = "f"), sep = " ", "\n")
  cat("Mean Width of Credible Interval for Control Prior", formatC(width_quantile_interval_mean, digits = 4, format = "f"), sep = " ", "\n")
  cat("% of times muc is in quantile interval is", formatC(quantile_interval_count_mean*100, digits = 4, format = "f"), sep = " ", "\n")
  cat("Bias of point estiamtor based on control prior", formatC(bias_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("Variance value of point estiamtor based on control prior", formatC(var_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("MSE of point estiamtor based on control prior", formatC(mse_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("total time for", R, "simulations is", formatC(timeend - timestart, digits = 4, format = "f"), sep = " ", "\n")
  return(list(power = power/R,
              EHSS = EHSS, 
              width_quantile_interval_mean = width_quantile_interval_mean, 
              quantile_interval_count_mean = quantile_interval_count_mean,
              bias_point_est = bias_point_est, var_point_est = var_point_est, 
              mse_point_est = mse_point_est,
              time_diff = timeend - timestart, 
              distr_at = distr_at, distr_bf = distr_bf))
  
  
}


get_control_prost <- function(Xc, Zc, xc, xh, nc, nh, N, pc) {
  
  # Logit and inverse logit functions
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  # Define the prior distribution for a0
  log_prior_a0 <- function(a0) {
    return(dbeta(a0, 0.1, 0.1, log = T))
  }
  
  # Define the posterior distribution for control arm using commensurate power prior
  log_posterior <- function(theta) {
    mu <- theta[1:pc]
    sigma2 <- exp(theta[pc + 1])
    tau2 <- exp(theta[pc + 2])
    a0 <- inv_logit(theta[pc + 3])
    
    # Likelihood
    log_lik_current <- dmvnorm(t(xc), mean = Xc %*% mu, sigma = tau2*(Zc %*% t(Zc)) + sigma2*diag(2*nc), log = T)
    log_lik_hist <- sum(dnorm(xh, mu[pc], sqrt(sigma2 + tau2), log = T))
    
    # Prior distributions
    log_prior_sigma <- -sigma2
    log_prior_tau <- -tau2
    log_prior_a0_val <- log_prior_a0(a0)
    
    integ_log <- -a0/(2*sigma2)*sum(xh^2) + nh*a0*mean(xh)^2/(2*sigma2) - log(sqrt(sigma2)/sqrt(nh*a0))
    jacob_log <- log(sigma2) + log(tau2) + log(a0*(1 - a0))
    
    # Posterior
    log_post <- log_lik_current + (a0 * log_lik_hist) + log_prior_sigma + 
      log_prior_a0_val + log_prior_tau + jacob_log #- integ_log
    
    if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
      return(-Inf)  # Return -Inf for invalid log-posterior values
    }
    return(log_post)
  }
  
  init_values <- c(rep(0, pc), log(1), log(1), logit(0.8))  # Initial values for MCMC sampling
  # Perform MCMC sampling for posterior
  samples <<- MCMCmetrop1R(log_posterior, theta.init = init_values, tune = 1,
                           mcmc = N, burnin = 1000, thin = 1, force.samp = T,
                           optim.method = "Nelder-Mead")
  
  # Extract posterior samples for mu, sigma, a0, theta0, and tau
  mu_samples <- samples[, 1:pc]
  sigma2_samples <- exp(samples[, pc + 1])
  tau2_samples <- inv_logit(samples[, pc + 2])
  
  # Calculate effective sample size (ESS)
  ess <- mean(sigma2_samples)/var(mu_samples[,pc])
  
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
  res <- matrix(rep(0, 2*n*p), nrow = 2*n)
  pairs_c <- t(combn(1:3, 2))
  for (i in 1:n){
    res[2*i - 1, pairs_c[(i %% nrow(pairs_c)) + 1, 1]] <- 1
    res[2*i, pairs_c[(i %% nrow(pairs_c)) + 1, 2]] <- 1
  }
  return(res)
}

## settings
nc <- 99
pc <- 3
uc <- c(2.41, 2.45, 2.59)
Zc <- gen_Z(nc)
Xc <- gen_X(nc, pc)
sig = 0.2
tau = 0.1
uh = 2.67
nh = 50

res <- run_simulation(sig, tau, uh, nh, uc, Xc, Zc, 5000, 10, 0.95)


