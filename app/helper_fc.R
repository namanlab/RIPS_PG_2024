library(shiny)
library(shinydashboard)
library(latex2exp)
library(tidyverse)
library(e1071)
library(stringr)
library(plotly)
library(reshape2)
library(bslib)
library(MCMCpack)
library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)
library(mvtnorm)

gen_Z <- function(n){
  bm <- diag(n)
  res <- matrix(rep(0, 2*n^2), nrow = 2*n)
  for (i in 1:(2*n)){
    print(bm[floor(i/2),])
    res[i,] <- bm[ceiling(i/2),]
  }
  return(res)
}

gen_X <- function(n, p, num_ones_last_col) {
  # Validate the num_ones_last_col argument
  if (num_ones_last_col < 0 || num_ones_last_col > 2 * n) {
    stop("num_ones_last_col should be between 0 and 2 * n.")
  }
  
  # Initialize the result matrix
  res <- matrix(0, nrow = 2 * n, ncol = p)
  pairs_c <- t(combn(1:p, 2))
  num_pairs <- nrow(pairs_c)
  
  # Generate pairs normally
  for (i in 1:n) {
    pair_idx <- (i %% num_pairs) + 1
    res[2 * i - 1, pairs_c[pair_idx, 1]] <- 1
    res[2 * i, pairs_c[pair_idx, 2]] <- 1
  }
  
  # Adjust the last column to have exactly num_ones_last_col 1s
  current_ones_last_col <- sum(res[, p])
  if (current_ones_last_col > num_ones_last_col) {
    # Too many 1s, remove some
    ones_indices <- which(res[, p] == 1)
    to_remove <- sample(ones_indices, current_ones_last_col - num_ones_last_col)
    res[to_remove, p] <- 0
    
    # Replace the removed 1s with 1s in other columns
    set.seed(42)
    for (idx in to_remove) {
      paired_row <- ifelse(idx %% 2 == 1, idx + 1, idx - 1)
      possible_cols <- setdiff(1:(p-1), which(res[paired_row, ] == 1))
      replace_with_col <- ifelse(length(possible_cols) == 1, possible_cols, sample(possible_cols, 1))
      res[idx, replace_with_col] <- 1
    }
  } 
  
  return(res)
}

get_control_prost <- function(Xc, Zc, xc, xh, nc, nh, N, pc) {
  
  xcc <- xc[which(Xc[,pc] == 1)]
  ncc <- length(xcc)
  gt <- getGT(xcc, xh, ncc, nh, N, pc)
  
  # Logit and inverse logit functions
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  # Define the posterior distribution for control arm using commensurate power prior
  log_posterior <- function(theta) {
    mu <- theta[1:pc]
    sigma2 <- exp(theta[pc + 1])
    tau2 <- exp(theta[pc + 2])
    
    # Likelihood
    log_lik_current <- dmvnorm(t(xc), mean = Xc %*% mu, sigma = tau2*(Zc %*% t(Zc)) + sigma2*diag(2*nc), log = T)
    log_lik_hist <- dnorm(mu[pc], mean(xh), sqrt((sigma2 + tau2)/(nh*gt)), log = T)
    
    # Prior distributions
    log_prior_sigma <- -sigma2
    log_prior_tau <- -tau2
    
    jacob_log <- log(sigma2) + log(tau2) 
    
    # Posterior
    log_post <- log_lik_current + log_lik_hist + log_prior_sigma + log_prior_tau 
    
    if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
      return(-Inf)  # Return -Inf for invalid log-posterior values
    }
    return(log_post)
  }
  
  init_mu <- rep(0, pc)
  for (i in 1:pc){init_mu[i] <- mean(xc[which(Xc[,i] == 1)])}
  init_values <- c(init_mu, log(1), log(1))  # Initial values for MCMC sampling
  # Perform MCMC sampling for posterior
  samples <<- MCMCmetrop1R(log_posterior, theta.init = init_values, tune = 1,
                           mcmc = N, burnin = 1000, thin = 1, force.samp = T,
                           optim.method = "Nelder-Mead")
  
  # Extract posterior samples for mu, sigma, a0, theta0, and tau
  mu_samples <- samples[, 1:pc]
  sigma2_samples <- exp(samples[, pc + 1])
  tau2_samples <- inv_logit(samples[, pc + 2])
  list(prost_samples = mu_samples, ess = gt*nh)
}

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
  return(list(a=a, b=b, c=c))
}

getGT <- function(xc, xh, nc, nh, N, pc){
  params <- decide_para(c=1, xh, nh, nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 5000)
  a <- params$a 
  b <- params$b 
  c <- params$c
  # calculate statistic and g(t)
  sp <- ((nh-1)*var(xh) + (nc-1)*var(xc))/(nh + nc - 2) # pooled variance
  cong_measure <- max(nh, nc)^(-1/4)*abs(mean(xh)-mean(xc))/(sqrt(sp/nh + sp/nc))
  gt <- 1/(1 + exp(a + b*(log(cong_measure))^c))
  return(gt)
}


# Define the Bayesian hypothesis testing functions
compute_fc_elastic <- function(Xc, Zc, xc, xh, N) {
  
  nh <- length(xh)
  pc <- ncol(Xc)
  nc <- ncol(Zc)
  
  res_inter <- get_control_prost(Xc, Zc, xc, xh, nc, nh, N, pc)
  muc <- res_inter$prost_samples
  
  res <- list()
  for (i in 1:(pc - 1)){
    curv = mean(muc[,i] <= muc[,pc])
    curv = max(curv, 1 - curv) # 2 sided
    res[[i]] = curv
  }
  print(res)
  return(res)
}

