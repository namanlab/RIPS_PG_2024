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


# Define the Bayesian hypothesis testing functions
compute_oc_elastic <- function(xc, xh, xt, N) {

  
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
  
  get_params <- function(xc, xh, xt){
    nh <- length(xh)
    nc <- length(xc)
    nt <- length(xt)
    params <- decide_para(c=1, xh, nh, nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 5000)
    return(params)
  }
  
  sample_poster <- function(x0, n0, xc, nc, gt, sim=20000, nburn=10000) {
    mean_x0 <- mean(x0)
    var_x0 <- var(x0)
    mean_xc <- mean(xc)
    var_xc <- var(xc)
    D <- var_x0 / (n0 * gt)
    alpha <- (1 + nc) / 2
    muc <- rep(mean_xc, sim)
    beta <- nc * (var_xc + (mean_xc - muc)^2) / 2
    sig <- rinvgamma(sim, alpha, beta)
    mu <- (nc * mean_xc * D + sig * mean_x0) / (nc * D + sig)
    var <- sig * D / (D * nc + sig)
    muc <- rnorm(sim, mu, sqrt(var))
    muc_post <- muc[(nburn + 1):sim]
    return(list(muc_post = muc_post))
  }
  
  get_control_prost <- function(params, xc, xh, N){
    nh <- length(xh)
    nc <- length(xc)
    a <- params$a 
    b <- params$b 
    c <- params$c
    # calculate statistic and g(t)
    sp <- ((nh-1)*var(xh) + (nc-1)*var(xc))/(nh + nc - 2) # pooled variance
    cong_measure <- max(nh, nc)^(-1/4)*abs(mean(xh)-mean(xc))/(sqrt(sp/nh + sp/nc))
    gt <- 1/(1 + exp(a + b*(log(cong_measure))^c))
    if(gt == 0) gt <- 0.00001 # numerical stability
    # posterior for control arm
    temp <- sample_poster(xh, nh, xc, nc, gt, sim = N + 10000, nburn=10000)
    res <- list(ess = gt*nh + nc, prost_samples = temp$muc_post) # 
    return(res)
  }
  
  # Control
  params <- get_params(xc, xh, xt)
  res_inter <- get_control_prost(params, xc, xh, N)
  muc <- res_inter$prost_samples
  
  # Treatment
  t.par1 <- mean(xt)
  t.par2 <- var(xt)/length(xt)
  mut <- rst(N, t.par1, sqrt(t.par2), length(xt) - 1)
  
  # Probability that treatment is superior to control
  pp <- max(mean(mut <= muc), 1 - mean(mut <= muc))
  return(c(pp, res_inter$ess))
  
  
}




compute_oc_normalized <- function(xc, xh, xt, N) {
  
  get_params <- function(xc, xh, xt) {
    list()
  }
  
  get_control_prost <- function(params, xc, xh, N) {
    
    nh <- length(xh)
    nc <- length(xc)
    
    # Define the prior distribution for a0
    log_prior_a0 <- function(a0) {
      return(dbeta(a0, 1, 1, log = T))
    }
    
    # Logit and inverse logit functions
    logit <- function(p) log(p / (1 - p))
    inv_logit <- function(x) exp(x) / (1 + exp(x))
    
    # Define the posterior distribution for control arm using commensurate power prior
    log_posterior <- function(theta) {
      mu <- theta[1]
      sigma2 <- exp(theta[2])
      a0 <- inv_logit(theta[3])
      
      # Likelihood
      log_lik_current <- sum(dnorm(xc, mu, sqrt(sigma2), log = TRUE))
      log_lik_hist <- sum(dnorm(xh, mu, sqrt(sigma2), log = TRUE))
      
      # Prior distributions
      log_prior_a0_val <- log_prior_a0(a0)
      log_prior_sigma <- -sigma2
      
      # Commensurate prior for θ
      integ_log <- -a0/(2*sigma2)*sum(xh^2) + nh*a0*mean(xh)^2/(2*sigma2) - log(sqrt(sigma2)/sqrt(nh*a0))
      jacob_log <- log(sigma2) + log(a0*(1 - a0))
      
      # Posterior
      log_post <- log_lik_current + (a0 * log_lik_hist) + log_prior_a0_val + log_prior_sigma - integ_log + jacob_log
      
      if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
        return(-Inf)  # Return -Inf for invalid log-posterior values
      }
      return(log_post)
    }
    
    init_values <- c(mean(xc), log(var(xc)), logit(0.5))  # Initial values for MCMC sampling
    # Perform MCMC sampling for posterior
    R.utils::captureOutput(expr={
      samples <<- MCMCmetrop1R(log_posterior, theta.init = init_values, tune = 1,
                               mcmc = N, burnin = 1000, thin = 1, force.samp = T,
                               optim.method = "Nelder-Mead")
    })
    
    # Extract posterior samples for mu, sigma, a0, theta0, and tau
    mu_samples <- samples[, 1]
    sigma_samples <- exp(samples[, 2])
    a0_samples <- inv_logit(samples[, 3])
    
    # Calculate effective sample size (ESS)
    ess <- mean(sigma_samples)/var(mu_samples)
    
    list(prost_samples = mu_samples, ess = ess, mss = mean(sigma_samples))
  }
  
  # Control
  params <- get_params(xc, xh, xt)
  res_inter <- get_control_prost(params, xc, xh, N)
  muc <- res_inter$prost_samples
  
  # Treatment
  t.par1 <- mean(xt)
  t.par2 <- var(xt)/length(xt)
  mut <- rst(N, t.par1, sqrt(t.par2), length(xt) - 1)
  
  # Probability that treatment is superior to control
  pp <- max(mean(mut <= muc), 1 - mean(mut <= muc))
  return(c(pp, res_inter$ess))
  
  
}

compute_oc_commensurate <- function(xc, xh, xt, N) {
  
  get_params <- function(xc, xh, xt) {
    list()
  }
  
  get_control_prost <- function(params, xc, xh, N) {
    
    nh <- length(xh)
    nc <- length(xc)
    
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
      mu <- theta[1]
      sigma2 <- exp(theta[2])
      a0 <- inv_logit(theta[3])
      theta0 <- theta[4]
      tau <- exp(theta[5])
      
      # Likelihood
      log_lik_current <- sum(dnorm(xc, mu, sqrt(sigma2), log = TRUE))
      log_lik_hist <- sum(dnorm(xh, theta0, sqrt(sigma2), log = TRUE))
      
      # Prior distributions
      log_prior_tau_val <- log_prior_tau(tau)
      log_prior_a0_val <- log_prior_a0(a0, tau)
      log_prior_sigma2 <- -sigma2
      
      # Commensurate prior for θ
      log_prior_mu <- -0.5 * tau * (mu - theta0)^2
      
      integ_log <- -a0/(2*sigma2)*sum(xh^2) + nh*a0*mean(xh)^2/(2*sigma2) - log(sqrt(sigma2)/sqrt(nh*a0))
      jacob_log <- log(sigma2) + log(a0*(1 - a0)) + log(tau)
      
      # Posterior
      log_post <- log_lik_current + (a0 * log_lik_hist) + log_prior_a0_val + log_prior_tau_val + log_prior_mu + log_prior_sigma2 - integ_log + jacob_log
      
      if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
        return(-Inf)  # Return -Inf for invalid log-posterior values
      }
      return(log_post)
    }
    
    init_values <- c(mean(xc), log(var(xc)), logit(0.5), mean(xh), log(1))  # Initial values for MCMC sampling
    # Perform MCMC sampling for posterior
    R.utils::captureOutput(expr={
      samples <- MCMCmetrop1R(log_posterior, theta.init = init_values, 
                              mcmc = N, burnin = 1000, thin = 1, force.samp = T)
    })
    
    # Extract posterior samples for mu, sigma, a0, theta0, and tau
    mu_samples <- samples[, 1]
    sigma_samples <- exp(samples[, 2])
    a0_samples <- inv_logit(samples[, 3])
    theta0_samples <- samples[, 4]
    tau_samples <- exp(samples[, 5])
    
    # Calculate effective sample size (ESS)
    num <- mean(sigma_samples)/(length(xc) + length(xh) + mean(sigma_samples)*mean(tau_samples))
    ess <- num/var(mu_samples)*(length(xc) + length(xh))
    
    list(prost_samples = mu_samples, ess = ess, mss = mean(sigma_samples))
  }
  
  # Control
  params <- get_params(xc, xh, xt)
  res_inter <- get_control_prost(params, xc, xh, N)
  muc <- res_inter$prost_samples
  
  # Treatment
  t.par1 <- mean(xt)
  t.par2 <- var(xt)/length(xt)
  mut <- rst(N, t.par1, sqrt(t.par2), length(xt) - 1)
  
  # Probability that treatment is superior to control
  pp <- max(mean(mut <= muc), 1 - mean(mut <= muc))
  return(c(pp, res_inter$ess))
  
  
}

compute_oc_robust_map <- function(xc, xh, xt, N) {
  
  get_params <- function(xc, xh, xt) {
    list()
  }
  
  get_control_prost <- function(params, xc, xh, N) {
    
    nh <- length(xh)
    nc <- length(xc)
    
    # Logit and inverse logit functions
    logit <- function(p) log(p / (1 - p))
    inv_logit <- function(x) exp(x) / (1 + exp(x))
    prob_w <- pnorm(mean(xc), mean = mean(xh), sd = sd(xh)/sqrt(nh))
    w = 2*min(prob_w, 1 - prob_w)
    
    # Define the posterior distribution for control arm using commensurate power prior
    log_posterior <- function(theta) {
      mu <- theta[1]
      sigma2 <- exp(theta[2])
      
      # Likelihood
      log_lik_current <- sum(dnorm(xc, mu, sqrt(sigma2), log = TRUE))
      log_lik_hist <- dnorm(mean(xh), mu, sd(xh)/sqrt(nh))
      log_lik_theta_hist <- log(w*log_lik_hist + (1 - w))
      
      # Prior distributions
      log_prior_sigma <- -sigma2
      
      # Commensurate prior for θ
      jacob_log <- log(sigma2)
      
      # Posterior
      log_post <- log_lik_current + log_lik_theta_hist + log_prior_sigma + jacob_log
      
      if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
        return(-Inf)  # Return -Inf for invalid log-posterior values
      }
      return(log_post)
    }
    
    init_values <- c(mean(xc), log(var(xc)))  # Initial values for MCMC sampling
    # Perform MCMC sampling for posterior
    R.utils::captureOutput(expr={
      samples <<- MCMCmetrop1R(log_posterior, theta.init = init_values, tune = 1,
                               mcmc = N, burnin = 1000, thin = 1, force.samp = T,
                               optim.method = "Nelder-Mead")
    })
    # Extract posterior samples for mu, sigma, a0, theta0, and tau
    mu_samples <- samples[, 1]
    sigma_samples <- exp(samples[, 2])
    
    # Calculate effective sample size (ESS)
    ess <- mean(sigma_samples)/var(mu_samples)
    
    list(prost_samples = mu_samples, ess = ess, mss = mean(sigma_samples))
  }
  
  # Control
  params <- get_params(xc, xh, xt)
  res_inter <- get_control_prost(params, xc, xh, N)
  muc <- res_inter$prost_samples
  
  # Treatment
  t.par1 <- mean(xt)
  t.par2 <- var(xt)/length(xt)
  mut <- rst(N, t.par1, sqrt(t.par2), length(xt) - 1)
  
  # Probability that treatment is superior to control
  pp <- max(mean(mut <= muc), 1 - mean(mut <= muc))
  return(c(pp, res_inter$ess))
  
}


compute_oc_elastic_power <- function(xc, xh, xt, N) {
  
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
  
  get_params <- function(xc, xh, xt) {
    list()
  }
  
  get_GT <- function(xc, xh, nc, nh, N){
    params <- decide_para(c=1, xh, nh, nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 5000)
    a <- params$a 
    b <- params$b 
    c <- params$c
    # calculate statistic and g(t)
    sp <- ((nh-1)*var(xh) + (nc-1)*var(xc))/(nh + nc - 2) # pooled variance
    cong_measure <- max(nh, nc)^(-1/4)*abs(mean(xh)-mean(xc))/(sqrt(sp/nh + sp/nc))
    gt <- 1/(1 + exp(a + b*(log(cong_measure))^c))
    if(gt == 0) gt <- 0.00001 # numerical stability
    return(gt)
  }
  
  get_control_prost <- function(params, xc, xh, N) {
    
    nh <- length(xh)
    nc <- length(xc)
    
    a0 <- get_GT(xc, xh, nc, nh, N)
    
    # Logit and inverse logit functions
    logit <- function(p) log(p / (1 - p))
    inv_logit <- function(x) exp(x) / (1 + exp(x))
    
    # Define the posterior distribution for control arm using commensurate power prior
    log_posterior <- function(theta) {
      mu <- theta[1]
      sigma2 <- exp(theta[2])
      
      # Likelihood
      log_lik_current <- sum(dnorm(xc, mu, sqrt(sigma2), log = TRUE))
      log_lik_hist <- sum(dnorm(xh, mu, sqrt(sigma2), log = TRUE))
      
      # Prior distributions
      log_prior_sigma <- -sigma2
      
      # Commensurate prior for θ
      jacob_log <- log(sigma2)
      
      # Posterior
      log_post <- log_lik_current + (a0 * log_lik_hist) + log_prior_sigma + jacob_log
      
      if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
        return(-Inf)  # Return -Inf for invalid log-posterior values
      }
      return(log_post)
    }
    
    init_values <- c(mean(xc), log(var(xc)))  # Initial values for MCMC sampling
    # Perform MCMC sampling for posterior
    R.utils::captureOutput(expr={
      samples <<- MCMCmetrop1R(log_posterior, theta.init = init_values, tune = 1,
                               mcmc = N, burnin = 1000, thin = 1, force.samp = T,
                               optim.method = "Nelder-Mead")
    })
    # Extract posterior samples for mu, sigma, a0, theta0, and tau
    mu_samples <- samples[, 1]
    sigma_samples <- exp(samples[, 2])
    
    # Calculate effective sample size (ESS)
    ess <- mean(sigma_samples)/var(mu_samples)
    
    list(prost_samples = mu_samples, ess = ess, mss = mean(sigma_samples))
  }
  
  # Control
  params <- get_params(xc, xh, xt)
  res_inter <- get_control_prost(params, xc, xh, N)
  muc <- res_inter$prost_samples
  
  # Treatment
  t.par1 <- mean(xt)
  t.par2 <- var(xt)/length(xt)
  mut <- rst(N, t.par1, sqrt(t.par2), length(xt) - 1)
  
  # Probability that treatment is superior to control
  pp <- max(mean(mut <= muc), 1 - mean(mut <= muc))
  return(c(pp, res_inter$ess))
  
  
}

# Define the frequentist hypothesis testing function
compute_oc_frequentist <- function(yc, yt) {
  n_c <- length(yc)
  n_t <- length(yt)
  mean_c <- mean(yc)
  mean_t <- mean(yt)
  var_c <- var(yc)
  var_t <- var(yt)
  
  pooled_var <- ((n_c - 1) * var_c + (n_t - 1) * var_t) / (n_c + n_t - 2)
  t_stat <- (mean_t - mean_c) / sqrt(pooled_var * (1/n_c + 1/n_t))
  df <- n_c + n_t - 2
  
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  
  return(p_value)
}











# library(writexl)
# set.seed(123)
# yc <- rnorm(30, mean = 10, sd = 2)
# yt <- rnorm(30, mean = 11, sd = 2)
# yh <- rnorm(120, mean = 10.5, sd = 2)
# data <- data.frame(yc = yc, yt = yt, yh = yh)
# write_xlsx(data, path = "sample_generated_data.xlsx")
