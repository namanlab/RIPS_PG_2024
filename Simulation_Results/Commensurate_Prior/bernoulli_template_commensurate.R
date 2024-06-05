#################
################# EP1 with commensurate power prior #######################
#################

library(MCMCpack)
library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)
# ADD ANY OTHER PACKAGES

run_simulation <- function(nt, nc, nh, pc, pt, ph, H = 1, N, R, cutoff){
  
  set.seed(10)
  xh <- rbinom(1, nh, ph) # historical data
  
  # Metrics to calculate
  rej_null <- 0 # number of rejections
  ess <- NULL # effective sample size
  timestart <- Sys.time() # computation time
  quantile_interval_count <- 0 # number of times muc is in quantile interval
  width_quantile_interval <- NULL # width of credible interval
  point_est <- NULL # point estimator based on control prior
  xc_act <- NULL # actual control samples
  xc_samp <- NULL # samples based on prosterior predictive
  
  # any params for models to be computed a-priori -- TO CHANGE DEPENDING ON METHOD
  params <- get_params(xh, nt, nc, nh, pc, pt, ph, H, N, R)
  
  for(trial in 1:R){
    
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rbinom(1, nc, pc)
    xt <- rbinom(1, nt, pt)
    
    # posterior for control arm -- TO CHANGE DEPENIDNG ON METHOD
    res_inter <- get_control_prost(params, xc, xh, nc, nh, H, N)
    pc_post <- res_inter$prost_samples
    
    # posterior for treatment arm
    pt_post <- rbeta(N, params$t_alpha + xt, params$t_beta + nt - xt)
    
    # METRICS:
    # probability that treatment is superior to control
    pp <- mean(pt_post > pc_post)
    # number of rejections of null
    if(pp >= cutoff){
      rej_null <- rej_null + 1
    }
    # effective sample size:
    ess <- c(ess, res_inter$ess)
    # Count if the true value is within the quantile interval
    quantile_interval <- quantile(pc_post, probs = c(0.025, 0.975))
    quantile_interval_count <- quantile_interval_count + ifelse((quantile_interval[1] <= pc) && (quantile_interval[2] >= pc), 1, 0)
    # Width of the quantile interval
    width_quantile_interval <- c(width_quantile_interval, quantile_interval[2] - quantile_interval[1])
    # point estimator
    val_point_est <- mean(pc_post)
    point_est <- c(point_est, val_point_est)
    # posterior predictive check:
    xc_act <- c(xc_act, rbinom(length(pc_post), nc, pc))
    xc_samp <- c(xc_samp, rbinom(length(pc_post), nc, pc_post))
  }
  timeend <- Sys.time()
  EHSS <- mean(ess)
  prob_rej <- rej_null/R
  quantile_interval_count_mean <-  quantile_interval_count/R
  bias_point_est <- mean(point_est) - uc
  var_point_est <- var(point_est)
  mse_point_est <- bias_point_est^2 + var_point_est
  width_quantile_interval_mean <- mean(width_quantile_interval)
  cat("probability of claiming efficacy is", prob_rej, "\n")
  cat("effective historical sample size is", formatC(EHSS, digits = 2, format = "f"), sep = " ", "\n")
  cat("Mean Width of Credible Interval for Control Prior", formatC(width_quantile_interval_mean, digits = 4, format = "f"), sep = " ", "\n")
  cat("% of times pc is in quantile interval is", formatC(quantile_interval_count_mean*100, digits = 4, format = "f"), sep = " ", "\n")
  cat("Bias of point estiamtor based on control prior", formatC(bias_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("Variance value of point estiamtor based on control prior", formatC(var_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("MSE of point estiamtor based on control prior", formatC(mse_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("total time for", R, "simulations is", formatC(timeend - timestart, digits = 4, format = "f"), sep = " ", "\n")
  plot_comp <- ggplot(data = tibble("Control Distribution" = xc_act, 
                                    "Prost Predictive" = xc_samp) %>% 
                        pivot_longer(1:2, names_to = "type", values_to = "val") ) +
    geom_histogram(aes(x = val, y = after_stat(density), fill = type), alpha = 0.7, binwidth = 1) +
    theme_bw() + scale_fill_manual(values = c("red", "blue"))
  return(list(prob.rej = prob_rej, width_quantile_interval_mean = width_quantile_interval_mean, 
              quantile_interval_count_mean = quantile_interval_count_mean,
              bias_point_est = bias_point_est, var_point_est = var_point_est, mse_point_est = mse_point_est,
              time_diff = timeend - timestart, plot_comp = plot_comp))
  
}


# can return a list of all params if needed to be used by model, return null if not needed
get_params <- function(xh, nt, nc, nh, pc, pt, ph, H, N, R){
  params <- list()
  params$t_alpha = 0.1
  params$t_beta = 0.1
  params$c_alpha = 0.1
  params$c_beta = 0.1
  return(params)
}


# Function to perform MCMC sampling for the posterior using power prior
get_control_prost <- function(params, xc, xh, nc, nh, H, N) {
  
  # Define the prior distribution for a0 depending on tau
  log_prior_a0 <- function(a0, tau) {
    g_tau <- function(tau) {
      return(max(tau, 1))
    }
    return((g_tau(tau) - 1) * log(a0))
  }
  
  # Define the prior distribution for tau
  log_prior_tau <- function(tau) {
    return(ifelse(tau <= 30 & tau >= 0, 1, 0))
  }
  
  # Logit and inverse logit functions
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  # Define the posterior distribution for control arm using commensurate power prior
  log_posterior <- function(theta) {
    p <- inv_logit(theta[1])
    a0 <- inv_logit(theta[2])
    theta0 <- inv_logit(theta[3])
    tau <- exp(theta[4])
    
    # Likelihood
    log_lik_current <- dbinom(xc, nc, p, log = TRUE)
    log_lik_hist <- dbinom(xh, nh, theta0, log = TRUE)
    
    # Prior distributions
    log_prior_tau_val <- log_prior_tau(tau)
    log_prior_a0_val <- log_prior_a0(a0, tau)
    log_prior_p0 <- dbeta(theta0, params$c_alpha, params$c_beta, log = TRUE)
    
    # Commensurate prior for p
    anew = theta0*1/tau
    bnew = (1 - theta0)*1/tau
    log_prior_p <- dbeta(p, anew, bnew, log = T)
    
    # Commensurate prior for p
    integ_log <- log(beta(a0*xh + 1, a0*(nh - xh) + 1))
    jacob_log <- log(p*(1 - p)) + log(a0*(1 - a0)) + log(theta0*(1 - theta0)) + log(tau)
    
    # Posterior
    log_post <- log_lik_current + (a0 * log_lik_hist) + log_prior_p + log_prior_tau_val + log_prior_a0_val - integ_log + jacob_log + log_prior_p0
    
    if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
      return(-Inf)  # Return -Inf for invalid log-posterior values
    }
    return(log_post)
  }
  init_values <- c(logit(mean(xc/nc)), logit(0.5), logit(mean(xh/nh)), log(2))
  # Perform MCMC sampling
  samples <- MCMCmetrop1R(log_posterior, theta.init = init_values, mcmc = N, burnin = 1000, thin = 10, force.samp = T)
  # Extract posterior samples for p
  logit_p_samples <- samples[, 1]
  p_samples <- inv_logit(logit_p_samples)
  
  # log prior for ess:
  log_prior <- function(theta) {
    p <- inv_logit(theta[1])
    a0 <- inv_logit(theta[2])
    theta0 <- inv_logit(theta[3])
    tau <- exp(theta[4])
    
    # Likelihood
    log_lik_hist <- dbinom(xh, nh, theta0, log = TRUE)
    
    # Prior distributions
    log_prior_tau_val <- log_prior_tau(tau)
    log_prior_a0_val <- log_prior_a0(a0, tau)
    log_prior_p0 <- dbeta(theta0, params$c_alpha, params$c_beta, log = TRUE)
    
    # Commensurate prior for p
    anew = theta0*1/tau
    bnew = (1 - theta0)*1/tau
    log_prior_p <- dbeta(p, anew, bnew, log = T)
    
    # Commensurate prior for p
    integ_log <- log(beta(a0*xh + 1, a0*(nh - xh) + 1))
    jacob_log <- log(p*(1 - p)) + log(a0*(1 - a0)) + log(theta0*(1 - theta0)) + log(tau)
    
    # Posterior
    log_post <- (a0 * log_lik_hist) + log_prior_p + log_prior_tau_val + log_prior_a0_val - integ_log + jacob_log + log_prior_p0
    
    if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
      return(-Inf)  # Return -Inf for invalid log-posterior values
    }
    return(log_post)
  }
  init_values <- c(logit(mean(xc/nc)), logit(0.5), logit(mean(xh/nh)), log(2))
  # Perform MCMC sampling
  samples <- MCMCmetrop1R(log_prior, theta.init = init_values, mcmc = N, burnin = 1000, thin = 10, force.samp = T)
  logit_p_samples <- samples[, 1]
  p_samples_prior <- inv_logit(logit_p_samples)
  ess <- mean(xh/nh)*(1 - mean(xh/nh))/(var(p_samples_prior))
  list(prost_samples = p_samples, ess = ess)
}

#########--------------------------------------------------------------#########
############################ NO CHANGE TO BE MADE  #############################
#########--------------------------------------------------------------#########



## settings
nc <- 25 # current control size
nt <- 50 # current treatment size
nh <- 50 # historical control size
pc <- 0.4 # true mean of control

# strong congruence between control and historical
res1 <- run_simulation(nt, nc, nh, pc, pt = 0.4, ph = 0.4, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res1$plot_comp
res2 <- run_simulation(nt, nc, nh, pc, pt = 0.6, ph = 0.4, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res2$plot_comp

# weak congruence between control and historical
res3 <- run_simulation(nt, nc, nh, pc, pt = 0.4, ph = 0.45, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res3$plot_comp
res4 <- run_simulation(nt, nc, nh, pc, pt = 0.6, ph = 0.45, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res4$plot_comp

# no congruence between control and historical
res5 <- run_simulation(nt, nc, nh, pc, pt = 0.4, ph = 0.6, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res5$plot_comp
res6 <- run_simulation(nt, nc, nh, pc, pt = 0.6, ph = 0.6, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res6$plot_comp

# Combine results into a list and Save the list as an RDS file
results <- list(res1 = res1, res2 = res2, res3 = res3, res4 = res4, res5 = res5, res6 = res6)
saveRDS(results, file = "simulation_results_elastic_prior_bernoulli.rds")

# Load the results: results <- readRDS("simulation_results.rds")