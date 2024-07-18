library(MCMCpack)     # For MCMC sampling
library(LaplacesDemon) # For inverse logit functions and more
library(invgamma)      # For inverse gamma distribution functions
library(tidyverse)     # For data manipulation and visualization
library(HDInterval)    # For credible intervals
library(mvtnorm)       # For multivariate normal distribution functions

# Function to run simulation for the Bayesian Dynamic Borrowing
run_simulation <- function(sig, tau, uh, nh, uc, Xc, Zc, N, R, cutoff) {
  # sig: variance of overall noise
  # tau: variance of individual noise
  # uh: historical mean
  # nh: historical sample size
  # uc: true values of control mean vector
  # Xc: design matrix for control group
  # Zc: design matrix for individual-specific effects
  # N: number of MCMC samples
  # R: number of simulation runs
  # cutoff: threshold for credible interval
  
  set.seed(42)  # Set seed for reproducibility
  xh <- rnorm(nh, uh, sqrt(sig^2 + tau^2))  # Historical data
  power <- 0  # Initialize power counter
  ess <- NULL  # Effective sample size
  timestart <- Sys.time()  # Start timer for computation time
  quantile_interval_count <- 0  # Counter for quantile interval inclusion
  width_quantile_interval <- NULL  # Width of credible intervals
  point_est <- NULL  # Point estimator values
  distr_bf <- NULL  # Distribution before fitting
  distr_at <- NULL  # Distribution after fitting
  pc <- ncol(Xc)  # Number of columns in design matrix Xc
  
  for (trial in 1:R) {
    set.seed(100 + trial)  # Set seed for each trial
    nc <- ncol(Zc)  # Number of columns in design matrix Zc
    overall_noise <- rnorm(2*nc, 0, sig^2)  # Overall noise
    indiv_noise <- rnorm(nc, 0, tau^2)  # Individual noise
    xc <- Xc %*% uc + Zc %*% indiv_noise + overall_noise  # Simulated control group data
    
    # Posterior sampling for control arm
    res_inter <- get_control_prost(Xc, Zc, xc, xh, nc, nh, N, pc)
    muc <- res_inter$prost_samples  # Posterior samples for mean
    print(apply(muc, 2, mean))  # Print mean of posterior samples
    
    res <- 0  # Initialize result counter
    for (i in 1:(pc - 1)) {
      curv <- mean(muc[,i] <= muc[,pc])  # Compute posterior probability
      curv <- max(curv, 1 - curv)  # Two-sided test
      res <- res + ifelse(curv >= cutoff, 1, 0)  # Compare with cutoff
    }
    power <- power + ifelse(res == pc - 1, 1, 0)  # Update power counter
    
    ess <- c(ess, res_inter$ess)  # Collect effective sample sizes
    
    # Credible interval checks
    quantile_interval <- quantile(muc[,pc], probs = c(0.025, 0.975))  # Quantile interval
    quantile_interval_count <- quantile_interval_count + 
      ifelse((quantile_interval[1] <= uc[pc]) && (quantile_interval[2] >= uc[pc]), 1, 0)
    width_quantile_interval <- c(width_quantile_interval, 
                                 quantile_interval[2] - quantile_interval[1])
    val_point_est <- mean(muc[,pc])  # Point estimator
    point_est <- c(point_est, val_point_est)
    
    # Distributions before and after fitting
    if (trial == R) {
      distr_at <- muc
      distr_bf <- matrix(rep(0, N*pc), nrow = N)
      for (i in 1:pc) {
        xt_cur <- xc[which(Xc[,i] == 1)]
        t.par1 <- mean(xt_cur)
        t.par2 <- var(xt_cur)/length(xt_cur)
        distr_bf[,i] <- rst(N, t.par1, sqrt(t.par2), length(xt_cur) - 1)
      }
    }
  }
  
  # Compute summary statistics
  timeend <- Sys.time()
  EHSS <- mean(ess)
  quantile_interval_count_mean <- quantile_interval_count/R
  bias_point_est <- mean(point_est) - uc[pc]
  var_point_est <- var(point_est)
  mse_point_est <- bias_point_est^2 + var_point_est
  width_quantile_interval_mean <- mean(width_quantile_interval)
  
  # Print results
  cat("Power", power/R, "\n")
  cat("Effective historical sample size is", formatC(EHSS, digits = 2, format = "f"), "\n")
  cat("Mean Width of Credible Interval for Control Prior", formatC(width_quantile_interval_mean, digits = 4, format = "f"), "\n")
  cat("% of times muc is in quantile interval is", formatC(quantile_interval_count_mean*100, digits = 4, format = "f"), "\n")
  cat("Bias of point estimator based on control prior", formatC(bias_point_est, digits = 4, format = "f"), "\n")
  cat("Variance of point estimator based on control prior", formatC(var_point_est, digits = 4, format = "f"), "\n")
  cat("MSE of point estimator based on control prior", formatC(mse_point_est, digits = 4, format = "f"), "\n")
  cat("Total time for", R, "simulations is", formatC(timeend - timestart, digits = 4, format = "f"), "\n")
  
  # Prepare data for plotting
  df_plot <- NULL
  for (i in 1:pc) {
    cur_df <- tibble(Before = distr_bf[,i], After = distr_at[,i]) %>% 
      pivot_longer(1:2, names_to = "type", values_to = "val") 
    cur_df$p <- i
    df_plot <- rbind(df_plot, cur_df)
  }
  
  return(list(power = power/R,
              EHSS = EHSS, 
              width_quantile_interval_mean = width_quantile_interval_mean, 
              quantile_interval_count_mean = quantile_interval_count_mean,
              bias_point_est = bias_point_est, var_point_est = var_point_est, 
              mse_point_est = mse_point_est,
              time_diff = timeend - timestart, 
              distr_df = df_plot))
}

# Function to get control posterior samples
get_control_prost <- function(Xc, Zc, xc, xh, nc, nh, N, pc) {
  # Xc: Design matrix for control group
  # Zc: Design matrix for individual-specific effects
  # xc: Simulated control group data
  # xh: Historical data
  # nc: Number of control group samples
  # nh: Number of historical samples
  # N: Number of MCMC samples
  # pc: Number of columns in design matrix Xc
  
  log_inverse_gamma <- function(sigma2) {
    if (sigma2 > 0) {
      return(dinvgamma(sigma2, shape = 2, scale = 2, log = TRUE))
    } else {
      return(-Inf)
    }
  }
  
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
    log_prior_sigma <- log_inverse_gamma(sigma2)
    log_prior_tau <- log_inverse_gamma(tau2)
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
  
  # Calculate effective sample size (ESS)
  ess <- mean(sigma2_samples)/var(mu_samples[,pc])
  
  list(prost_samples = mu_samples, ess = ess)
}



#########--------------------------------------------------------------#########
#################################### MODEL #####################################
#########--------------------------------------------------------------#########

# Function to generate Z matrix
gen_Z <- function(n){
  bm <- diag(n)  # Create a diagonal matrix of size n
  res <- matrix(rep(0, 2*n^2), nrow = 2*n)  # Initialize result matrix
  for (i in 1:(2*n)){
    res[i,] <- bm[ceiling(i/2),]  # Assign values to result matrix based on bm matrix
  }
  return(res)
}

# Function to generate X matrix
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


## settings
nc <- 99
pc <- 3
uc <- c(2.6, 2.6, 2.6)
Zc <- gen_Z(nc)
Xc <- gen_X(nc, pc, round(2*nc/pc))
sig = 0.2
uh = 2.6
nh = 50
res = run_simulation(sig, 0.1, uh, nh, uc, Xc, Zc, 5000, 1, 0.95)
res$distr_df %>%
  ggplot(aes(x = val, fill = type)) + geom_density(alpha = 0.5) +
  facet_wrap(~p) +
  theme_bw() +
  labs(fill = "Type")

# Tau Analysis:
final_df_1 <- NULL
final_df_2 <- NULL

for (i in seq(0, 0.4, 0.02)){
  cat("\n==========\nTau:", i, "\n==========\n")
  res <- run_simulation(sig, i, uh, nh, uc, Xc, Zc, 5000, 5, 0.95)
  temp_df_1 = data.frame(tau = i, pow = res$power, ess = res$EHSS, mse = res$mse_point_est)
  final_df_1 <- rbind(final_df_1, temp_df_1)
  temp_df_2 <- res$distr_df
  temp_df_2$tau = i
  final_df_2 <- rbind(final_df_2, temp_df_2)
}

write.csv(final_df_1, "results/commensurate_results_tau1_updated.csv")
write.csv(final_df_2, "results/commensurate_results_tau2_updated.csv")



nc <- 99
pc <- 3
Zc <- gen_Z(nc)
sig = 0.2
tau = 0.1
nh = 50
final_df <- NULL
delta1 <- seq(0, 0.2, 0.02)
delta2 <- seq(0, 0.2, 0.02)
nc_seq <- seq(66, 18, -12)
for (nc_val in nc_seq){
  for (i in delta1){
    cat("\n==========\nProcessing nc:", nc_val, " delta1:", i, "\n==========\n")
    for (j in delta2){
      set.seed(42)
      temp_uc <- c(rep(2.6, pc - 1), 2.6 + rnorm(1, i, 0.02))
      temp_uh <-  2.6 + j
      temp_Xc <- gen_X(nc, pc, nc_val)
      res1 <- run_simulation(sig, tau, temp_uh, nh, temp_uc, temp_Xc, Zc, 5000, 10, 0.95)
      temp_df <- data.frame(nc = nc_val, delta1 = i, delta2 = j, pow = res1$power,
                            ess = res1$EHSS)
      final_df <- rbind(final_df, temp_df)
      
      print(res1$distr_df %>%
              ggplot(aes(x = val, fill = type)) + geom_density(alpha = 0.5) +
              facet_wrap(~p) +
              theme_bw() +
              labs(fill = "Type"))
      
      # Checkpointing
      if (nrow(final_df) %% 110 == 0) {
        write.csv(final_df, file = paste0("results/commensurate_checkpoint_fc_nc_", nrow(final_df), ".csv"))
      }
    }
  }
}
# 
write.csv(final_df, "results/commensurate_results_nc_fc_temp.csv")


