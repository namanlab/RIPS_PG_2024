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
    print(apply(muc, 2, mean))
    
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
  
  df_plot <- NULL
  for (i in 1:pc){
    cur_df <- tibble(Before = distr_bf[,i], After = distr_at[,i]) %>% 
      pivot_longer(1:2, names_to = "type", values_to = "val") 
    cur_df$p = i
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


get_control_prost <- function(Xc, Zc, xc, xh, nc, nh, N, pc) {
  
  xcc <- xc[which(Xc[,pc] == 1)]
  ncc <- length(xcc)
  gt <- getGT(xcc, xh, ncc, nh, N, pc)
  
  # Logit and inverse logit functions
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  log_inverse_gamma <- function(sigma2) {
    if (sigma2 > 0) {
      return(dinvgamma(sigma2, shape = 2, scale = 2, log = TRUE))
    } else {
      return(-Inf)
    }
  }
  
  
  # Define the posterior distribution for control arm using commensurate power prior
  log_posterior <- function(theta) {
    mu <- theta[1:pc]
    sigma2 <- exp(theta[pc + 1])
    tau2 <- exp(theta[pc + 2])
    
    # Likelihood
    log_lik_current <- dmvnorm(t(xc), mean = Xc %*% mu, sigma = tau2*(Zc %*% t(Zc)) + sigma2*diag(2*nc), log = T)
    log_lik_hist <- dnorm(mu[pc], mean(xh), sqrt((sigma2 + tau2)/(nh*gt)), log = T)
    
    # Prior distributions
    log_prior_sigma <- log_inverse_gamma(sigma2)
    log_prior_tau <- log_inverse_gamma(tau2)
    # log_prior_sigma <- -sigma2
    # log_prior_tau <- -tau2
    
    jacob_log <- log(sigma2) + log(tau2) 
    
    # Posterior
    log_post <- log_lik_current + log_lik_hist + log_prior_sigma + log_prior_tau + jacob_log 
    
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
# res = run_simulation(sig, 0.1, uh, nh, uc, Xc, Zc, 5000, 1, 0.95)
# res$distr_df %>%
#   ggplot(aes(x = val, fill = type)) + geom_density(alpha = 0.5) +
#   facet_wrap(~p) +
#   theme_bw() +
#   labs(fill = "Type")
# 
# # Tau Analysis:
# final_df_1 <- NULL
# final_df_2 <- NULL
# 
# for (i in seq(0, 0.4, 0.02)){
#   cat("\n==========\nTau:", i, "\n==========\n")
#   res <- run_simulation(sig, i, uh, nh, uc, Xc, Zc, 5000, 5, 0.95)
#   temp_df_1 = data.frame(tau = i, pow = res$power, ess = res$EHSS, mse = res$mse_point_est)
#   final_df_1 <- rbind(final_df_1, temp_df_1)
#   temp_df_2 <- res$distr_df
#   temp_df_2$tau = i
#   final_df_2 <- rbind(final_df_2, temp_df_2)
# }
# 
# write.csv(final_df_1, "results/elastic_results_tau1_updated.csv")
# write.csv(final_df_2, "results/elastic_results_tau2_updated.csv")




# library(foreach)
# library(doParallel)
# library(doSNOW)
# 
# numCores <- detectCores() - 2
# cl <- makeCluster(numCores)
# registerDoParallel(cl)
# 
# final_df <- NULL
# delta1 <- seq(0, 0.2, 0.02)
# delta2 <- seq(0, 0.2, 0.02)
# nc_seq <- seq(5, 35, 5)
# checkpoint_interval <- 150  # Set your checkpoint interval here
# checkpoint_counter <- 0
# 
# final_df <- foreach(nc_val = nc_seq, .combine = rbind, .packages = c("MCMCpack", "LaplacesDemon", "invgamma", "tidyverse", "HDInterval", "mvtnorm")) %:%
#   foreach(i = delta1, .combine = rbind) %dopar% {
#     cat("\n==========\nProcessing nc:", nc_val, " delta1:", i, "\n==========\n")
#     temp_results <- foreach(j = delta2, .combine = rbind, .packages = c("MCMCpack", "LaplacesDemon", "invgamma", "tidyverse", "HDInterval", "mvtnorm")) %dopar% {
#       set.seed(42)
#       temp_uc <- rep(2.6, pc) + c(rnorm(pc - 1, i, 0.05), 0)
#       temp_uh <-  2.6 + j
#       temp_Xc <- gen_X(nc, pc, nc_val)
#       res1 <- run_simulation(sig, 0.1, temp_uh, nh, temp_uc, temp_Xc, Zc, 5000, 10, 0.95)
#       temp_df <- data.frame(nc = nc_val, delta1 = i, delta2 = j, pow = res1$power, ess = res1$EHSS)
#       
#       # Return temp_df for combining
#       temp_df
#     }
#     
#     # Save checkpoint within the outer loop
#     checkpoint_counter <<- checkpoint_counter + nrow(temp_results)
#     if (checkpoint_counter >= checkpoint_interval) {
#       write.csv(temp_results, file = paste0("results/elastic_checkpoint_nc_", Sys.time(), ".csv"), row.names = FALSE)
#       checkpoint_counter <<- 0  # Reset the checkpoint counter
#     }
#     
#     temp_results  # Return the results for this iteration
#   }
# 
# write.csv(final_df, "results/elastic_results_nc_fc.csv", row.names = FALSE)
# 
# 
# stopCluster(cl)




final_df <- NULL
delta1 <- seq(0, 0.2, 0.02)
delta2 <- seq(0, 0.2, 0.02)
nc_seq <- seq(5, 35, 5)
for (nc_val in nc_seq){
  for (i in delta1){
    cat("\n==========\nProcessing nc:", nc_val, " delta1:", i, "\n==========\n")
    for (j in delta2){
      set.seed(42)
      temp_uc <- rep(2.6, pc) + c(rnorm(pc - 1, i, 0.05), 0)
      temp_uh <-  2.6 + j
      temp_Xc <- gen_X(nc, pc, nc_val)
      res1 <- run_simulation(sig, 0.1, temp_uh, nh, temp_uc, temp_Xc, Zc, 5000, 10, 0.95)
      temp_df <- data.frame(nc = nc_val, delta1 = i, delta2 = j, pow = res1$power,
                            ess = res1$EHSS)
      final_df <- rbind(final_df, temp_df)
      
      # Checkpointing
      if (nrow(final_df) %% 110 == 0) {
        write.csv(final_df, file = paste0("results/elastic_checkpoint_fc_nc_", nrow(final_df), ".csv"))
      }
    }
  }
}

write.csv(final_df, "results/elastic_results_nc_fc_temp.csv")


