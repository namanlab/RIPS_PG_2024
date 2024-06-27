library(MCMCpack)
library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)

run_simulation <- function(nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H = 1, N, R, cutoff){
  
  set.seed(42)
  xh_list <- lapply(1:length(nh), function(i) rnorm(nh[i], uh[i], sigh[i])) # list of historical data
  xh <- do.call(c, xh_list)
  nh <- length(xh)
  sigh <- sd(xh)
  uh <- mean(xh)
  
  # Metrics to calculate
  rej_null <- 0 # number of rejections
  ess <- NULL # effective sample size
  timestart <- Sys.time() # computation time
  quantile_interval_count <- 0 # number of times muc is in quantile interval
  width_quantile_interval <- NULL # width of credible interval
  point_est <- NULL # point estimator based on control prior
  xc_act <- NULL # actual control samples
  xc_samp <- NULL # samples based on posterior predictive
  
  # for plot of distributions (only calculated in last sim)
  distr_plot_prost <- NULL
  distr_plot_prost_method <- NULL
  
  # any params for models to be computed a-priori -- TO CHANGE DEPENDING ON METHOD
  params <- get_params(xh, nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H, N, R)
  
  for(trial in 1:R){
    
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rnorm(nc, uc, sigc)
    xt <- rnorm(nt, ut, sigt)
    
    # posterior for control arm -- TO CHANGE DEPENDING ON METHOD
    res_inter <- get_control_prost(params, xc, xh, nc, nh, H, N)
    muc <- res_inter$prost_samples
    
    # posterior for treatment arm
    t.par1 <- mean(xt)
    t.par2 <- var(xt)/nt
    mut <- rst(N, t.par1, sqrt(t.par2), nt - 1)
    
    # METRICS:
    # probability that treatment is superior to control
    pp <- mean(mut <= muc)
    # number of rejections of null
    if(pp >= cutoff){
      rej_null <- rej_null + 1
    }
    # effective sample size:
    ess <- c(ess, res_inter$ess)
    # Count if the true value is within the quantile interval
    quantile_interval <- quantile(muc, probs = c(0.025, 0.975))
    quantile_interval_count <- quantile_interval_count + ifelse((quantile_interval[1] <= uc) && (quantile_interval[2] >= uc), 1, 0)
    # Width of the quantile interval
    width_quantile_interval <- c(width_quantile_interval, quantile_interval[2] - quantile_interval[1])
    # point estimator
    val_point_est <- mean(muc)
    point_est <- c(point_est, val_point_est)
    # posterior predictive check:
    xc_act <- c(xc_act, rnorm(length(muc), uc, sigc))
    xc_samp <- c(xc_samp, rnorm(length(muc), muc, sigc))
    
    # plot of prost distr:
    if (trial == R){
      # for plot of distributions (only calculated in last sim)
      all_c <- c(xc, xh)
      t.par1 <- mean(all_c)
      t.par2 <- var(all_c)/length(all_c)
      muc_def <- rst(N, t.par1, sqrt(t.par2), length(all_c) - 1)
      distr_plot_prost <- c(distr_plot_prost, muc_def)
      distr_plot_prost_method <- c(distr_plot_prost_method, muc)
    }
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
  cat("% of times muc is in quantile interval is", formatC(quantile_interval_count_mean*100, digits = 4, format = "f"), sep = " ", "\n")
  cat("Bias of point estiamtor based on control prior", formatC(bias_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("Variance value of point estiamtor based on control prior", formatC(var_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("MSE of point estiamtor based on control prior", formatC(mse_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("total time for", R, "simulations is", formatC(timeend - timestart, digits = 4, format = "f"), sep = " ", "\n")
  plot_comp <- ggplot(data = tibble("Control Distribution" = xc_act, "Prost Predictive" = xc_samp) %>% 
                        pivot_longer(1:2, names_to = "type", values_to = "val") ) +
    geom_density(aes(x = val, fill = type), alpha = 0.3) +
    theme_bw() + scale_fill_manual(values = c("yellow", "blue"))
  plot_density <- ggplot(data = tibble("All Historical Data" = distr_plot_prost, "Adjusting with Method" = distr_plot_prost_method) %>% 
                           pivot_longer(1:2, names_to = "type", values_to = "val") ) +
    geom_density(aes(x = val, fill = type), alpha = 0.3) +
    theme_bw() + scale_fill_manual(values = c("yellow", "blue")) +
    geom_vline(xintercept = uc, linetype = "dashed")
  return(list(prob_rej = prob_rej, EHSS = EHSS, 
              width_quantile_interval_mean = width_quantile_interval_mean, 
              quantile_interval_count_mean = quantile_interval_count_mean,
              bias_point_est = bias_point_est, var_point_est = var_point_est, 
              mse_point_est = mse_point_est,
              time_diff = timeend - timestart, 
              plot_comp = plot_comp, plot_density = plot_density))
  
}


get_params <- function(xh, nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H, N, R) {
  list()
}

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
  
  init_values <- c(mean(xc), log(sd(xc)), logit(0.5), mean(xh), log(1))  # Initial values for MCMC sampling
  # Perform MCMC sampling for posterior
  samples <- MCMCmetrop1R(log_posterior, theta.init = init_values, 
                          mcmc = N, burnin = 1000, thin = 1, force.samp = T)
  
  # Extract posterior samples for mu, sigma, a0, theta0, and tau
  mu_samples <- samples[, 1]
  sigma_samples <- exp(samples[, 2])
  a0_samples <- inv_logit(samples[, 3])
  theta0_samples <- samples[, 4]
  tau_samples <- exp(samples[, 5])
  
  # Calculate effective sample size (ESS)
  ess <- mean(sigma_samples)/var(mu_samples)
  
  list(prost_samples = mu_samples, ess = ess)
}


#########--------------------------------------------------------------#########
#################################### MODEL #####################################
#########--------------------------------------------------------------#########




## settings
nc <- 30 # current control size
nt <- 29 # current treatment size
nh <- c(20, 25, 29, 24) # historical control size
sigc <- 0.153 # control sd
sigt <- 0.17 # treatment sd
sigh <- c(0.09, 0.09, 0.33, 0.22) # historical sd
uc <- 1.26 + 1.33 # true mean of control

final_df <- NULL
delta1 <- seq(-1, 1, 0.1)
delta2 <- seq(-1, 1, 0.1)
for (i in delta1){
  print(i)
  for (j in delta2){
    ut <- 1.08 + 1.33 + i
    set.seed(42)
    uh <-  c(1.24 + 1.62, 1.21 + 1.2, 1.05 + 1.73, 1.18 + 1.45) + rnorm(4, j, 0.05)
    res1 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H = 1, N = 10000, R = 100, cutoff = 0.95) 
    temp_df <- data.frame(delta1 = i, delta2 = j, pow = res1$prob_rej, ess = res1$EHSS)
    final_df <- rbind(final_df, temp_df)
  }
}




write.csv(final_df, "results/commensurate_results.csv")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

final_df <- read.csv("results/commensurate_results.csv")
library(plotly)

# ESS
delta1 <- seq(-1, 1, length.out = 21)
delta2 <- seq(-1, 1, length.out = 21)
ess <- matrix(final_df$ess, nrow = 21, ncol = 21, byrow = TRUE)
plot_ly(
  x = ~delta2, y = ~delta1, z = ~ess,
  type = 'surface'
) %>% layout(
  scene = list(
    xaxis = list(title = "Delta 2"),
    yaxis = list(title = "Delta 1"),
    zaxis = list(title = "ESS")
  ),
  title = "3D Plot of ESS vs Delta1 and Delta2"
)
plot_ly(
  x = ~delta2, y = ~delta1, z = ~ess,
  type = 'heatmap'
) %>% layout(
  scene = list(
    xaxis = list(title = "Delta 2"),
    yaxis = list(title = "Delta 1"),
    zaxis = list(title = "ESS")
  ),
  title = "3D Plot of ESS vs Delta1 and Delta2"
)


# Power
pow <- matrix(final_df$pow, nrow = 21, ncol = 21, byrow = TRUE)
plot_ly(
  x = ~delta2, y = ~delta1, z = ~pow,
  type = 'surface'
)  %>% layout(
  scene = list(
    xaxis = list(title = "Delta 2"),
    yaxis = list(title = "Delta 1"),
    zaxis = list(title = "Power")
  ),
  title = "3D Plot of Power vs Delta1 and Delta2"
)
plot_ly(
  x = ~delta2, y = ~delta1, z = ~pow,
  type = 'heatmap'
) %>% layout(
  scene = list(
    xaxis = list(title = "Delta 2"),
    yaxis = list(title = "Delta 1"),
    zaxis = list(title = "Power")
  ),
  title = "3D Plot of Power vs Delta1 and Delta2"
)





