library(MCMCpack)
library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)

run_simulation <- function(nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H = 1, N, R, cutoff){
  
  set.seed(42)
  xh <- rnorm(nh, uh, sigh) # historical data
  
  # Metrics to calculate
  rej_null <- 0 # number of rejections
  ess <- NULL # effective sample size
  timestart <- Sys.time() # computation time
  quantile_interval_count <- 0 # number of times muc is in quantile interval
  width_quantile_interval <- NULL # width of credible interval
  point_est <- NULL # point estimator based on control prior
  xc_act <- NULL # actual control samples
  xc_samp <- NULL # samples based on prosterior predictive
  
  # for plot of distributions (only calculated in last sim)
  distr_plot_prost <- NULL
  distr_plot_prost_method <- NULL
  
  # any params for models to be computed a-priori -- TO CHANGE DEPENIDNG ON METHOD
  params <- get_params(xh, nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H, N, R)
  
  for(trial in 1:R){
    
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rnorm(nc, uc, sigc)
    xt <- rnorm(nt, ut, sigt)
    
    # posterior for control arm -- TO CHANGE DEPENIDNG ON METHOD
    res_inter <- get_control_prost(params, xc, xh, nc, nh, H, N)
    muc <- res_inter$prost_samples
    
    # posterior for treatment arm
    t.par1 <- mean(xt)
    t.par2 <- var(xt)/nt
    mut <- rst(N, t.par1, sqrt(t.par2), nt - 1)
    
    # METRICS:
    # probability that treatment is superior to control
    pp <- mean(mut > muc)
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

# Function to perform MCMC sampling for the posterior using power prior
get_control_prost <- function(params, xc, xh, nc, nh, H, N) {
  all_c <- c(xc, xh)
  t.par1 <- mean(all_c)
  t.par2 <- var(all_c)/length(all_c)
  mu_samples <- rst(N, t.par1, sqrt(t.par2), length(all_c) - 1)
  ess <- var(c(xc, xh))/var(mu_samples)
  list(prost_samples = mu_samples, ess = ess)
}

#########--------------------------------------------------------------#########
############################ NO CHANGE TO BE MADE  #############################
#########--------------------------------------------------------------#########




## settings
nc <- 25 # current control size
nt <- 50 # current treatment size
nh <- 50 # historical control size
sigc <- 1 # control sd
sigt <- 1 # treatment sd
sigh <- 1 # historical sd
uc <- 1 # true mean of control

# strong congruence between control and historical
res1 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1, uh = 1, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res1$plot_comp
res1$plot_density
res2 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1.5, uh = 1, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res2$plot_comp
res2$plot_density

# weak congruence between control and historical
res3 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1, uh = 1.2, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res3$plot_comp
res3$plot_density
res4 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1.5, uh = 1.2, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res4$plot_comp
res4$plot_density

# no congruence between control and historical
res5 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1, uh = 1.5, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res5$plot_comp
res5$plot_density
res6 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1.5, uh = 1.5, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res6$plot_comp
res6$plot_density

# Combine results into a list and Save the list as an RDS file
results <- list(res1 = res1, res2 = res2, res3 = res3, res4 = res4, res5 = res5, res6 = res6)
saveRDS(results, file = "../results/simulation_results_MAP_prior_normal.rds")

# Load the results: results <- readRDS("simulation_results.rds")



# Delta vs MSE Plot
# delta <- seq(0, 0.5, by = 0.01)
# mse_vals <- NULL
# bias_vals <- NULL
# var_vals <- NULL
# for (i in delta){
#   sim <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1, uh = 1 + i, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
#   mse_vals <- c(mse_vals, sim$mse_point_est)
#   bias_vals <- c(bias_vals, sim$bias_point_est)
#   var_vals <- c(var_vals, sim$var_point_est)
# }
# ggplot(data = data.frame(delta = delta, mse = mse_vals)) +
#   geom_line(aes(x = delta, y = mse_vals)) +
#   theme_bw()
