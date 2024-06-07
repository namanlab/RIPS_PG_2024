library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)

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
  
  # for plot of distributions (only calculated in last sim)
  distr_plot_prost <- NULL
  distr_plot_prost_method <- NULL
  
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
    
    # plot of prost distr:
    if (trial == R){
      # for plot of distributions (only calculated in last sim)
      pc_def <-rbeta(N, params$c_alpha + xc + xh, params$c_beta + nc + nh - xc - xh)
      distr_plot_prost <- c(distr_plot_prost, pc_def)
      distr_plot_prost_method <- c(distr_plot_prost_method, pc_post)
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
  plot_density <- ggplot(data = tibble("All Historical Data" = distr_plot_prost, "Adjusting with Method" = distr_plot_prost_method) %>% 
                           pivot_longer(1:2, names_to = "type", values_to = "val") ) +
    geom_density(aes(x = val, fill = type), alpha = 0.3) +
    theme_bw() + scale_fill_manual(values = c("yellow", "blue")) +
    geom_vline(xintercept = pc, linetype = "dashed")
  return(list(prob.rej = prob_rej, width_quantile_interval_mean = width_quantile_interval_mean, 
              quantile_interval_count_mean = quantile_interval_count_mean,
              bias_point_est = bias_point_est, var_point_est = var_point_est, mse_point_est = mse_point_est,
              time_diff = timeend - timestart, plot_comp = plot_comp, plot_density = plot_density))
  
}


# can return a list of all params if needed to be used by model, return null if not needed
get_params <- function(xh, nt, nc, nh, pc, pt, ph, H, N, R){
  params <- decide_para(c=1, xh, nh, nc, gamma=0.24, q1=0.9, q2=0.1, small = 0.01, large = 0.99, R = 50000)
  params$t_alpha = 0.1
  params$t_beta = 0.1
  params$c_alpha = 0.1
  params$c_beta = 0.1
  return(params)
}

# returns a list of:
# i) prost_samples: N samples from prosterior distribution
# ii) ess: effective sample size for the prior used
get_control_prost <- function(params, xc, xh, nc, nh, H, N){
  a <- params$a 
  b <- params$b 
  c <- params$c
  c_alpha <- params$c_alpha
  c_beta <- params$c_beta
  # calculate g(t)
  phat <- (xh + xc)/(nh + nc)
  O <- cbind(c(xh, xc), c(nh-xh, nc-xc))
  E <- cbind(c(nh, nc) * phat, c(nh, nc) * (1 - phat))
  cong_measure <- max(nh, nc)^(-1/4)*sum((O - E)^2/E)
  gt <- 1/(1+exp(a + b*(log(cong_measure))^c))
  if(gt < 0.01) gt <- 0.01 # numerical stability
  # posterior for control arm
  c_pa1 <- (c_alpha + xh)*gt + xc
  c_pa2 <- (c_beta + nh - xh)*gt + nc - xc
  pc_post <- rbeta(N, c_pa1, c_pa2)
  res <- list(ess = gt*nh, prost_samples = pc_post)
  return(res)
}



#-----Function to decide tuning parameter a and b in logistic elastic function------------#
# Inputs:
# c: pre-specified tuning parameter controlling the shape of elastic function. Here, we fix it by 1
# x0: historical data
# n0: sample size for historical data
# nc: sample size for control arm
# gamma: clinically highly meaningful difference,  -- 
# q1: q1th percentile of congruence measure T in the homogeneous case, 0.95 --
# q2: q2th percentile of congruence measure T in the heterogeneous case, 0.02 --
# small: value of elastic function in the heterogeneous case, 0.01 --
# large: value of elastic function in the homogeneous case, 0.99 --
# R: the number of simulations
# 
# Outputs:
# a: tuning parameter in elastic function 
# b: tuning parameter in elastic function)
decide_para <- function(x0, n0, nc, gamma, q1, q2, c, small, large, R=50000){
  set.seed(10)
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
  return(list(a=a, b=b, c = c))
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
res1$plot_density
res2 <- run_simulation(nt, nc, nh, pc, pt = 0.6, ph = 0.4, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res2$plot_comp
res2$plot_density

# weak congruence between control and historical
res3 <- run_simulation(nt, nc, nh, pc, pt = 0.4, ph = 0.45, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res3$plot_comp
res3$plot_density
res4 <- run_simulation(nt, nc, nh, pc, pt = 0.6, ph = 0.45, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res4$plot_comp
res4$plot_density

# no congruence between control and historical
res5 <- run_simulation(nt, nc, nh, pc, pt = 0.4, ph = 0.6, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
res5$plot_comp
res5$plot_density
res6 <- run_simulation(nt, nc, nh, pc, pt = 0.6, ph = 0.6, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res6$plot_comp
res6$plot_density

# Combine results into a list and Save the list as an RDS file
results <- list(res1 = res1, res2 = res2, res3 = res3, res4 = res4, res5 = res5, res6 = res6)
saveRDS(results, file = "../../results/simulation_results_elastic_prior_bernoulli.rds")

# Load the results: results <- readRDS("simulation_results.rds")




# Delta vs MSE Plot
delta <- seq(0, 0.2, by = 0.02)
mse_vals <- NULL
for (i in delta){
  sim <- run_simulation(nt, nc, nh, pc, pt = 0.4, ph = 0.4 + i, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
  mse_vals <- c(mse_vals, sim$mse_point_est)
}
ggplot(data = data.frame(delta = delta, mse = mse_vals)) +
  geom_line(aes(x = delta, y = mse_vals)) +
  theme_bw() 
