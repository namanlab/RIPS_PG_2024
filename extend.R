library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)
library(RBesT)

pm = c(0, 1, 1, 1) # (mu0, n0, a0, b0)

run_simulation <- function(nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, 
                           H = 1, N, R, cutoff){
  
  set.seed(42)
  xh_list <- lapply(1:length(nh), function(i) rnorm(nh[i], uh[i], sigh[i])) # list of historical data
  
  # Metrics to calculate
  rej_null <- NULL # prob of rejections
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
  params <- get_params(xh_list, nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H, N, R)
  
  for(trial in 1:R){
    
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rnorm(nc, uc, sigc)
    xt <- rnorm(nt, ut, sigt)
    
    # posterior for control arm -- TO CHANGE DEPENDING ON METHOD
    res_inter <- get_control_prost(params, xc, xh_list, nc, nh, H, N)
    muc <- res_inter$prost_samples
    
    # posterior for treatment arm
    t.par1 <- (sum(xt) + pm[2]*pm[1])/(length(xt) + pm[2])
    b1 <- 1/2*(2*pm[4] + (length(xt) - 1)*var(xt) +
                 length(xt)*pm[2]*(mean(xt) - pm[1])^2/(length(xt) + pm[2]))
    a1 <- pm[3] + length(xt)/2
    t.par2 <- b1/(a1*(length(xt) + pm[2]))
    mut <- rst(N, t.par1, sqrt(t.par2), length(xt)/2 + pm[3])
    
    # METRICS:
    # probability that treatment is superior to control
    pp <- mean(mut <= muc)
    # number of rejections of null
    rej_null <- c(rej_null, pp)
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
      b1 <- 1/2*(2*pm[4] + (sum(sapply(xh_list, function(xh) (length(xh) - 1) * var(xh))) +
                              sum(sapply(xh_list, function(xh, i) length(xh) * pm[2] * (mean(xh) - pm[1])^2 / (length(xh) + pm[2])))))
      a1 <- pm[3] + sum(sapply(xh_list, length))/2
      sig <- rinvgamma(N, a1, b1)
      mu <- (nc*mean(xc) + sum(sapply(1:length(xh_list), function(i) nh[i] * mean(xh_list[[i]]))) + pm[2]*pm[1]) / (nc + sum(nh) + pm[2])
      var <- sig/(1/nc + 1/sum(1/nh) + 1/pm[2])
      muc_def <- rnorm(N, mu, sqrt(var)) 
      distr_plot_prost <- c(distr_plot_prost, muc_def)
      distr_plot_prost_method <- c(distr_plot_prost_method, muc)
    }
  }
  timeend <- Sys.time()
  EHSS <- mean(ess)
  prob_rej <- mean(rej_null)
  quantile_interval_count_mean <-  quantile_interval_count/R
  bias_point_est <- mean(point_est) - uc
  var_point_est <- var(point_est)
  mse_point_est <- bias_point_est^2 + var_point_est
  width_quantile_interval_mean <- mean(width_quantile_interval)
  cat("probability of claiming efficacy is", prob_rej, "\n")
  cat("effective historical sample size is", formatC(EHSS, digits = 2, format = "f"), sep = " ", "\n")
  # cat("Mean Width of Credible Interval for Control Prior", formatC(width_quantile_interval_mean, digits = 4, format = "f"), sep = " ", "\n")
  # cat("% of times muc is in quantile interval is", formatC(quantile_interval_count_mean*100, digits = 4, format = "f"), sep = " ", "\n")
  # cat("Bias of point estiamtor based on control prior", formatC(bias_point_est, digits = 4, format = "f"), sep = " ", "\n")
  # cat("Variance value of point estiamtor based on control prior", formatC(var_point_est, digits = 4, format = "f"), sep = " ", "\n")
  cat("MSE of point estiamtor based on control prior", formatC(mse_point_est, digits = 4, format = "f"), sep = " ", "\n")
  # cat("total time for", R, "simulations is", formatC(timeend - timestart, digits = 4, format = "f"), sep = " ", "\n")
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


# can return a list of all params if needed to be used by model, return null if not needed
get_params <- function(xh_list, nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H, N, R){
  params <- lapply(1:length(xh_list), function(i) decide_para(c=1, xh_list[[i]], nh[i], nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 50000))
  return(params)
}

# returns a list of:
# i) prost_samples: N samples from posterior distribution
# ii) ess: effective sample size for the prior used
get_control_prost <- function(params, xc, xh_list, nc, nh, H, N){
  xh <- do.call(c, xh_list)
  b <- 1/2*(2*pm[4] + (nc-1)*var(xc) + sum(sapply(xh_list, function(xh) (length(xh) - 1) * var(xh))) +
              sum(sapply(1:length(xh_list), function(i) length(xh_list[[i]]) * pm[2] * (mean(xh_list[[i]]) - pm[1])^2 / (length(xh_list[[i]]) + pm[2]))))
  a <- pm[3] + (nc + sum(nh))/2
  sig <- rinvgamma(N, a, b)
  mu <- (nc*mean(xc) + sum(sapply(1:length(xh_list), function(i) nh[i] * mean(xh_list[[i]]))) + pm[2]*pm[1]) / (nc + sum(nh) + pm[2])
  var <- sig/(1/nc + 1/sum(1/nh) + 1/pm[2])
  muc <- rnorm(N, mu, sqrt(var)) 
  ess <- sum(sapply(params, function(x) x$sample_size))
  return(list(prost_samples = muc, ess = ess))
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
      t[i,j] <- abs(u0-mean(xc))/(sqrt(sp/n0 + sp/nc)) # *max(n0, nc)^(-1/4)
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


## settings
nc <- 30 # current control size
nt <- 29 # current treatment size
nh <- c(20, 25, 29, 24) # historical control size
sigc <- 0.153 # control sd
sigt <- 0.17 # treatment sd
sigh <- c(0.09, 0.09, 0.33, 0.22) # historical sd
uc <- 1.26 + 1.33 # true mean of control
uh <- c(1.24 + 1.62, 1.21 + 1.2, 1.05 + 1.73, 1.18 + 1.45) # true mean of historical
ut <- 1.08 + 1.33 # true mean of treatment

# strong congruence between control and historical
for (i in seq(nc, 2)){
  print(i)
  res1 <- run_simulation(nt, i, nh, sigc, sigt, sigh, uc, ut, uh, H = 1, N = 10000, R = 100, cutoff = 0.95) 
}

# res1$plot_comp
# res1$plot_density
