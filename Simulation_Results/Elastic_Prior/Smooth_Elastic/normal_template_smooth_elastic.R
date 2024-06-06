#################
################# EP1 with smooth elastic function #######################
#################

library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)
library(RBesT)

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
  return(list(prob.rej = prob_rej, width_quantile_interval_mean = width_quantile_interval_mean, 
              quantile_interval_count_mean = quantile_interval_count_mean,
              bias_point_est = bias_point_est, var_point_est = var_point_est, mse_point_est = mse_point_est,
              time_diff = timeend - timestart, plot_comp = plot_comp, distr_plot_prost = distr_plot_prost,
              distr_plot_prost_method = distr_plot_prost_method, plot_density = plot_density))
  
}


# can return a list of all params if needed to be used by model, return null if not needed
get_params <- function(xh, nt, nc, nh, sigc, sigt, sigh, uc, ut, uh, H, N, R){
  params <- decide_para(c=1, xh, nh, nc, gamma=1, q1=0.95, q2=0.02, small = 0.01, large = 0.99, R = 50000)
  return(params)
}

# returns a list of:
# i) prost_samples: N samples from prosterior distribution
# ii) ess: effective sample size for the prior used
get_control_prost <- function(params, xc, xh, nc, nh, H, N){
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
  res <- list(ess = gt*nh, prost_samples = temp$muc_post) # 
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
# b: tuning parameter in elastic function
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

#-----Function to sample posterior of mean value for control arm------------#
# Inputs:
# x0: historical data
# n0: sample size for historical data
# xc: control data
# nc: sample size for control arm
# gt: value of elastic function at current congruence measure t between historical and control data
# sim: number of simulated trial
# nburn: number of simulated trial in burn-in process
# 
# Outputs:
# muc_post: posterior of mean value for control arm
sample_poster <- function(x0, n0, xc, nc, gt, sim=20000, nburn=10000){
  muc_post <- numeric(sim-nburn)
  muc_ini <- mean(xc)
  muc <- muc_ini
  for (s in 1:sim) {
    # sample control variance from posterior
    alpha <- (1 + nc)/2
    beta <- nc*(var(xc) + (mean(xc) - muc)^2)/2
    sig <- rinvgamma(1, alpha, beta)
    # sample mean from posterior
    D <- var(x0)/(n0*gt)
    mu <- (nc*mean(xc)*D + sig*mean(x0))/(nc*D + sig)
    var <- sig*D/(D*nc + sig)
    muc <- rnorm(1, mu, sqrt(var)) 
    if(s > nburn){
      muc_post[s-nburn] <- muc
    }
  }
  return(list(muc_post=muc_post))
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
res15plot_density
res6 <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1.5, uh = 1.5, H = 1, N = 10000, R = 100, cutoff = 0.95) # false null
res6$plot_comp
res6$plot_density

# Combine results into a list and Save the list as an RDS file
results <- list(res1 = res1, res2 = res2, res3 = res3, res4 = res4, res5 = res5, res6 = res6)
saveRDS(results, file = "../../results/ssimulation_results_elastic_prior_normal.rds")

# Load the results: results <- readRDS("simulation_results.rds")



# Delta vs MSE Plot
delta <- seq(0, 0.5, by = 0.05)
mse_vals <- NULL
for (i in delta){
  sim <- run_simulation(nt, nc, nh, sigc, sigt, sigh, uc, ut = 1, uh = 1 + i, H = 1, N = 10000, R = 100, cutoff = 0.95) # true null
  mse_vals <- c(mse_vals, sim$mse_point_est)
}
ggplot(data = data.frame(delta = delta, mse = mse_vals)) +
  geom_line(aes(x = delta, y = mse_vals)) +
  theme_bw() 