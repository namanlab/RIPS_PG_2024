library(MCMCpack)
library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)
library(mvtnorm)

## settings
nc <- 99
pc <- 3
uc <- c(2.6, 2.6, 2.6)
uh = 2.6
Zc <- gen_Z(nc)
Xc <- gen_X(nc, pc, round(2*nc/pc))
sig = 0.2
tau = 0.1
nh = 50

set.seed(42)
xh <- rnorm(nh, uh, sqrt(sig^2 + tau^2))
pc <- ncol(Xc)
set.seed(100)
nc <- ncol(Zc)
overall_noise <- rnorm(2*nc, 0, sig^2)
indiv_noise <- rnorm(nc, 0, tau^2)
xc <- Xc %*% uc + Zc %*% indiv_noise + overall_noise

log_inverse_gamma <- function(sigma2) {
  if (sigma2 > 0) {
    return(dinvgamma(sigma2, shape = 2, scale = 2, log = TRUE))
  } else {
    return(-Inf)
  }
}

gt = 1
log_posterior <- function(theta) {
  mu <- theta[1:pc]
  sigma2 <- theta[pc + 1]
  tau2 <- theta[pc + 2]
  
  # Likelihood
  log_lik_current <- dmvnorm(t(xc), mean = Xc %*% mu, sigma = tau2*(Zc %*% t(Zc)) + sigma2*diag(2*nc), log = T)
  log_lik_hist <- dnorm(mu[pc], mean(xh), sqrt((sigma2 + tau2)/(nh*gt)), log = T)
  
  # Prior distributions
  log_prior_sigma <- log_inverse_gamma(sigma2)
  log_prior_tau <- log_inverse_gamma(tau2)
  
  # Posterior
  log_post <- log_lik_current + log_lik_hist + log_prior_sigma + log_prior_tau 
  
  if (is.na(log_post) || is.nan(log_post) || is.infinite(log_post)) {
    return(-Inf)  # Return -Inf for invalid log-posterior values
  }
  return(log_post)
}

init_mu <- rep(0, pc)
for (i in 1:pc){init_mu[i] <- mean(xc[which(Xc[,i] == 1)])}
init_values <- c(init_mu, 1, 1) 
exp(log_posterior(init_values))



f <- function(x){
  theta <- init_values
  theta[5] <- x
  exp(log_posterior(theta))
}
x_values <- seq(0, 1, length.out = 100)
y_values <- sapply(x_values, f)
data <- data.frame(x = x_values, y = y_values)
ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Plot of 1D Function",
       x = "x",
       y = "f(x)") +
  theme_minimal()




f <- function(theta){exp(log_posterior(theta))}

library(cubature)
# Define a wrapper function to integrate with respect to a specific parameter
integrate_with_respect_to <- function(f, position, lower_limit, upper_limit, other_params) {
  # Create an integrand function for cubature
  integrand <- function(x) {
    theta <- other_params
    theta[position] <- x
    f(theta)
  }
  
  # Perform the integration
  result <- adaptIntegrate(integrand, lowerLimit = lower_limit, upperLimit = upper_limit)
  result$value
}

theta <- rep(0, pc + 2)  # 5-dimensional vector with all elements set to 0
position <- pc + 2  # Integrate with respect to the 3rd parameter
lower_limit <- -Inf  # Define the lower limit of integration
upper_limit <- Inf  # Define the upper limit of integration

# Perform the integration
integrated_value <- integrate_with_respect_to(f, position, lower_limit, upper_limit, theta)
print(integrated_value)


library(rjags)

get_control_prost <- function(Xc, Zc, xc, xh, nc, nh, N, pc) {
  
  xcc <- xc[which(Xc[,pc] == 1)]
  ncc <- length(xcc)
  gt <- getGT(xcc, xh, ncc, nh, N, pc)
  
  # JAGS model as a string
  jags_model <- "
  model {
    for (i in 1:nc) {
      xc[i] ~ dnorm(mu[i], tau)
      mu[i] <- inprod(Xc[i,], beta)
    }
    
    for (j in 1:pc) {
      beta[j] ~ dnorm(0, 0.001)
    }
    
    tau <- pow(sigma, -2)
    sigma ~ dunif(0, 100)
    
    for (k in 1:nh) {
      xh[k] ~ dnorm(beta[pc], tau_hist)
    }
    
    tau_hist <- pow(tau2, -1)
    tau2 ~ dunif(0, 100)
  }
  "
  
  # Data list for JAGS
  data_jags <- list(
    xc = xc,
    Xc = as.matrix(Xc),
    nc = nc,
    pc = pc,
    xh = xh,
    nh = nh
  )
  
  # Initial values for the parameters
  init_values <- function() {
    list(beta = rep(0, pc), sigma = 1, tau2 = 1)
  }
  
  # Parameters to monitor
  params <- c("beta", "sigma", "tau2")
  
  # Run the JAGS model
  model <- jags.model(textConnection(jags_model), data = data_jags, inits = init_values, n.chains = 3)
  update(model, 1000)  # Burn-in
  
  samples <- coda.samples(model, variable.names = params, n.iter = N)
  samples_matrix <- as.matrix(samples)
  
  # Extract posterior samples for beta, sigma, and tau2
  mu_samples <- samples_matrix[, grep("beta", colnames(samples_matrix))]
  sigma_samples <- samples_matrix[, "sigma"]
  tau2_samples <- samples_matrix[, "tau2"]
  
  list(prost_samples = mu_samples, ess = gt * nh)
}


