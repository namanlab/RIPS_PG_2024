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
source("helper_oc.R")
source("helper_fc.R")

render_p1 <- function(data, var) {
  df_cur <- data
  
  y_var <- switch(var, "Power" = "pow", "PESS" = "ess", "MSE" = "mse")
  if (var == "PESS"){y_label = "ESS"}
  else if (var == "MSE"){y_label = "MSE"}
  else {y_label = "Power"}
  p_cur <- ggplot(df_cur, aes(x = tau, y = !!sym(y_var))) +
    labs(x = TeX("$\\tau$"), y = y_label) +
    geom_line() + geom_point() +
    theme_minimal() +
    theme(axis.title = element_text(size = 13))
  
  ggplotly(p_cur, height = 400)
}

render_p2 <- function(data, feat) {
  p_cur <- data %>% 
    filter(p == feat, tau %in% c(0.1, 0.2, 0.3, 0.4)) %>%
    ggplot(aes(x = val, fill = type)) + geom_density(alpha = 0.5) +
    facet_wrap(~tau, nrow = 3, scales = "free_x") + theme_bw() +
    labs(x = TeX("$\\theta_p$"), fill = "Type") +
    theme(axis.text.x = element_text(size = 4))
  
  ggplotly(p_cur)
}

render_p3 <- function(data, tau_val) {
  p_cur <- data %>% 
    filter(tau == tau_val) %>%
    ggplot(aes(x = val, fill = type)) + geom_density(alpha = 0.5) +
    facet_wrap(~p, nrow = 1, scales = "free_x") + theme_bw() +
    labs(x = TeX("$\\theta_p$"), fill = "Type") +
    theme(axis.text.x = element_text(size = 7))
  
  ggplotly(p_cur)
}


render_tile_plot <- function(data, var, title, nc_val) {
  df_cur <- data %>% filter(nc == nc_val)
  
  fill_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess", "Power"="pow")
  fill_label <- ifelse(var == "PESS", "ESS", "Power")
  if (var == "PESS"){color_scale <- scale_fill_gradient(low = "lightblue", high = "darkblue")}
  else {color_scale <- scale_fill_gradient(low = "red", high = "green", limits = c(0, 1))}
  p_cur <- ggplot(df_cur, aes(x = delta1, y = delta2, fill = !!sym(fill_var))) +
    labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
         y = TeX("$\\delta_2 = \\mu_h - \\mu_c$"),
         fill = fill_label) +
    color_scale +
    geom_tile() +
    theme_minimal() +
    theme(axis.title = element_text(size = 13),
          legend.position = "bottom") +
    xlim(0, 0.2) + ylim(0, 0.2)
  
  ggplotly(p_cur)
}

render_line_chart <- function(data, var, delta1_val, delta2_val) {
  y_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess", "Power"="pow")
  if (var == "PESS"){y_label <- "ESS"}
  else {y_label <- "Power"}
  p_cur <- ggplot(data %>% filter(delta1 == delta1_val, delta2 == delta2_val), aes(x = nc, y = !!sym(y_var))) +
    labs(x = "nc", y = y_label) +
    geom_line() + geom_point() +
    theme_minimal() +
    theme(axis.title = element_text(size = 13))
  
  ggplotly(p_cur)
}

render_scatter3d <- function(data, var) {
  fill_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess", "Power"="pow")
  fig <- plot_ly(data, x = ~delta1, y = ~delta2, z = as.formula(paste("~", fill_var)), 
                 color = ~nc, colors = viridis::viridis(length(unique(data$nc))),
                 type = 'scatter3d', mode = 'markers')
  fig <- fig %>%
    layout(scene = list(
      xaxis = list(title = TeX("delta_1")),
      yaxis = list(title = TeX("delta_2")),
      zaxis = list(title = 'Power')
    ))
  fig
}


# Helper function to read data from the uploaded file
read_data_OC <- function(file) {
  data <- read_excel(file)
  yc <- data[[1]]
  yt <- data[[2]]
  yh <- as.vector(as.matrix(data[, 3:ncol(data)]))
  list(yc = yc, yt = yt, yh = yh)
}

# Helper function to compute all Bayesian probabilities
compute_bayesian_probs_OC <- function(yc, yh, yt, N) {
  list(
    elastic_power = compute_oc_elastic_power(yc, yh, yt, N),
    elastic = compute_oc_elastic(yc, yh, yt, N),
    normalized = compute_oc_normalized(yc, yh, yt, N),
    commensurate = compute_oc_commensurate(yc, yh, yt, N),
    robust_map = compute_oc_robust_map(yc, yh, yt, N)
  )
}

# Helper function to generate result UI
generate_result_ui <- function(method_name, bayesian_prob, frequentist_p_value) {
  withMathJax(HTML(paste0(
    "\n\n\n",
    "<p><strong>", method_name, ":</strong></p>",
    "<p>Hypothesis: \\(\\mu_c = \\mu_t \\) vs \\(\\mu_c \\neq \\mu_t\\)</p>",
    "<p>Bayesian Probability of Rejection: ", format(bayesian_prob[1], digits = 4), "</p>",
    "<p>Prior Effective Sample Size: ", format(bayesian_prob[2], digits = 4), "</p>",
    "<p>Result: ", ifelse(bayesian_prob[1] >= 1 - 0.05, "Reject", "Fail to Reject"), "</p>",
    "\n",
    "<p><strong>Frequentist Approach:</strong></p>",
    "<p>Hypothesis: \\(\\mu_c = \\mu_t \\) vs \\(\\mu_c \\neq \\mu_t\\) </p>",
    "<p>Frequentist P-value: ", format(frequentist_p_value, digits = 4), "</p>",
    "<p>Result: ", ifelse(frequentist_p_value < 0.05, "Reject", "Fail to Reject"), "</p>"
  )))
}






# Helper function to read data from the uploaded file
read_data_FC <- function(file1, file2) {
  
  data <- read_excel(file1)
  data[[1]] = as.numeric(as.factor(data[[1]]))
  data[[2]] = as.numeric(as.factor(data[[2]]))
  data[[3]] = as.numeric(as.factor(data[[3]]))
  colnames(data) <- c("ids", "sides", "treat_grp", "y")
  ids <- data[[1]]
  sides <- data[[2]]
  treat_grp <- data[[3]]
  nc = length(unique(ids))
  pc = length(unique(treat_grp))
  Zc = gen_Z(nc)
  Xc <- matrix(0, nrow = 2 * nc, ncol = pc)
  xc <- rep(0, 2*nc)
  for (i in 1:n) {
    temp_df <- data %>% filter(ids == i) %>% arrange(sides)
    cur_t <- temp_df %>% pull(treat_grp)
    cur_y <- temp_df %>% pull(y)
    Xc[2 * i - 1, cur_t[1]] <- 1
    Xc[2 * i, cur_t[2]] <- 1
    xc[2*i - 1] <- cur_y[1]
    xc[2*i] <- cur_y[2]
  }
  
  data_hist <- read_excel(file2)
  xh <- as.vector(as.matrix(data[, 1:ncol(data_hist)]))
  
  list(Xc = Xc, Zc = Zc, xc = xc, xh = xh)
}

# Helper function to compute all Bayesian probabilities
compute_bayesian_probs_FC <- function(Xc, Zc, xc, xh, N) {
  list(
    elastic = compute_fc_elastic(Xc, Zc, xc, xh, N)
  )
}

generate_result_ui_FC <- function(method_name, bayesian_probs) {
  # Generate UI elements for each Bayesian probability result
  lapply(seq_along(bayesian_probs), function(i) {
    withMathJax(HTML(paste0(
      "\n\n\n",
      "<p><strong>", method_name, " (Treatment ", i, "):</strong></p>",
      "<p>Hypothesis: \\(\\mu_c = \\mu_t \\) vs \\(\\mu_c \\neq \\mu_t\\)</p>",
      "<p>Bayesian Probability of Rejection: ", format(bayesian_probs[[i]], digits = 4), "</p>",
      "<p>Result: ", ifelse(bayesian_probs[[i]][1] >= 1 - 0.05, "Reject", "Fail to Reject"), "</p>"
    )))
  })
}




