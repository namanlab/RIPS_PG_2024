library(tidyverse)
library(e1071)
library(metR)
library(stringr)
library(latex2exp)

data <- read.csv("results/s.csv")

df_f <- expand.grid(delta1 = seq(0, 1, length.out = 100),
                    delta2 = seq(0, 1, length.out = 100))
for (k in 1:length(nc_seq)){
  cur_data <- data %>% filter(nc == nc_seq[k]) %>% na.omit()
  cur_data$class <- ifelse(cur_data$pow > 0.5, 1, 0)
  svm_model <- svm(class ~ delta1 + delta2, data = cur_data, kernel = "radial")
  df_f$pred <- predict(svm_model, grid)
  colnames(df_f)[k + 2] = str_c("nc_", nc_seq[k])
}

df_f %>% pivot_longer(3:ncol(df_f), names_to = "nc", values_to = "pred") %>%
  ggplot() +
  geom_contour(aes(x = delta1, y = delta2, z = as.numeric(pred),
                   color = nc), breaks = 0.5) +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$")) +
  theme_bw() +
  theme(axis.title = element_text(size = 13))  +
  scale_color_viridis_d() +
  xlim(c(0, 1)) + ylim(0, 1)

df_f %>% pivot_longer(3:ncol(df_f), names_to = "nc", values_to = "pred") %>%
  ggplot(aes(x = delta1, y = delta2, z = pred)) +
  geom_contour_fill() +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$"),
       fill = "nc") +
  theme_bw() +
  theme(axis.title = element_text(size = 13)) +
  xlim(0, 1) + ylim(0, 1) +
  facet_grid(~nc)



data %>% 
  ggplot(aes(x = delta1, y = delta2, z = pow1)) +
  geom_contour_fill() +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$"),
       fill = "nc") +
  theme_bw() +
  theme(axis.title = element_text(size = 13)) +
  xlim(0, 0.2) + ylim(0, 0.2) +
  facet_grid(~nc)

# ------------- contour plot

df_f <- expand.grid(delta1 = seq(0, 1, length.out = 100),
                    delta2 = seq(0, 1, length.out = 100))
cur_data <- data %>% filter(nc == 5) 
cur_data$class <- ifelse(cur_data$pow1 > 0.95, 1, 0)
svm_model <- svm(class ~ delta1 + delta2, data = cur_data, kernel = "radial")
df_f$pred <- predict(svm_model, grid)

ggplot() +
  geom_point(data = data %>% filter(nc == 5), 
             aes(x = delta1, y = delta2, color = ifelse(pow1 > 0.95, 1, 0))) +
  geom_contour(data = df_f, aes(x = delta1, y = delta2, z = as.numeric(pred)), breaks = 0.95)  +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$"),
       fill = "nc") +
  theme_bw() +
  theme(axis.title = element_text(size = 13)) +
  xlim(0, 0.2) + ylim(0, 0.2)

#-------------- tile plot

data %>% filter(nc == 5) %>%
  ggplot() +
  geom_tile(aes(x = delta1, y = delta2, fill = pow2))  +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
     y = TeX("$\\delta_2 = \\mu_h - \\mu_c$"),
     fill = "nc") +
  theme_bw() +
  theme(axis.title = element_text(size = 13)) +
  xlim(0, 0.2) + ylim(0, 0.2) +
  scale_fill_gradient(low = "red", high = "green")

#-------------- contour fill plot

data %>% filter(nc == 5) %>%
  ggplot() +
  geom_contour_fill(aes(x = delta1, y = delta2, z = pow2))  +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$"),
       fill = "nc") +
  theme_bw() +
  theme(axis.title = element_text(size = 13)) +
  xlim(0, 0.2) + ylim(0, 0.2) +
  scale_fill_gradient(low = "red", high = "green")


#--------- 3d surface

library(plotly)
tex <- function(x) {
  structure(x, class = "LaTeX")
}
fig <- plot_ly(x = ~seq(0, 0.2, 0.02), y = ~seq(0, 0.2, 0.02), 
               z = ~matrix(data$pow1, nrow = 11), 
               color = ~nc, colors = viridis::viridis(length(unique(data$nc))),
               type = 'surface')

fig <- fig %>%
  layout(scene = list(
    xaxis = list(title = tex("\\delta_1 = \\mu_t - \\mu_c")),
    yaxis = list(title = tex("\\delta_2 = \\mu_h - \\mu_c")),
    zaxis = list(title = 'Power')
  ))
fig



#----------3d scatter

library(plotly)
fig <- plot_ly(data, x = ~delta1, y = ~delta2, z = ~pow1, color = ~nc, colors = viridis::viridis(length(unique(data$nc))),
               type = 'scatter3d', mode = 'markers')
fig <- fig %>%
  layout(scene = list(
    xaxis = list(title = TeX("$\\delta_1 = \\mu_t - \\mu_c$")),
    yaxis = list(title = TeX("$\\delta_2 = \\mu_h - \\mu_c$")),
    zaxis = list(title = 'Power')
  ))
fig







library(tidyverse)
library(plotly)
library(stringr)
library(latex2exp)
library(cowplot)

final_df_1 <- read_csv("results/map_results.csv") 
final_df_2 <- read_csv("results/normalized_results.csv")
final_df_3 <- read_csv("results/commensurate_results.csv")
final_df_4 <- read_csv("results/elastic_results.csv")

# ESS
ess_df <- final_df_1 %>% select("delta1","delta2","MAP" = "ess")
ess_df$Normalized_Power = final_df_2$ess
ess_df$Commensurate_Power = final_df_3$ess
ess_df$Elastic = final_df_4$ess
ess_df %>% pivot_longer(MAP:Elastic, names_to = "Method", values_to = "ESS") %>%
  ggplot(aes(x = delta2, y = ESS, color = Method)) +
  geom_line() + theme_bw() +
  labs(x = TeX("$\\delta_2 = \\mu_h - \\mu_c$")) + geom_vline(xintercept = 0) +
  scale_y_log10() +
  theme(axis.title = element_text(size = 13))

# Power
pow_df <- final_df_1 %>% select("delta1","delta2","MAP" = "pow")
pow_df$Normalized_Power = final_df_2$pow
pow_df$Commensurate_Power = final_df_3$pow
pow_df$Elastic = final_df_4$pow
pow_df %>% pivot_longer(MAP:Elastic, names_to = "Method", values_to = "Power") %>% 
  ggplot(aes(x = delta1, y = delta2, fill = Power)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "red", high = "green") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$")) +
  theme_bw() +
  facet_wrap(~Method) +
  theme(axis.title = element_text(size = 13))

