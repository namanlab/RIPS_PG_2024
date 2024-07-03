library(tidyverse)
library(plotly)
library(stringr)
library(latex2exp)
library(cowplot)

final_df_1 <- read_csv("results/map_results.csv") 
final_df_2 <- read_csv("results/normalized_results.csv")
final_df_3 <- read_csv("results/commensurate_results.csv")
final_df_4 <- read_csv("results/elastic_results.csv")
final_df_5 <- read_csv("results/elastic_power_results.csv")

# ESS
ess_df <- final_df_1 %>% select("delta1","delta2","MAP" = "ess")
ess_df$Normalized_Power = final_df_2$ess
ess_df$Commensurate_Power = final_df_3$ess
ess_df$Elastic = final_df_4$ess
ess_df$Elastic_Power = final_df_5$ess
ess_df %>% pivot_longer(MAP:Elastic_Power, names_to = "Method", values_to = "ESS") %>%
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
pow_df$Elastic_Power = final_df_5$pow
pow_df %>% pivot_longer(MAP:Elastic_Power, names_to = "Method", values_to = "Power") %>% 
  ggplot(aes(x = delta1, y = delta2, fill = Power)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "red", high = "green") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$")) +
  theme_bw() +
  facet_wrap(~Method) +
  theme(axis.title = element_text(size = 13))

