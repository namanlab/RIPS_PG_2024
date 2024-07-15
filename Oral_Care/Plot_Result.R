library(tidyverse)  # Load the tidyverse package for data manipulation and visualization
library(plotly)  # Load the plotly package for interactive plots
library(stringr)  # Load the stringr package for string manipulation
library(latex2exp)  # Load the latex2exp package for LaTeX expressions
library(cowplot)  # Load the cowplot package for arranging plots

# Read the CSV files into data frames
final_df_1 <- read_csv("results/map_results.csv") 
final_df_2 <- read_csv("results/normalized_results.csv")
final_df_3 <- read_csv("results/commensurate_results.csv")
final_df_4 <- read_csv("results/elastic_results.csv")

# ESS plot
ess_df <- final_df_1 %>% select("delta1","delta2","MAP" = "ess")  # Select columns from final_df_1
ess_df$Normalized_Power = final_df_2$ess  # Add a column from final_df_2
ess_df$Commensurate_Power = final_df_3$ess  # Add a column from final_df_3
ess_df$Elastic = final_df_4$ess  # Add a column from final_df_4
ess_df %>% 
  pivot_longer(MAP:Elastic, names_to = "Method", values_to = "ESS") %>%  # Reshape the data frame
  ggplot(aes(x = delta2, y = ESS, color = Method)) +  # Create a ggplot object
  geom_line() + theme_bw() +  # Add a line plot with a black and white theme
  labs(x = TeX("$\\delta_2 = \\mu_h - \\mu_c$")) + geom_vline(xintercept = 0) +  # Add x-axis label and vertical line at x = 0
  scale_y_log10() +  # Use a logarithmic scale for the y-axis
  theme(axis.title = element_text(size = 13))  # Set the size of axis titles

# Power plot
pow_df <- final_df_1 %>% select("delta1","delta2","MAP" = "pow")  # Select columns from final_df_1
pow_df$Normalized_Power = final_df_2$pow  # Add a column from final_df_2
pow_df$Commensurate_Power = final_df_3$pow  # Add a column from final_df_3
pow_df$Elastic = final_df_4$pow  # Add a column from final_df_4
pow_df %>% 
  pivot_longer(MAP:Elastic, names_to = "Method", values_to = "Power") %>%  # Reshape the data frame
  ggplot(aes(x = delta1, y = delta2, fill = Power)) +  # Create a ggplot object
  geom_tile(color = "black") +  # Add a tiled plot with black borders
  scale_fill_gradient(low = "red", high = "green") +  # Use a color gradient for fill
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +  # Add horizontal and vertical lines at y = 0 and x = 0
  labs(x = TeX("$\\delta_1 = \\mu_t - \\mu_c$"), 
       y = TeX("$\\delta_2 = \\mu_h - \\mu_c$")) +  # Add x-axis and y-axis labels
  theme_bw() +  # Use a black and white theme
  facet_wrap(~Method) +  # Create multiple plots based on the Method variable
  theme(axis.title = element_text(size = 13))  # Set the size of axis titles
