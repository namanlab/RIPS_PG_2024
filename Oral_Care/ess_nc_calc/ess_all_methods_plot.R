library(latex2exp)
# Load datasets for different methods
commensurate <- read.csv("results/commensurate_results_nc.csv")
elastic <- read.csv("results/elastic_results_nc.csv")
map <- read.csv("results/rMAP_results_nc.csv")
norm <- read.csv("results/normalized_results_nc.csv")
elastic_pow <- read.csv("results/elastic_power_results_nc.csv")

clab <- rep("Commensurate", 726)
elab <- rep("Elastic", 726)
mlab <- rep("Robust MAP", 726)
nlab <- rep("Normalised", 726)
eplab <- rep("Elastic Power", 726)

commensurate$method <- clab
elastic$method <- elab
map$method <- mlab
norm$method <- nlab
elastic_pow$method <- eplab

df <- rbind(commensurate, elastic, map, norm, elastic_pow)
df2 <- df %>%
  filter(nc == 30, delta1 == 0)

# library(ggplot2)
# ggplot(df, aes(x = delta2, y = ess, color = method)) +
#   geom_line() +
#   labs(
#     x = TeX("$\\delta_2 = \\mu_h - \\mu_c$"),
#     y = "ESS"
#   )


library(plotly)
library(dplyr)
library(htmlwidgets)

plot <- plot_ly(data = df2, 
                x = ~delta2, 
                y = ~ess, 
                color = ~method, 
                type = 'scatter', 
                mode = 'lines') %>%
  layout(
    xaxis = list(title = TeX("$\\delta_2 = \\theta_h - \\theta_c$")),
    yaxis = list(title = "ESS")
  )

plot <- config(plot, mathjax = 'cdn')

saveWidget(plot, "plot.html")
browseURL("plot.html")


# Delta2 = history - control
# Delta1 = treatment - control















