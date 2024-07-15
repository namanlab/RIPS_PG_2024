setwd("~/Desktop/RIPS_PG_2024")
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::activate()


install.packages(c("shiny", "shinydashboard", "latex2exp", "tidyverse", "e1071", 
                   "stringr", "plotly", "reshape2", "bslib", "readxl", 
                   "MCMCpack", "LaplacesDemon", "invgamma", "HDInterval", "mvtnorm"))
install.packages("viridis")
install.packages("R.utils")

renv::snapshot()
