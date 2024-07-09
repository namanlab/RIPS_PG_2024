library(shiny)
library(shinydashboard)
library(latex2exp)
library(tidyverse)
library(e1071)
library(stringr)
library(plotly)
library(reshape2)

# Sample data to use in the visualizations
data_elastic_tau_1 <- read.csv("results/elastic_results_tau1.csv")
data_elastic_tau_2 <- read.csv("results/elastic_results_tau2.csv")
data_normalized_tau_1 <- read.csv("results/elastic_results_tau1.csv")
data_normalized_tau_2 <- read.csv("results/elastic_results_tau2.csv")




  

# 🥺
# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "BDB Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Elastic", tabName = "elastic", icon = icon("dashboard")),
      menuItem("Normalized Power", tabName = "normalized_power", icon = icon("dashboard")),
      selectInput("tau", "Select tau:", choices = unique(data_elastic_tau_1$tau)),
      selectInput("feature", "Select Feature:", choices = unique(data_elastic_tau_2$p)),
      selectInput("var", "Select Variable:", choices = c("Power", "PESS", "MSE"))
    )
  ),
  dashboardBody(
    withMathJax(),
    tabItems(
      tabItem(tabName = "elastic",
              h3("Elastic Method"),
              fluidRow(
                box(plotlyOutput("Plot1Elastic"), width = 6),
                box(plotlyOutput("Plot2Elastic"), width = 6)
              ),
              fluidRow(
                box(plotlyOutput("Plot3Elastic"), width = 12)
              )
      ),
      tabItem(tabName = "normalized_power",
              h3("Normalized Power Method"),
              fluidRow(
                box(plotlyOutput("Plot1NormalizedPower"), width = 6),
                box(plotlyOutput("Plot2NormalizedPower"), width = 6)
              ),
              fluidRow(
                box(plotlyOutput("Plot3NormalizedPower"), width = 12)
              )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
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
  
  
  
  output$Plot1Elastic <- renderPlotly({ render_p1(data_elastic_tau_1, input$var) })
  output$Plot2Elastic <- renderPlotly({ render_p2(data_elastic_tau_2, input$feature) })
  output$Plot3Elastic <- renderPlotly({ render_p3(data_elastic_tau_2, input$tau) })
  
  output$Plot1NormalizedPower <- renderPlotly({ render_p1(data_normalized_tau_1, input$var) })
  output$Plot2NormalizedPower <- renderPlotly({ render_p2(data_normalized_tau_2, input$feature) })
  output$Plot3NormalizedPower <- renderPlotly({ render_p3(data_normalized_tau_2, input$tau) })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
