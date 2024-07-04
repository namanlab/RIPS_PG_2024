library(shiny)
library(shinydashboard)
library(latex2exp)
library(tidyverse)
library(e1071)
library(stringr)
library(plotly)

# Sample data to use in the visualizations
data_elastic_power <- read.csv("results/elastic_power_results_nc.csv")
data_elastic <- read.csv("results/elastic_power_results_nc.csv")
data_commensurate <- read.csv("results/elastic_power_results_nc.csv")
data_normalized <- read.csv("results/elastic_power_results_nc.csv")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Respected Honourable Dr. Wang: Free Lunch Please 🥺"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Elastic Power", tabName = "elastic_power", icon = icon("dashboard")),
      menuItem("Elastic", tabName = "elastic", icon = icon("dashboard")),
      menuItem("Normalized Power", tabName = "normalized_power", icon = icon("dashboard")),
      menuItem("Commensurate Power", tabName = "commensurate_power", icon = icon("dashboard")),
      selectInput("nc", "Select nc:", choices = unique(data_elastic_power$nc)),
      sliderInput("delta1", "Select delta1:", min = 0, max = 0.2, value = 0.1, step = 0.02),
      sliderInput("delta2", "Select delta2:", min = 0, max = 0.2, value = 0.1, step = 0.02),
      selectInput("var", "Select Variable:", choices = c("Pow1", "Pow2", "PESS"))
    )
  ),
  dashboardBody(
    withMathJax(), # Ensure MathJax is included
    tabItems(
      tabItem(tabName = "elastic_power",
              h3("Elastic Power Method"),
              fluidRow(
                box(plotlyOutput("tilePlotElasticPower"), width = 6),
                box(plotlyOutput("lineChartElasticPower"), width = 6)
              ),
              fluidRow(
                box(plotlyOutput("scatter3dElasticPower"), width = 12)
              )
      ),
      tabItem(tabName = "elastic",
              h3("Elastic Method"),
              fluidRow(
                box(plotlyOutput("tilePlotElastic"), width = 6),
                box(plotlyOutput("lineChartElastic"), width = 6)
              ),
              fluidRow(
                box(plotlyOutput("scatter3dElastic"), width = 12)
              )
      ),
      tabItem(tabName = "normalized_power",
              h3("Normalized Power Method"),
              fluidRow(
                box(plotlyOutput("tilePlotNormalizedPower"), width = 6),
                box(plotlyOutput("lineChartNormalizedPower"), width = 6)
              ),
              fluidRow(
                box(plotlyOutput("scatter3dNormalizedPower"), width = 12)
              )
      ),
      tabItem(tabName = "commensurate_power",
              h3("Commensurate Power Method"),
              fluidRow(
                box(plotlyOutput("tilePlotCommensuratePower"), width = 6),
                box(plotlyOutput("lineChartCommensuratePower"), width = 6)
              ),
              fluidRow(
                box(plotlyOutput("scatter3dCommensuratePower"), width = 12)
              )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  render_tile_plot <- function(data, var, title) {
    df_cur <- data %>% filter(nc == input$nc)
    
    fill_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess")
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
  
  render_line_chart <- function(data, var) {
    df_cur <- data %>%
      filter(delta1 == input$delta1, delta2 == input$delta2)
    
    y_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess")
    y_label <- ifelse(var == "PESS", "ESS", "Power")
    
    p_cur <- ggplot(df_cur, aes(x = nc, y = !!sym(y_var))) +
      labs(x = "nc", y = y_label) +
      geom_line() + geom_point() +
      theme_minimal() +
      theme(axis.title = element_text(size = 13))
    
    ggplotly(p_cur)
  }
  
  render_scatter3d <- function(data, var) {
    fill_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess")
    df_cur <- data 
    fig <- plot_ly(df_cur, x = ~delta1, y = ~delta2, z = as.formula(paste("~", fill_var)), 
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
  
  output$tilePlotElasticPower <- renderPlotly({ render_tile_plot(data_elastic_power, input$var, "Elastic Power Method") })
  output$lineChartElasticPower <- renderPlotly({ render_line_chart(data_elastic_power, input$var) })
  output$scatter3dElasticPower <- renderPlotly({ render_scatter3d(data_elastic_power, input$var) })
  
  output$tilePlotElastic <- renderPlotly({ render_tile_plot(data_elastic, input$var, "Elastic Method") })
  output$lineChartElastic <- renderPlotly({ render_line_chart(data_elastic, input$var) })
  output$scatter3dElastic <- renderPlotly({ render_scatter3d(data_elastic, input$var) })
  
  output$tilePlotNormalizedPower <- renderPlotly({ render_tile_plot(data_normalized, input$var, "Normalized Power Method") })
  output$lineChartNormalizedPower <- renderPlotly({ render_line_chart(data_normalized, input$var) })
  output$scatter3dNormalizedPower <- renderPlotly({ render_scatter3d(data_normalized, input$var) })
  
  output$tilePlotCommensuratePower <- renderPlotly({ render_tile_plot(data_commensurate, input$var, "Commensurate Power Method") })
  output$lineChartCommensuratePower <- renderPlotly({ render_line_chart(data_commensurate, input$var) })
  output$scatter3dCommensuratePower <- renderPlotly({ render_scatter3d(data_commensurate, input$var) })
}

# Run the application 
shinyApp(ui = ui, server = server)
