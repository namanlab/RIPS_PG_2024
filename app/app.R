library(shiny)
library(shinydashboard)
library(latex2exp)
library(tidyverse)
library(e1071)
library(stringr)
library(plotly)
library(reshape2)

# Sample data to use in the visualizations
data_elastic_tau_1 <- read.csv("results_fc/elastic_results_tau1.csv")
data_elastic_tau_2 <- read.csv("results_fc/elastic_results_tau2.csv")
data_normalized_tau_1 <- read.csv("results_fc/elastic_results_tau1.csv")
data_normalized_tau_2 <- read.csv("results_fc/elastic_results_tau2.csv")

data_elastic_power <- read.csv("results_oc/elastic_power_results_nc.csv")
data_elastic <- read.csv("results_oc/elastic_results_nc.csv")
data_commensurate <- read.csv("results_oc/commensurate_results_nc.csv")
data_normalized <- read.csv("results_oc/normalized_results_nc.csv")
data_rmap <- read.csv("results_oc/rMAP_results_nc.csv")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "BDB Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Face Cream", tabName = "face_cream", icon = icon("smile")),
      menuItem("Oral Care", tabName = "oral_care", icon = icon("tooth")),
      selectInput("dataset", "Select Dataset:", choices = c("Face Cream", "Oral Care"))
    ),
    uiOutput("dynamic_sidebar")
  ),
  dashboardBody(
    withMathJax(),
    uiOutput("dynamic_body")
  )
)

# Define server logic
server <- function(input, output, session) {
  
  
  
  output$dynamic_sidebar <- renderUI({
    if (input$dataset == "Face Cream") {
      sidebarMenu(
        selectInput("tau", "Select tau:", choices = unique(data_elastic_tau_1$tau)),
        selectInput("feature", "Select Feature:", choices = unique(data_elastic_tau_2$p)),
        selectInput("var", "Select Variable:", choices = c("Power", "PESS", "MSE"))
      )
    } else {
      sidebarMenu(
        selectInput("nc", "Select nc:", choices = unique(data_elastic_power$nc)),
        sliderInput("delta1", "Select delta1:", min = 0, max = 0.2, value = 0.1, step = 0.02),
        sliderInput("delta2", "Select delta2:", min = 0, max = 0.2, value = 0.1, step = 0.02),
        selectInput("var", "Select Variable:", choices = c("Pow1", "Pow2", "PESS"))
      )
    }
  })
  
  output$dynamic_body <- renderUI({
    if (input$dataset == "Face Cream") {
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
      
    } else {
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
        ),
        tabItem(tabName = "robust_map",
                h3("Robust MAP Method"),
                fluidRow(
                  box(plotlyOutput("tilePlotRMAP"), width = 6),
                  box(plotlyOutput("lineChartRMAP"), width = 6)
                ),
                fluidRow(
                  box(plotlyOutput("scatter3dRMAP"), width = 12)
                )
        )
      )
    }
  })
  
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
    y_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess")
    if (var == "PESS"){y_label <- "ESS"}
    else {y_label <- "Power"}
    p_cur <- ggplot(data %>% filter(delta1 == input$delta1, delta2 == input$delta2), aes(x = nc, y = !!sym(y_var))) +
      labs(x = "nc", y = y_label) +
      geom_line() + geom_point() +
      theme_minimal() +
      theme(axis.title = element_text(size = 13))
    
    ggplotly(p_cur)
  }
  
  render_scatter3d <- function(data, var) {
    fill_var <- switch(var, "Pow1" = "pow1", "Pow2" = "pow2", "PESS" = "ess")
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
  
  output$tilePlotRMAP <- renderPlotly({ render_tile_plot(data_rmap, input$var, "Robust MAP Method") })
  output$lineChartRMAP <- renderPlotly({ render_line_chart(data_rmap, input$var) })
  output$scatter3dRMAP <- renderPlotly({ render_scatter3d(data_rmap, input$var) })
}

# Run the application 
shinyApp(ui = ui, server = server)
