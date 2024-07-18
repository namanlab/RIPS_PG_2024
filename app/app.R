source("helper_oc.R")
source("helper_fc.R")
source("helper.R")


library(shiny)
library(shinydashboard)
library(latex2exp)
library(tidyverse)
library(e1071)
library(stringr)
library(plotly)
library(reshape2)
library(bslib)
library(readxl)
library(MCMCpack)
library(LaplacesDemon)
library(invgamma)
library(tidyverse)
library(HDInterval)
library(mvtnorm)
library(viridis)
library(R.utils)

# Sample data to use in the visualizations
data_elastic_tau_1 <- read.csv("results_fc/elastic_results_tau1_updated.csv") 
data_elastic_tau_2 <- read.csv("results_fc/elastic_results_tau2_updated.csv")
data_elastic_fc <- read.csv("results_fc/elastic_results_nc_fc.csv")

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
      selectInput("dataset", "Select Dataset:", choices = c("Oral Care", "Face Cream"))
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
  
  observeEvent(input$submit_OC, {
    req(input$datafile)
    
    data <- read_data_OC(input$datafile$datapath)
    yc <- data$yc
    yt <- data$yt
    yh <- data$yh
    
    # Compute results for each method
    bayesian_probs <- compute_bayesian_probs_OC(yc, yh, yt, input$R_samples)
    print(bayesian_probs)
    
    # Frequentist method
    frequentist_p_value <- compute_oc_frequentist(yc, yt)
    
    output$results_elastic_power_OC <- renderUI({
      generate_result_ui("Elastic Power Method", bayesian_probs$elastic_power, frequentist_p_value)
    })
    output$results_elastic_OC <- renderUI({
      generate_result_ui("Elastic Method", bayesian_probs$elastic, frequentist_p_value)
    })
    output$results_normalized_OC <- renderUI({
      generate_result_ui("Normalized Power Method", bayesian_probs$normalized, frequentist_p_value)
    })
    output$results_commensurate_OC <- renderUI({
      generate_result_ui("Commensurate Power Method", bayesian_probs$commensurate, frequentist_p_value)
    })
    output$results_robust_map_OC <- renderUI({
      generate_result_ui("Robust MAP Method", bayesian_probs$robust_map, frequentist_p_value)
    })
  })
  
  
  observeEvent(input$submit_sim_OC, {
    req(input$datafile)
    
    data <- read_excel(input$datafile$datapath)
    yc <- data[[1]]
    yt <- data[[2]]
    
    H <- input$H
    yh_list <- lapply(1:H, function(i) {
      muh <- input[[paste0("muh_", i)]]
      sigma2_h <- input[[paste0("sigma2_h_", i)]]
      nh <- input[[paste0("nh_", i)]]
      rnorm(nh, muh, sqrt(sigma2_h))
    })
    yh <- unlist(yh_list)
    
    # Compute results for each method
    bayesian_probs <- compute_bayesian_probs_OC(yc, yh, yt, input$R_samples)
    print(bayesian_probs)
    
    # Frequentist method
    frequentist_p_value <- compute_oc_frequentist(yc, yt)
    
    output$results_elastic_power_OC <- renderUI({
      generate_result_ui("Elastic Power Method", bayesian_probs$elastic_power, frequentist_p_value)
    })
    output$results_elastic_OC <- renderUI({
      generate_result_ui("Elastic Method", bayesian_probs$elastic, frequentist_p_value)
    })
    output$results_normalized_OC <- renderUI({
      generate_result_ui("Normalized Power Method", bayesian_probs$normalized, frequentist_p_value)
    })
    output$results_commensurate_OC <- renderUI({
      generate_result_ui("Commensurate Power Method", bayesian_probs$commensurate, frequentist_p_value)
    })
    output$results_robust_map_OC <- renderUI({
      generate_result_ui("Robust MAP Method", bayesian_probs$robust_map, frequentist_p_value)
    })
  })
  
  observeEvent(input$simulate_OC, {
    muc <- input$muc
    mut <- input$mut
    sigma2_c <- input$sigma2_c
    sigma2_t <- input$sigma2_t
    nc <- input$nc
    nt <- input$nt
    H <- input$H
    
    yh_list <- lapply(1:H, function(i) {
      muh <- input[[paste0("muh_", i)]]
      sigma2_h <- input[[paste0("sigma2_h_", i)]]
      nh <- input[[paste0("nh_", i)]]
      rnorm(nh, muh, sqrt(sigma2_h))
    })
    yc <- rnorm(nc, muc, sqrt(sigma2_c))
    yt <- rnorm(nt, mut, sqrt(sigma2_t))
    yh <- unlist(yh_list)
    
    # Compute results for each method
    bayesian_probs <- compute_bayesian_probs_OC(yc, yh, yt, input$R_samples)
    print(bayesian_probs)
    
    # Frequentist method
    frequentist_p_value <- compute_oc_frequentist(yc, yt)
    
    output$results_elastic_power_OC <- renderUI({
      generate_result_ui("Elastic Power Method", bayesian_probs$elastic_power, frequentist_p_value)
    })
    output$results_elastic_OC <- renderUI({
      generate_result_ui("Elastic Method", bayesian_probs$elastic, frequentist_p_value)
    })
    output$results_normalized_OC <- renderUI({
      generate_result_ui("Normalized Power Method", bayesian_probs$normalized, frequentist_p_value)
    })
    output$results_commensurate_OC <- renderUI({
      generate_result_ui("Commensurate Power Method", bayesian_probs$commensurate, frequentist_p_value)
    })
    output$results_robust_map_OC <- renderUI({
      generate_result_ui("Robust MAP Method", bayesian_probs$robust_map, frequentist_p_value)
    })
  })
    
  
  observeEvent(input$submit_FC, {
    req(input$datafile)
    req(input$datafile_hist)
    data <- read_data_FC(input$datafile$datapath, input$datafile_hist$datapath)
    Zc <- data$Zc
    Xc <- data$Xc
    xh <- data$xh
    xc <- data$xc
    
    # Compute results for each method
    bayesian_probs <- compute_bayesian_probs_FC(Xc, Zc, xc, xh, input$R_samples)
    bayesian_probs2 <- compute_bayesian_probs_FC_nohist(Xc, Zc, xc, xh, input$R_samples)
    
    output$results_elastic_FC <- renderUI({
      generate_result_ui_FC("Elastic Method", bayesian_probs$elastic, bayesian_probs2$elastic)
    })
    
  })
  
  observeEvent(input$simulate_FC, {
    
    nc <- input$n
    pc <- input$p_feats
    
    uc_list <- lapply(1:pc, function(i) {
      input[[paste0("mut_", i)]]
    })
    uc <- unlist(uc_list)
    
    Zc <- gen_Z(nc)
    Xc <- gen_X(nc, pc, input$nc)
    sig = sqrt(input$sigma2)
    tau = sqrt(input$tau2)
    
    H <- input$H
    yh_list <- lapply(1:H, function(i) {
      muh <- input[[paste0("muh_", i)]]
      sigma2_h <- input[[paste0("sigma2_h_", i)]]
      nh <- input[[paste0("nh_", i)]]
      rnorm(nh, muh, sqrt(sigma2_h))
    })
    xh <- unlist(yh_list)
    
    overall_noise <- rnorm(2*nc, 0, sig^2)
    indiv_noise <- rnorm(nc, 0, tau^2)
    xc <- Xc %*% uc + Zc %*% indiv_noise + overall_noise
    
    # Compute results for each method
    bayesian_probs <- compute_bayesian_probs_FC(Xc, Zc, xc, xh, input$R_samples)
    bayesian_probs2 <- compute_bayesian_probs_FC_nohist(Xc, Zc, xc, xh, input$R_samples)
    
    output$results_elastic_FC <- renderUI({
      generate_result_ui_FC("Elastic Method", bayesian_probs$elastic, bayesian_probs2$elastic)
    })
    
  })
  
  
  # This is to get the desired menuItem selected initially. 
  # selected=T seems not to work with a dynamic sidebarMenu.
  observeEvent(session, {
    updateTabItems(session, "tabs", selected = "initial")
  })
  
  output$dynamic_sidebar <- renderUI({
    if (input$dataset == "Face Cream") {
      sidebarMenu(
        selectInput("anal_type", "Use:", choices = c("Power Analysis", "Tau Analysis", "Import Data", "Run Simulation")),
        uiOutput("dynamic_sidebar_sub_fc")
      )
    } else {
      sidebarMenu(
        selectInput("anal_type", "Use:", choices = c("Power Analysis", "Import All Data", "Run Simulation", "Import Current Data")),
        uiOutput("dynamic_sidebar_sub_oc")
      )
    }
  })
  
  
  
  output$dynamic_sidebar_sub_oc <- renderUI({
    if (input$anal_type == "Power Analysis") {
      sidebarMenu(
        selectInput("nc", "Select nc:", choices = unique(data_elastic_power$nc)),
        sliderInput("delta1", "Select delta1:", min = 0, max = 0.2, value = 0.1, step = 0.02),
        sliderInput("delta2", "Select delta2:", min = 0, max = 0.2, value = 0.1, step = 0.02),
        selectInput("var", "Select Variable:", choices = c("Pow1", "Pow2", "PESS"))
      )
    } else if (input$anal_type == "Import All Data") {
      sidebarMenu(
        numericInput("R_samples", "No. of Samples from Prosterior Distr.", value = 10000),
        tags$div(title=paste0("Data Import Format\n\n",
                              "When uploading an Excel file, please ensure the data is organized in the following format:\n\n",
                              "1. First Column (yc): This column should contain the data for the control group (`yc`). Each value represents a sample from the control group.\n\n",
                              "2. Second Column (yt): This column should contain the data for the treatment group (`yt`). Each value represents a sample from the treatment group.\n\n",
                              "3. Subsequent Columns (yh1, yh2, ... yhn): These columns should contain data for the historical datasets (`yh`). Each column represents a separate historical dataset, with values representing samples from each respective dataset.\n"),
          fileInput("datafile", "Upload Excel File", accept = c(".xlsx"))
        ),
        actionButton("submit_OC", "Submit")
      )
    } else if (input$anal_type == "Import Current Data") {
      sidebarMenu(
        numericInput("R_samples", "No. of Samples from Prosterior Distr.", value = 10000),
        tags$div(title=paste0("Data Import Format\n\n",
                              "When uploading an Excel file, please ensure the data is organized in the following format:\n\n",
                              "1. First Column (yc): This column should contain the data for the control group (`yc`). Each value represents a sample from the control group.\n\n",
                              "2. Second Column (yt): This column should contain the data for the treatment group (`yt`). Each value represents a sample from the treatment group.\n\n"),
                 fileInput("datafile", "Upload Excel File", accept = c(".xlsx"))
        ),
        numericInput("H", "Number of historical datasets", value = 1, min = 1),
        uiOutput("historical_data_inputs"),
        actionButton("submit_sim_OC", "Submit & Simulate Historical")
      )
    } else {
      sidebarMenu(
        numericInput("R_samples", "No. of Samples from Prosterior Distr.", value = 10000),
        numericInput("muc", "Mean of yc", value = 2.5),
        numericInput("mut", "Mean of yt", value = 2.4),
        numericInput("sigma2_c", "Variance of yc", value = 0.22),
        numericInput("sigma2_t", "Variance of yt", value = 0.23),
        numericInput("nc", "Number of samples for yc", value = 30),
        numericInput("nt", "Number of samples for yt", value = 29),
        numericInput("H", "Number of historical datasets", value = 1, min = 1),
        uiOutput("historical_data_inputs"),
        actionButton("simulate_OC", "Simulate Data")
      )
    }
  })
  
  observe({
    H <- input$H
    set.seed(42)
    output$historical_data_inputs <- renderUI({
      lapply(1:H, function(i) {
        tagList(
          numericInput(paste0("muh_", i), paste0("Mean of yh_", i), value = 2.45 + rnorm(1, 0, 0.05)),
          numericInput(paste0("sigma2_h_", i), paste0("Variance of yh_", i), value = 0.3),
          numericInput(paste0("nh_", i), paste0("Number of samples for yh_", i), value = 31),
          hr() # Optional: Adds a horizontal line to separate each set of inputs
        )
      })
    })
  })
  
  output$dynamic_sidebar_sub_fc <- renderUI({
    if (input$anal_type == "Power Analysis") {
      sidebarMenu(
        selectInput("nc", "Select nc:", choices = unique(data_elastic_fc$nc)),
        sliderInput("delta1", "Select delta1:", min = 0, max = 0.2, value = 0.1, step = 0.02),
        sliderInput("delta2", "Select delta2:", min = 0, max = 0.2, value = 0.1, step = 0.02),
        selectInput("var", "Select Variable:", choices = c("Power", "PESS"))
      )
    } else if (input$anal_type == "Tau Analysis") {
      sidebarMenu(
        selectInput("tau", "Select tau:", choices = unique(data_elastic_tau_1$tau)),
        selectInput("feature", "Select Feature:", choices = unique(data_elastic_tau_2$p)),
        selectInput("var", "Select Variable:", choices = c("Power", "PESS", "MSE"))
      )
    } else if (input$anal_type == "Import Data") {
      sidebarMenu(
        numericInput("R_samples", "No. of Samples from Prosterior Distr.", value = 10000),
        tags$div(title=paste0("Data Import Format\n\n",
                              "When uploading an Excel file, please ensure the data is organized in the following format:\n\n",
                              "1. First Column (id): This column should contain the face id (unique identifier for the face/indivdual the sample comes from)\n\n",
                              "2. Second Column (side): This column should contain either L (left) or R (right) denoting the side of the face from which the sample came.\n\n",
                              "3. Third Column (treatment_grp): This column should contain the treatment id/name (unique identifier for the treatment/control applied on the sample)\n\n",
                              "4. Fourth Column (y): This column should contain the measured value/observation (after adjusting for the baseline).\n\n"),
                 fileInput("datafile", "Upload Excel File (Current Data)", accept = c(".xlsx", ".csv"))
        ), 
        tags$div(title=paste0("Data Import Format\n\n",
                              "When uploading an Excel file, please ensure the data is organized in the following format:\n\n",
                              "1. Historical Data (yh1, yh2, ... yhn): These columns should contain data for the historical datasets (`yh`). Each column represents a separate historical dataset, with values representing samples from each respective dataset.\n"),
                 fileInput("datafile_hist", "Upload Excel File (Historical Data)", accept = c(".xlsx", ".csv"))
        ),
        actionButton("submit_FC", "Submit")
      )
    } else {
      sidebarMenu(
        numericInput("R_samples", "No. of Samples from Prosterior Distr.", value = 10000),
        numericInput("p_feats", "Number of Treatments (Including Control)", value = 3, min = 3),
        uiOutput("feature_data_inputs"),
        numericInput("sigma2", "Between Sample variance", value = 0.2^2),
        numericInput("tau2", "Within Sample Variance", value = 0.2^2),
        numericInput("n", "Number of Samples", value = 99),
        numericInput("nc", "Number of Control Faces", value = 29),
        numericInput("H", "Number of historical datasets", value = 1, min = 1),
        uiOutput("historical_data_inputs"),
        actionButton("simulate_FC", "Simulate Data")
      )
    }
  })
  
  observe({
    p <- input$p_feats
    output$feature_data_inputs <- renderUI({
      set.seed(42)
      lapply(1:p, function(i) {
        if (i == p){tagList(
          numericInput(paste0("mut_", i), paste0("Mean of control"), value = 2.45 - rnorm(1, 0.1, 0.05)),
          hr() # Optional: Adds a horizontal line to separate each set of inputs
        )} else {tagList(
          numericInput(paste0("mut_", i), paste0("Mean of treatment ", i), value = 2.4 + rnorm(1, 0, 0.05)),
          hr() # Optional: Adds a horizontal line to separate each set of inputs
        )}
      })
    })
  })
  
  
  # Dynamic body with navset_card_underline
  output$dynamic_body <- renderUI({
    if (input$dataset == "Face Cream") {
      
      if (input$anal_type == "Tau Analysis") {
        navset_tab(
          nav_panel(
            "Elastic Method",
            fluidRow(
              box(plotlyOutput("Plot1Elastic"), width = 6),
              box(plotlyOutput("Plot2Elastic"), width = 6)
            ),
            fluidRow(
              box(plotlyOutput("Plot3Elastic"), width = 12)
            )
          )
        )
      } else if (input$anal_type == "Power Analysis") {
        navset_tab(
          nav_panel(
            "Elastic Method",
            fluidRow(
              box(plotlyOutput("tilePlotElasticFC"), width = 6),
              box(plotlyOutput("lineChartElasticFC"), width = 6)
            ),
            fluidRow(
              box(plotlyOutput("scatter3dElasticFC"), width = 12)
            )
          )
        )
      } else {
        navset_tab(
          nav_panel("Elastic Method", br(), uiOutput("results_elastic_FC"))
        )
      }
      
    } else {
      if (input$anal_type == "Power Analysis") {
        navset_tab(
          nav_panel(
            "Elastic Power Method",
            fluidRow(
              box(plotlyOutput("tilePlotElasticPower"), width = 6),
              box(plotlyOutput("lineChartElasticPower"), width = 6)
            ),
            fluidRow(
              box(plotlyOutput("scatter3dElasticPower"), width = 12)
            )
          ),
          nav_panel(
            "Elastic Method",
            fluidRow(
              box(plotlyOutput("tilePlotElastic"), width = 6),
              box(plotlyOutput("lineChartElastic"), width = 6)
            ),
            fluidRow(
              box(plotlyOutput("scatter3dElastic"), width = 12)
            )
          ),
          nav_panel(
            "Normalized Power Method",
            fluidRow(
              box(plotlyOutput("tilePlotNormalizedPower"), width = 6),
              box(plotlyOutput("lineChartNormalizedPower"), width = 6)
            ),
            fluidRow(
              box(plotlyOutput("scatter3dNormalizedPower"), width = 12)
            )
          ),
          nav_panel(
            "Commensurate Power Method",
            fluidRow(
              box(plotlyOutput("tilePlotCommensuratePower"), width = 6),
              box(plotlyOutput("lineChartCommensuratePower"), width = 6)
            ),
            fluidRow(
              box(plotlyOutput("scatter3dCommensuratePower"), width = 12)
            )
          ),
          nav_panel(
            "Robust MAP Method",
            fluidRow(
              box(plotlyOutput("tilePlotRMAP"), width = 6),
              box(plotlyOutput("lineChartRMAP"), width = 6)
            ),
            fluidRow(
              box(plotlyOutput("scatter3dRMAP"), width = 12)
            )
          )
        )
      } else {
        navset_tab(
          nav_panel("Elastic Power Method", br(), uiOutput("results_elastic_power_OC")),
          nav_panel("Elastic Method", br(), uiOutput("results_elastic_OC")),
          nav_panel( "Normalized Power Method", br(), uiOutput("results_normalized_OC")),
          nav_panel("Commensurate Power Method", br(), uiOutput("results_commensurate_OC")),
          nav_panel( "Robust MAP Method", br(), uiOutput("results_robust_map_OC"))
        )
      }
    }
  })
  
  
  
  
  output$Plot1Elastic <- renderPlotly({ render_p1(data_elastic_tau_1, input$var) })
  output$Plot2Elastic <- renderPlotly({ render_p2(data_elastic_tau_2, input$feature) })
  output$Plot3Elastic <- renderPlotly({ render_p3(data_elastic_tau_2, input$tau) })
  
  output$tilePlotElasticPower <- renderPlotly({ render_tile_plot(data_elastic_power, input$var, "Elastic Power Method", input$nc) })
  output$lineChartElasticPower <- renderPlotly({ render_line_chart(data_elastic_power, input$var, input$delta1, input$delta2) })
  output$scatter3dElasticPower <- renderPlotly({ render_scatter3d(data_elastic_power, input$var) })
  
  output$tilePlotElastic <- renderPlotly({ render_tile_plot(data_elastic, input$var, "Elastic Method", input$nc) })
  output$lineChartElastic <- renderPlotly({ render_line_chart(data_elastic, input$var, input$delta1, input$delta2) })
  output$scatter3dElastic <- renderPlotly({ render_scatter3d(data_elastic, input$var) })
  
  output$tilePlotNormalizedPower <- renderPlotly({ render_tile_plot(data_normalized, input$var, "Normalized Power Method", input$nc) })
  output$lineChartNormalizedPower <- renderPlotly({ render_line_chart(data_normalized, input$var, input$delta1, input$delta2) })
  output$scatter3dNormalizedPower <- renderPlotly({ render_scatter3d(data_normalized, input$var) })
  
  output$tilePlotCommensuratePower <- renderPlotly({ render_tile_plot(data_commensurate, input$var, "Commensurate Power Method", input$nc) })
  output$lineChartCommensuratePower <- renderPlotly({ render_line_chart(data_commensurate, input$var, input$delta1, input$delta2) })
  output$scatter3dCommensuratePower <- renderPlotly({ render_scatter3d(data_commensurate, input$var) })
  
  output$tilePlotRMAP <- renderPlotly({ render_tile_plot(data_rmap, input$var, "Robust MAP Method", input$nc) })
  output$lineChartRMAP <- renderPlotly({ render_line_chart(data_rmap, input$var, input$delta1, input$delta2) })
  output$scatter3dRMAP <- renderPlotly({ render_scatter3d(data_rmap, input$var) })
  
  output$tilePlotElasticFC <- renderPlotly({ render_tile_plot(data_elastic_fc, input$var, "Elastic Method", input$nc) })
  output$lineChartElasticFC <- renderPlotly({ render_line_chart(data_elastic_fc, input$var, input$delta1, input$delta2) })
  output$scatter3dElasticFC <- renderPlotly({ render_scatter3d(data_elastic_fc, input$var) })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
