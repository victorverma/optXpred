library(shiny)
library(DT)
library(plotly)

source("/home/sstoev/Dropbox/doc/GitHub/optXpred/opt_ext_pred.R")
source("/home/sstoev/Dropbox/doc/GitHub/optXpred/tail_dependence.R")

to_pareto_scale <- function(y){
  n =length(y);
  return(1/(1-rank(y)/(n+1)))
}

to_uniform_scale <- function(y){
  n =length(y);
  return(rank(y)/(n+1))
}



compute_metrics <- function(obs, pred) {
  if (length(obs) != length(pred)) stop("Observed and predicted vectors must be the same length.")
  
  obs <- as.integer(obs)
  pred <- as.integer(pred)
  
  TP <- sum(obs == 1 & pred == 1)
  TN <- sum(obs == 0 & pred == 0)
  FP <- sum(obs == 0 & pred == 1)
  FN <- sum(obs == 1 & pred == 0)
  
  sensitivity <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  specificity <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  TSS <- sensitivity + specificity - 1
  
  precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  recall <- sensitivity
  F1 <- if (is.na(precision) || is.na(recall) || (precision + recall == 0)) NA else {
    2 * precision * recall / (precision + recall)
  }
  
  list(
    TSS = TSS,
    prec = precision,
    recall = recall,
    missed = if ((TP + FN) == 0) NA else FN / (TP + FN),
    alarm = mean(pred),
    F1 = F1,
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN
  )
}

ui <- fluidPage(
  titlePanel("optXpred: Optimal Extremal Homogeneous Prediction"),
  
  p("This Shiny app illustrates the optimal homogeneous prediction methodology developed by Benjamin Bobbia and Stilian Stoev. It 
   allows the user to upload training and test datasets, select a response variable and predictor variables, 
    and compute an estimate of the optimal homogeneous predictor via quantile random forests using the R-package 'ranger'.  
    The predictions are graphed over the test data set and various metrics are displayed in a table such as TSS (True Skill Statistic), 
    Precision, F1 Score, and Miss Rate. The empirical tail-dependence of the response and the predictor are illustrated in another tab.
    NOTE: By default the App loads a training and testing data frames from daily maxima of 
    X-ray flux, and 23 other predictors of the so-called SHARP time series derived from active regions 
    of the Sun compiled and curated by Victor Verma. The prediction of Solar flares is challenging, so the please do not be surprised with
    the low scores.  Try it with your own data and let us know how it goes! (c) {sstoev,vkverma}@umich.edu and Benjamin.BOBBIA@isae-supaero.fr."),
  sidebarLayout(
    sidebarPanel(
      fileInput("train_file", "Upload Training Data (CSV)", accept = ".csv"),
      fileInput("test_file", "Upload Test Data (CSV)", accept = ".csv"),
      #sliderInput("p_range", "Precision plot range (p):", min = 0, max = 1, step = 0.005, value = c(0.6, 0.99)),
      radioButtons("extrapolate_tail", "GPD tail extrapolation", choices = c(TRUE, FALSE), selected = FALSE),
      radioButtons("verbose", "Verbose output", choices = c(TRUE, FALSE), selected = TRUE),
      radioButtons("log_scale", "Y-axis scale for predictions plot:", choices = c("Linear", "Log"), selected = "Linear"),
      numericInput("prop", "Proportion of extremes to consider", value = 0.9, step = 0.01, min = 0, max = 1),
      numericInput("alarm_rate", "Alarm rate (probability):", min = 0, max = 1, value = 0.2, step = 0.01),
      numericInput("event_rate", "Event rate (probability):", min = 0, max = 1, value = 0.2, step = 0.01),
      hr(),
      uiOutput("response_selector"),
      uiOutput("predictor_selector"),
      hr(),
      p("Select predictor variables manually or select a response and click the button below for automatic variable selection."),
      numericInput("keep_threshold", "Prop threshold for automatic variable selection",
                   value = 0.95, step = 0.01, min = 0, max = 1),
      numericInput("n_trees","Number of trees in the forest", value = 2000,min=100,max=20000),
      actionButton("suggest_vars", "Suggest Predictors (via Partial Tail Dependence)"),
      actionButton("run_model", "Run Model")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Predictions Plot",
                 p("Color coding:"),
                 tags$ul(
                   tags$li(strong("Blue:"), " Correct predictions (True Positives)"),
                   tags$li(strong("Red:"), " Missed events (False Negatives)"),
                   tags$li(strong("Black:"), " Other cases (True Negatives or False Positives)")
                 ),
                 plotlyOutput("predictions_plot")),
        tabPanel("True vs Predicted Scatter",
                 p("True Y vs Predicted Y scatter plot. Alarm and event thresholds shown."),
                 plotlyOutput("scatter_plot")),
        tabPanel("Metrics Table", DTOutput("metrics_table")),
        tabPanel("Tail Dependence", plotlyOutput("emp_precision")),
        tabPanel("Download Predictions", 
                 p("Download the predicted values alongside the true values 
                   for the selected response variable."),
                 downloadButton("download_predictions", "Download CSV")),
        tabPanel("Console Output", verbatimTextOutput("console_output"))
      )
    )
  )
)

server <- function(input, output, session) {
  data <- reactiveValues(train = NULL, test = NULL, preds = NULL, console_log = NULL)
  
  observe({
    data$train <- if (!is.null(input$train_file)) {
      read.csv(input$train_file$datapath)
    } else {
      read.csv("~/ShinyApps/optXpred/data/train_flux_24_2015_2017.csv")
    }
  })
  
  observe({
    data$test <- if (!is.null(input$test_file)) {
      read.csv(input$test_file$datapath)
    } else {
      read.csv("~/ShinyApps/optXpred/data/test_flux_24_2018_2024.csv")
    }
  })
  
  output$response_selector <- renderUI({
    req(data$train)
    selectInput("response_var", "Select Response Variable (Y):", choices = names(data$train), multiple = FALSE)
  })
  
  output$predictor_selector <- renderUI({
    req(data$train, input$response_var)
    selectInput("predictor_vars", "Select Predictor Variables (X):",
                choices = setdiff(names(data$train), input$response_var),
                multiple = TRUE)
  })
  
  observeEvent(input$run_model, {
    req(data$train, data$test, input$response_var, input$predictor_vars)
    
    train_mat <- as.matrix(data$train)
    test_mat <- as.matrix(data$test)
    all_vars <- colnames(train_mat)
    
    idx_Y <- match(input$response_var, all_vars)
    idx_X <- match(input$predictor_vars, all_vars)
    
    data$console_log <- capture.output({
      cat("Running opt_pred with response:", input$response_var, "\n")
      cat("Using predictors:", paste(input$predictor_vars, collapse = ", "), "\n")
      
      out <- opt_pred(train_mat[, idx_Y], train_mat[, idx_X], test_mat[, idx_X], 
                      prob = input$prop,
                      verbose = as.logical(input$verbose),
                      extrapolate_tail_CDF = as.logical(input$extrapolate_tail),
                      n_trees = input$n_trees)
      data$preds <- out$y_pred
    })
  })
  
  metrics_data <- reactive({
    req(data$preds, data$test, input$response_var)
    true_y <- data$test[[input$response_var]]
    pred_y <- data$preds
    
    obs <- as.integer(true_y > quantile(true_y, probs = 1 - input$event_rate, na.rm = TRUE))
    pred <- as.integer(pred_y > quantile(pred_y, probs = 1 - input$alarm_rate, na.rm = TRUE))
    
    compute_metrics(obs, pred)
  })
  
  output$metrics_table <- renderDT({
    req(metrics_data())
    m <- metrics_data()
    df <- data.frame(
      Metric = c("TSS", "Precision", "Recall", "F1 Score", "Miss Rate", "Alarm Rate", "TP", "FP", "FN", "TN"),
      Value = c(m$TSS, m$prec, m$recall, m$F1, m$missed, m$alarm, m$TP, m$FP, m$FN, m$TN)
    )
    datatable(df, options = list(dom = 't'))
  })
  
  output$console_output <- renderPrint({
    req(data$console_log)
    cat(paste(data$console_log, collapse = "\n"))
  })
  
  output$emp_precision <- renderPlotly({
    req(data$preds, data$test, input$response_var)
    true_y <- data$test[[input$response_var]]
    pred_y <- data$preds
    #p <- seq(from = input$p_range[1], to = input$p_range[2], length.out = 100)
    p <- seq(from = 0.75, to = 0.99, length.out = 100)
    out_tdep <- tdep(true_y, pred_y, p = p, se_return = TRUE)
    L <- out_tdep$L
    se <- out_tdep$se
    
    plot_ly() %>%
      add_lines(x = p, y = L, name = "Empirical Precision", line = list(color = "blue")) %>%
      add_lines(x = p, y = L + se * 1.96, name = "Upper 95%", line = list(dash = "dash",color="green")) %>%
      add_lines(x = p, y = L - se * 1.96, name = "Lower 95%", line = list(dash = "dash",color="green")) %>%
      layout(
        title = paste0("Tail Dependence for ", input$response_var),
        xaxis = list(title = "p"),
        yaxis = list(title = "Tail Dep (Precision for a Balanced Predictor)", range = c(0, 1))
      )
  })
  
  output$predictions_plot <- renderPlotly({
    req(data$preds, data$test, input$response_var)
    true_y <- (data$test[[input$response_var]])
    pred_y <- (data$preds)
    n_test <- length(true_y)
    obs_idx <- 1:n_test
    
    event_thresh <- quantile(true_y, probs = 1 - input$event_rate, na.rm = TRUE)
    alarm_thresh <- quantile(pred_y, probs = 1 - input$alarm_rate, na.rm = TRUE)
    
    event_flag <- true_y > event_thresh
    alarm_flag <- pred_y > alarm_thresh
    
    color_vec <- ifelse(event_flag & alarm_flag, "Correct",
                        ifelse(event_flag & !alarm_flag, "Missed", "Other"))
    
    plot_ly(
      x = obs_idx,
      y = true_y,
      type = "scatter",
      mode = "markers",
      color = color_vec,
      colors = c("Correct" = "blue", "Missed" = "red", "Other" = "black"),
      marker = list(size = 6)
    ) %>%
      layout(
        title = paste0("Predictions for ", input$response_var),
        xaxis = list(title = "Observation Index"),
        yaxis = list(
          title = input$response_var,
          type = ifelse(input$log_scale == "Log", "log", "linear")
        ),
        legend = list(title = list(text = "Classification"))
      )
  })
  output$scatter_plot <- renderPlotly({
    req(data$preds, data$test, input$response_var)
    
    true_y <- to_uniform_scale(data$test[[input$response_var]])
    pred_y <- to_uniform_scale(data$preds)
    
    event_thresh <- quantile(true_y, probs = 1 - input$event_rate, na.rm = TRUE)
    alarm_thresh <- quantile(pred_y, probs = 1 - input$alarm_rate, na.rm = TRUE)
    
    event_flag <- true_y > event_thresh
    alarm_flag <- pred_y > alarm_thresh
    
    # Classify points
    point_class <- ifelse(event_flag & alarm_flag, "TP",
                          ifelse(event_flag & !alarm_flag, "Missed",
                                 ifelse(!event_flag & alarm_flag, "False", "Other")))
    
    # Get indices for each class
    idx_TP <- which(point_class == "TP")
    idx_Missed <- which(point_class == "Missed")
    idx_False <- which(point_class == "False")
    idx_Other <- which(point_class == "Other")
    
    plt <- plot_ly()
    
    if (length(idx_TP) > 0) {
      plt <- plt %>% add_trace(
        x = pred_y[idx_TP],
        y = true_y[idx_TP],
        type = "scatter",
        mode = "markers",
        marker = list(color = "blue", size = 6),
        name = "True Positives (TP)"
      )
    }
    
    if (length(idx_Missed) > 0) {
      plt <- plt %>% add_trace(
        x = pred_y[idx_Missed],
        y = true_y[idx_Missed],
        type = "scatter",
        mode = "markers",
        marker = list(color = "red", size = 6),
        name = "Missed Events (FN)"
      )
    }
    
    if (length(idx_False) > 0) {
      plt <- plt %>% add_trace(
        x = pred_y[idx_False],
        y = true_y[idx_False],
        type = "scatter",
        mode = "markers",
        marker = list(color = "orange", size = 6),
        name = "False Alarms (FP)"
      )
    }
    
    if (length(idx_Other) > 0) {
      plt <- plt %>% add_trace(
        x = pred_y[idx_Other],
        y = true_y[idx_Other],
        type = "scatter",
        mode = "markers",
        marker = list(color = "black", size = 6),
        name = "True Negatives (TN)"
      )
    }
    
    # Horizontal event threshold line (Y)
    plt <- plt %>% add_lines(
      x = c(min(pred_y, na.rm = TRUE), max(pred_y, na.rm = TRUE)),
      y = c(event_thresh, event_thresh),
      line = list(color = "green", dash = "dash"),
      name = "Event Threshold (Y)"
    )
    
    # Vertical alarm threshold line (Predicted)
    plt <- plt %>% add_lines(
      x = c(alarm_thresh, alarm_thresh),
      y = c(min(true_y, na.rm = TRUE), max(true_y, na.rm = TRUE)),
      line = list(color = "purple", dash = "dash"),
      name = "Alarm Threshold (Pred)"
    )
    
    # Final layout
    plt <- plt %>%
      layout(
        title = paste0("True Y vs Predicted Y (", input$response_var, ")"),
        xaxis = list(title = "Predicted Y (Uniform Scale)"),
        yaxis = list(title = "True Y (Uniform Scale)"),
        legend = list(title = list(text = "Point Classification + Thresholds"))
      )
    
    plt
  })
  
  output$download_predictions <- downloadHandler(
    filename = function() {
      paste0("predictions_", input$response_var, ".csv")
    },
    content = function(file) {
      req(data$preds, data$test, input$response_var)
      df <- data.frame(
        True_Y = data$test[[input$response_var]],
        True_Y_Pareto_scale = to_pareto_scale(data$test[[input$response_var]]),
        Predicted_Y_Pareto_Scale = to_pareto_scale(data$preds)
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$suggest_vars, {
    req(data$train, input$response_var)
    
    all_vars_minus_response <- names(data$train);
    idx_response = which(all_vars_minus_response==input$response_var);
    all_vars_minus_response= all_vars_minus_response[-idx_response];
    
    keep = tdep_partial_select(y = data$train[[input$response_var]],
                               x = as.matrix(data$train[,-idx_response]), 
                               input$keep_threshold)
    
    suggested_vars = all_vars_minus_response[keep];
    
    # Update the predictor selection input
    updateSelectInput(session, "predictor_vars", selected = suggested_vars)
  })
  
  
}

shinyApp(ui, server)
