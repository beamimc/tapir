
exonLevelUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      class = "border-bottom p-3",   
      
      column(width = 6,
             tags$h4("Significant DTU transcripts"),
             DTOutput(ns("dtu_table"))
      ),
      column(width = 6,
             tags$div(class = "card mb-3",
                      tags$div(class = "card-header bg-primary text-white", 
                               "Isoform Structures with selected exons"),
                      tags$div(class = "card-body",
                               plotlyOutput(ns("exon_level_plot"), width = "100%"))
                      )
              )
    ),
    # Tabs for downstream windows and another plot
    tabsetPanel(
      id = ns("exon_tabs"),  # optional id
      tabPanel(
        "Window Plots",
        fluidRow(
          column(
            width = 6,
            tags$div(class = "card mb-3",
                     tags$div(class = "card-header bg-primary text-white", 
                              "Up and downstream sequences"),
                     tags$div(class = "card-body",
                              plotOutput(ns("window_summary_plot_down"), width = "100%"))
            )
          ),
          column(
            width = 6,
            tags$div(class = "card mb-3",
                     tags$div(class = "card-header bg-primary text-white", 
                              "Up and downstream sequences"),
                     tags$div(class = "card-body",
                              plotOutput(ns("window_summary_plot_nonreg"), width = "100%"))
            )
          )
        )
      ),
      tabPanel(
        "Single Plot",
        fluidRow(
          class = "border-bottom p-3",   
          column(
            width = 12,
            tags$div(class = "card mb-3",
                     tags$div(class = "card-header bg-primary text-white", 
                              "Alternative View"),
                     tags$div(class = "card-body",
                              plotOutput(ns("plot_window_comparison"), width = "100%"))
            )
          )
        )
      )
    )
  )
}
    
    
exonLevelServer <- function(id, exons, dtu_df, x_flat, sig_res) {
  moduleServer(id, function(input, output, session) {

    
    selected_gene <- reactive({
      sel <- input$dtu_table_rows_selected
      if (!is.null(sel) && length(sel) == 1) {
        dtu_df()[["Symbol"]][sel]
      } else {
        dtu_df()[["Symbol"]][1]
      }
    })
    
    output$dtu_table <- renderDT({
      datatable(dtu_df(), selection = "single", options = list(pageLength = 5))
    })
    
    # Reactive: Downregulated exons and windows
    downstream_data <- reactive({
      req(x_flat())
      downreg_exons <- detect_downreg_exonsv2(x_flat())
      validate(need(length(downreg_exons) > 0, "No downregulated exons found"))
      list(
        exons = downreg_exons,
        downstream = get_downstream_from_GRanges(downreg_exons),
        upstream = get_upstream_from_GRanges(downreg_exons)
      )
    })
    
    # Reactive: Nonregulated exons and windows
    nonreg_data <- reactive({
      req(x_flat())
      downreg_exons <- detect_downreg_exonsv2(x_flat())
      nonreg_exons <- get_nonreg_exons(x_flat(), downreg_exons)
      validate(need(length(nonreg_exons) > 0, "No nonregulated exons found"))
      list(
        exons = nonreg_exons,
        downstream = get_downstream_from_GRanges(nonreg_exons),
        upstream = get_upstream_from_GRanges(nonreg_exons)
      )
    })
    
    output$exon_level_plot <- renderPlotly({
      req(selected_gene(), downstream_data())
      plot_downreg_exons(exons, selected_gene(), sig_res(), downstream_data()$exons)
    })
    
    output$window_summary_plot_down <- renderPlot({
      req(downstream_data())
      plot_updownstream_windows(
        downstream_data()$upstream,
        downstream_data()$downstream,
        exon_label = "downreg exon"
      )
    })
    
    output$window_summary_plot_nonreg <- renderPlot({
      req(nonreg_data())
      plot_updownstream_windows(
        nonreg_data()$upstream,
        nonreg_data()$downstream,
        exon_label = "nonreg exon"
      )
    })
    
    output$plot_window_comparison <- renderPlot({
      req(nonreg_data())
      req(downstream_data())
      plot_window_comparison(
        downstream_data()$upstream,
        nonreg_data()$upstream
      )
    })
    
  })
}
