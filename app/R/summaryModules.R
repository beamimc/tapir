summaryStatsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 4,
        tags$h4("Genes with DTU isoforms"),
        checkboxInput(ns("select_all"), "Select All Genes", value = FALSE),
        DTOutput(ns("gene_table")),
        br(),
        actionButton(ns("run_go"), "Run GO", class = "btn btn-primary mt-2")
      ),
      
      column(
        width = 8,
        tags$div(
          class = "card mb-3",
          tags$div(
            class = "card-header bg-primary text-white", 
            "Gene Ontology"
          ),
          tags$div(
            class = "card-body",
            withSpinner(
              plotOutput(ns("go_plot"), width = "100%"),
              type = 4,
              color = "#0d6efd"
            )
          )
        )
      )
    )
  )
}

summaryStatsServer <- function(id, dtu_df, sig_res) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive gene table
    gene_df <- reactive({
      dtu_df() |> dplyr::distinct(Symbol)
    })
    
    # Reactive: genes selected by user
    selected_gene_list <- reactive({
      genes <- gene_df()
      
      # compute selection indices
      sel <- if (isTRUE(input$select_all)) {
        seq_len(nrow(genes))  # all row indices
      } else {
        input$gene_table_rows_selected
      }
      
      # return selected gene symbols
      if (!is.null(sel) && length(sel) > 0) {
        genes$Symbol[sel]
      } else {
        genes$Symbol[1]  # fallback
      }
    })
    
    
    # Output: DT table
    output$gene_table <- renderDT({
      datatable(
        gene_df(),
        selection = "multiple",
        options = list(pageLength = 5)
      )
    })

    
    proxy <- dataTableProxy(ns("gene_table"))
    
    observeEvent(input$select_all, {
      genes <- gene_df()
      all_rows <- seq_len(nrow(genes))
      
      if (isTRUE(input$select_all)) {
        selectRows(proxy, all_rows)
      } else {
        selectRows(proxy, NULL)
      }
    })
    
    proxy <- dataTableProxy("gene_table", session = session)
    
    observeEvent(input$select_all, {
      rows <- if (isTRUE(input$select_all)) {
        seq_len(nrow(gene_df()))   # all rows
      } else {
        NULL                       # clear
      }
      selectRows(proxy, rows)     
    })
    
    observe({
      cat("Final selected genes:\n")
      print(selected_gene_list())
    })
    # Store selected genes on button press
    stored_genes <- reactiveVal(NULL)
    observeEvent(input$run_go, {
      req(selected_gene_list())
      stored_genes(selected_gene_list())
    })
    
    # Render GO plot
    output$go_plot <- renderPlot({
      req(stored_genes())
      go_plot(stored_genes())
    })
  })
}
