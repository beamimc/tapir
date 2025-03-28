server <- function(input, output, session) {
  
  # --- Initialize selected_gene() reactiveVal
  selected_gene <- reactiveVal(NULL)
  
  # --- Populate selectizeInput choices (assumes gene_ids exists globally or is passed in)
  updateSelectizeInput(session, "selected_gene", 
                       choices = gene_ids, 
                       server = TRUE)
  
  # --- When the selectizeInput changes, update selected_gene()
  observeEvent(input$selected_gene, {
    selected_gene(input$selected_gene)
  })
  
  # --- When table row is selected, update selected_gene() and sync selectizeInput
  observeEvent(input$gene_table_rows_selected, {
    selected_row <- input$gene_table_rows_selected
    if (!is.null(selected_row)) {
      gene <- gene_data$Symbol[selected_row]
      selected_gene(gene)
    }
  })
  
  # --- Reactive expressions now depend on selected_gene() instead of input$selected_gene
  mean_diffs_DTU <- reactive({
    req(selected_gene())
    calc_mean_diff_DTU(selected_gene())
  })
  
  prop <- reactive({
    req(selected_gene())
    calc_prop(selected_gene())
  })
  
  pvals <- reactive({
    req(selected_gene())
    get_pvals(selected_gene())
  })
  
  plot_path <- reactive({
    req(selected_gene())
    temp_file <- "www/temp.png"
    plot_gene_txs(selected_gene(), temp_file, mean_diffs_DTU(), pvals())  # Generate plot
    return(temp_file)
  })
  
  # --- Update gene description when gene is selected
  output$gene_description <- renderText({
    req(selected_gene())
    gene_info <- get_description(selected_gene())
    if (length(gene_info) > 0) {
      return(gene_info)
    } else {
      return("Select a gene to see description.")
    }
  })
  
  # --- Reactive table for gene selection
  output$gene_table <- renderDT({
    datatable(gene_data, selection = "single", 
              options = list(pageLength = 5))
  })
  
  # --- GO table (recalculate when gene changes)
  output$go_table <- renderDT({
    req(selected_gene())
    go_data <- get_GO(selected_gene())
    datatable(go_data, 
              options = list(pageLength = 5))
  })
  
  # --- Render main gene plot image
  output$gene_plot <- renderImage({
    req(plot_path())
    list(
      src = plot_path(),
      contentType = "image/png",
      width = "100%",
      height = "auto"
    )
  }, deleteFile = FALSE)
  
  # --- Render barplot
  output$barplot <- renderPlotly({
    req(mean_diffs_DTU(), pvals())
    barplot_meandifs(mean_diffs_DTU(), pvals())
  })
  
  # --- Render lineplot
  output$lineplot <- renderPlotly({
    req(prop())
    line_plot_txp_comparison(prop())
  })
}
