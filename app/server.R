server <- function(input, output, session) {
  applied_fdr <- reactiveVal(0.05)
  observeEvent(input$apply_fdr, {
    applied_fdr(input$fdr_threshold)  # apply only when button is clicked
  })
  
  # Define sig_res reactively based on the applied FDR
  sig_res <- reactive({
    get_sig_res(applied_fdr())
  })
  
  filtered_dtu_df <- reactive({
    dtu_df <- get_dtu_df(sig_res())  # use the cached reactive value
    dtu_df
  })
  
  x_flat <- reactive({
    x_flat <- get_x_flat(sig_res())  # reuse same reactive sig_res()
    x_flat
  })
  
  gene_ids <- reactive({
    sort(unique(sig_res()$symbol))
  })
  
  
  isoformAnalysisServer("isoform", gene_ids, filtered_dtu_df, sig_res)
  exonLevelServer("exon", filtered_dtu_df, x_flat, sig_res)
  summaryStatsServer("summary", filtered_dtu_df, sig_res)
  
}


