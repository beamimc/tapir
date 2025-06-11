server <- function(input, output, session, data) {
  
  se <- data$se
  exons <- data$exons
  condition_choices_df <- data$condition_choices
  
  applied_fdr <- reactiveVal(0.05)
  observeEvent(input$apply_fdr, {
    applied_fdr(input$fdr_threshold)  # apply only when button is clicked
  })
  
  
  # Initialize cd1 choices
  observe({
    if (nrow(condition_choices_df) > 0) {
      all_cd1 <- unique(condition_choices_df$cd1)
      updateRadioButtons(session, "cd1", choices = all_cd1, selected = all_cd1[1])
    }
  })
  
  
  # Dynamically update cd2 based on selected cd1
  observeEvent(input$cd1, {
    req(input$cd1)
    valid_cd2 <- condition_choices_df |> 
      filter(cd1 == input$cd1) |> 
      pull(cd2) |> 
      unique()
    
    updateRadioButtons(session, "cd2", choices = valid_cd2, selected = valid_cd2[1])
  })
  
  # On Apply: check if selected cd1/cd2 form a valid pair initialize with first row
  # Set initial default contrast
  default_cd1 <- condition_choices_df$cd1[1]
  default_cd2 <- condition_choices_df$cd2[1]
  
  selected_conditions <- reactiveVal(
    list(cd1 = default_cd1, cd2 = default_cd2)
  )
  
  observeEvent(input$apply_pair, {
    req(input$cd1, input$cd2)
    valid <- any(condition_choices_df$cd1 == input$cd1 & 
                   condition_choices_df$cd2 == input$cd2)
    
    if (valid) {
      selected_conditions(list(cd1 = input$cd1, cd2 = input$cd2))
      cat("Applied pair:", input$cd1, "vs", input$cd2, "\n")
    } else {
      showNotification("Invalid contrast pair selected", type = "error")
    }
  })
  output$current_comparison <- renderText({
    conds <- selected_conditions()
    req(conds$cd1, conds$cd2)
    paste(conds$cd2, "/", conds$cd1)
  })
  
  
  
  
  # Define sig_res reactively based on the applied FDR
  sig_res <- reactive({
    conds <- selected_conditions()
    req(conds$cd1, conds$cd2)
    get_sig_res(se, applied_fdr(), conds$cd1, conds$cd2)
  })
  
  filtered_dtu_df <- reactive({
    dtu_df <- get_dtu_df(sig_res())  # use the cached reactive value
    dtu_df
  })
  
  x_flat <- reactive({
    x_flat <- get_x_flat(exons, sig_res())  # reuse same reactive sig_res()
    x_flat
  })
  # 
  # gene_ids <- reactive({
  #   sort(unique(sig_res()$symbol))
  # })
  # 
  
  isoformAnalysisServer("isoform", se, exons, filtered_dtu_df, sig_res)
  exonLevelServer("exon",exons, filtered_dtu_df, x_flat, sig_res)
  summaryStatsServer("summary", filtered_dtu_df, sig_res)
  
}


