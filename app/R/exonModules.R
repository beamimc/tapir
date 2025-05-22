
exonLevelUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6,
             tags$h3("Significant DTU transcripts"),
             DTOutput(ns("exon_gene_table"))
      ),
      column(6,
             plotOutput(ns("window_summary_plot"))
      )
    ),
    fluidRow(
      column(6,
             plotOutput(ns("exon_level_plot"))   # <-- new plotOutput
      )
    )
  )
}
detect_downreg_exonsv2 <- function(x_flat) {
  # x_flat <- flat_sig_exons
  x_flat <- x_flat |>
    mutate(key = paste0(isoform, "-", exon_rank))
  # 1) split + vs – strand
  plus_exons  <- x_flat %>% filter(sign ==  1)
  neg_exons <- x_flat %>% filter(sign ==  -1)
  
  candidates <- plus_exons %>%
    filter_by_non_overlaps_directed(neg_exons) %>%
    mutate(SE = TRUE) |>
    filter(internal == TRUE)
  
  left_keys <- paste0(candidates$isoform, "-", candidates$exon_rank-1)
  left_exons <- x_flat |>
    filter(key %in% left_keys)
  
  right_keys <- paste0(candidates$isoform, "-", candidates$exon_rank+1)
  right_exons <- x_flat |>
    filter(key %in% right_keys)
  
  candidates <-  candidates |>
    mutate(left_and_right =
             left_exons %in% neg_exons &
             right_exons %in% neg_exons
    )
  downreg_exons  <- candidates |> filter(left_and_right == TRUE)
  
  
}

get_downstream_from_GRanges <- function(GRanges,
                                        width_upstream=100
){
  
  upstr_exons <- GRanges %>%
    flank_downstream(width = width_upstream) 
  # The result will be another GRanges object that still contains 158 ranges,
  # but each range now represents the upstream flanking region of the corresponding exon. 
  
  df <- get_sliding_windows(upstr_exons)
  
  return(df)
}


get_upstream_from_GRanges <- function(GRanges,
                                        width_upstream=100
){
  
  upstr_exons <- GRanges %>%
    flank_upstream(width = width_upstream) 
  # The result will be another GRanges object that still contains 158 ranges,
  # but each range now represents the upstream flanking region of the corresponding exon. 
  
  df <- get_sliding_windows(upstr_exons)
  
  return(df)
}


exonLevelServer <- function(id, gene_data) {
  moduleServer(id, function(input, output, session) {
    output$exon_gene_table <- renderDT({
      datatable(gene_data, selection = "single", options = list(pageLength = 5))
    })
    
    output$exon_level_plot <- renderPlot({
      sel <- input$exon_gene_table_rows_selected
      req(sel)
      gene <- gene_data$Symbol[sel]
      exon_level_plot(gene)
    })
    

      downreg_exons <- detect_downreg_exonsv2(x_flat)
      downreg_exons_wind_downstream <- get_downstream_from_GRanges(downreg_exons)
      downreg_exons_wind_upstream <- get_upstream_from_GRanges(downreg_exons)
      
      output$window_summary_plot <- renderPlot({
        plot_updownstream_windows(downreg_exons_wind_upstream, downreg_exons_wind_downstream,
                            exon_label="downreg exon")
        
    })
  })
}
