server <- function(input, output, session) {

  
  iso_mod <- isoformAnalysisServer("isoform", gene_ids, dtu_df)
  exonLevelServer("exon", exon_ids)
  summaryStatsServer("summary")
  
}


