server <- function(input, output, session) {
  # … your view‐saving logic here …
  iso_mod <- isoformAnalysisServer("isoform", gene_ids, dtu_df)
  exonLevelServer("exon", exon_ids)
  summaryStatsServer("summary")
  
  # sync global_genes ↔ iso_mod$selected_gene() as before…
}


