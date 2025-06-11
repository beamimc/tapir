
plot_isoforms_wiggle <- function(exons, gene_of_interest, sig_res){
  sig_res <- sig_res %>%
    mutate(
      score = (1 - empirical_pval) * sign(estimates),  # Compute the score
      computed_color = custom_pal(score)               # Apply your custom palette
    )
  
  
  #tx from the gene present in the se
  # get the sig tx 
  sig_tx <- sig_res |> filter(symbol== gene_of_interest) |> select(isoform_id)

  exons_toshow <-  exons[sig_tx$isoform_id]
  
  ## create transcript_annotation for coloring. needs at least transcript_id and strand
  transcript_annotations <- tibble(
    transcript_id = names(exons_toshow),
    strand = purrr::map_chr(exons_toshow, function(gr) {
      strands <- as.character(strand(gr))
      unique_strands <- unique(strands)
      if (length(unique_strands) == 1) {
        unique_strands
      } else {
        "*"
      }
    })
  )
  
  transcript_annotations <- transcript_annotations %>%
    left_join(
      sig_res %>% select(isoform_id, computed_color),
      by = c("transcript_id" = "isoform_id")
    ) %>%
    mutate(
      color_by = ifelse(is.na(computed_color), "#F6BD60", computed_color)
    ) %>%
    select(-computed_color)  # optional: remove helper column
  
  p <- wiggleplotr::plotTranscripts(exons_toshow, 
                               transcript_annotations = transcript_annotations, 
                               rescale_introns = TRUE)
  p2 <- p +
    scale_x_continuous(
      expand = expansion(mult = c(0.35, 0.25))  # 5% padding each side
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0.1, 0.1))    # if you want vertical padding too
    )
  
  ggplotly(p2)
}

plot_downreg_exons <- function(exons, gene_of_interest, sig_res, downreg_exons){
  sig_res <- sig_res %>%
    mutate(
      score = (1 - empirical_pval) * sign(estimates),  # Compute coloring score
      computed_color = custom_pal(score)               # Apply custom palette
    )

  # get the sig tx 
  sig_tx <- sig_res |> filter(symbol== gene_of_interest) |> select(isoform_id)
  
  gene_id <-  sig_res |> filter(symbol== gene_of_interest) |> select(gene_id) 
  gene_id <- gene_id |> unique()
  
  exons_toshow <-  exons[sig_tx$isoform_id]
  
  
  filtered <- downreg_exons %>%
    filter(gene == gene_id[[1]])
  if (length(filtered) > 0) {
    exons_toshow[["downreg"]] <- filtered
  }
  
  # exons_toshow[["downreg"]] <- downreg_exons |> filter(gene == gene_of_interest)
  
  transcript_annotations <- tibble(
    transcript_id = names(exons_toshow),
    strand = purrr::map_chr(exons_toshow, function(gr) {
      strands <- as.character(strand(gr))
      unique_strands <- unique(strands)
      if (length(unique_strands) == 1) {
        unique_strands
      } else {
        "*"
      }
    })
  )
  
  
  transcript_annotations <- transcript_annotations %>%
    left_join(
      sig_res %>% select(isoform_id, computed_color),
      by = c("transcript_id" = "isoform_id")
    ) %>%
    mutate(
      color_by = ifelse(is.na(computed_color), "#F6BD60", computed_color)
    ) %>%
    select(-computed_color)  # optional: remove helper column
  
  p <- wiggleplotr::plotTranscripts(exons_toshow, 
                               transcript_annotations = transcript_annotations, 
                               rescale_introns = TRUE)
  
  p2 <- p +
    scale_x_continuous(
      expand = expansion(mult = c(0.35, 0.25))  # 5% padding each side
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0.1, 0.1))    # if you want vertical padding too
    )
  
  ggplotly(p2)  
}
