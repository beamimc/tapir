
plot_isoforms_wiggle <- function(gene_of_interest){
  sig_res <- sig_res %>%
    mutate(
      score = (1 - empirical_pval) * sign(estimates),  # Compute the score
      computed_color = custom_pal(score)               # Apply your custom palette
    )
  
  
  #tx from the gene present in the se
  # get the sig tx 
  sig_tx <- sig_res |> filter(symbol== gene_of_interest) |> select(isoform_id)
  
  exons_toshow <-  exons[sig_tx$isoform_id]
  
  # make a named vector of signs
  signs <- setNames(sig_res$sign, sig_res$isoform_id)
  # pull out the signs in the same order as our GRangesList names
  the_signs <- signs[names(exons_toshow)]
  # build the new names
  new_names <- ifelse(
    !is.na(the_signs) & the_signs == 1,
    paste0("+", names(exons_toshow)),
    paste0("-", names(exons_toshow))
  )
  # assign back
  names(exons_toshow) <- new_names
  
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
    mutate(
      transcript_id_aux = sub("^[+-]", "", transcript_id)  # Clean for joining
    )
  
  transcript_annotations <- transcript_annotations %>%
    left_join(
      sig_res %>% select(isoform_id, computed_color),
      by = c("transcript_id_aux" = "isoform_id")
    ) %>%
    mutate(
      color_by = ifelse(is.na(computed_color), "green", computed_color)
    ) %>%
    select(-computed_color)  # optional: remove helper column
  
  p <- wiggleplotr::plotTranscripts(exons_toshow, 
                               transcript_annotations = transcript_annotations, 
                               rescale_introns = TRUE)
  ggplotly(p)
  
  
}
plot_downreg_exons <- function(gene_of_interest){
  
  
  
  #tx from the gene present in the se
  present_tx <-  names(se)[ rowData(se)$gene_id ==gene_of_interest ]
  # get the sig tx 
  sig_tx <- sig_res |> filter(gene_id== gene_of_interest) |> select(isoform_id)
  
  exons_toshow <-  exons[sig_tx$isoform_id]
  
  # make a named vector of signs
  signs <- setNames(sig_res$sign, sig_res$isoform_id)
  # pull out the signs in the same order as our GRangesList names
  the_signs <- signs[names(exons_toshow)]
  # build the new names
  new_names <- ifelse(
    !is.na(the_signs) & the_signs == 1,
    paste0("+", names(exons_toshow)),
    paste0("-", names(exons_toshow))
  )
  # assign back
  names(exons_toshow) <- new_names
  
  exons_toshow[["downreg"]] <- downreg_exons |> filter(gene == gene_of_interest)
  
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
    mutate(
      transcript_id_aux = sub("^[+-]", "", transcript_id)  # Clean for joining
    )
  
  transcript_annotations <- transcript_annotations %>%
    left_join(
      sig_res %>% select(isoform_id, computed_color),
      by = c("transcript_id_aux" = "isoform_id")
    ) %>%
    mutate(
      color_by = ifelse(is.na(computed_color), "green", computed_color)
    ) %>%
    select(-computed_color)  # optional: remove helper column
  
  wiggleplotr::plotTranscripts(exons_toshow, 
                               transcript_annotations = transcript_annotations, 
                               rescale_introns = TRUE)
  
  
}
