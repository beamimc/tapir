
plot_isoforms_wiggle <- function(gene_of_interest, sig_res){
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
  
  cat("Names of exons_toshow:", names(exons_toshow), "\n")
  cat("the_signs:", the_signs, "\n")
  
  
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

plot_downreg_exons <- function(gene_of_interest, sig_res, downreg_exons){
  sig_res <- sig_res %>%
    mutate(
      score = (1 - empirical_pval) * sign(estimates),  # Compute the score
      computed_color = custom_pal(score)               # Apply your custom palette
    )
  
  #tx from the gene present in the se
  # get the sig tx 
  sig_tx <- sig_res |> filter(symbol== gene_of_interest) |> select(isoform_id)
  
  gene_id <-  sig_res |> filter(symbol== gene_of_interest) |> select(gene_id) 
  gene_id <- gene_id |> unique()
  
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
    mutate(
      transcript_id_aux = sub("^[+-]", "", transcript_id)  # Clean for joining
    )
  
  transcript_annotations <- transcript_annotations %>%
    left_join(
      sig_res %>% select(isoform_id, computed_color),
      by = c("transcript_id_aux" = "isoform_id")
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
