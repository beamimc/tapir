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
  return(downreg_exons)
  
  
}

get_nonreg_exons <- function(x_flat, downreg_exons){
  x_flat <- x_flat |>
    mutate(key = paste0(isoform, "-", exon_rank))
  right_keys <- paste0(downreg_exons$isoform, "-", downreg_exons$exon_rank+1)
  nonreg_exons <- x_flat |>
    filter(key %in% right_keys)
  
  # get all the exons not in downreg_exons
  non_spliced_exons <- x_flat |>
    filter(!(exon_id %in% downreg_exons$exon_id)) 
  
  
  #keep only the nearest exon for nonreg_exons 
  nonreg_auxid <-pair_nearest(downreg_exons, non_spliced_exons)$exon_id.y #only 22 unique -- many SE exons have the same nearest
  nonreg_exons <- x_flat |> filter(exon_id %in% nonreg_auxid) |> unique()
  return(nonreg_exons)
  
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
