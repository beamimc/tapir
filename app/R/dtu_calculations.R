
my_theme <- bs_theme(version = 5, bootswatch = "yeti")
width_upstream <- 100

# # Palette for negative values: from red (-1) to gray (0)
# pal_neg <- col_numeric(
#   palette = c("#FF2900", "#ffd4cc", "#E0E0E0"),
#   domain = c(-1, .8)
# )
# 
# # Palette for positive values: from gray (0) to blue (1)
# pal_pos <- col_numeric(
#   palette = c("#E0E0E0", "#cae3ff", "#4FA4FF"),
#   domain = c(0.8, 1)
# )
# 

# Palette for negative values: from red (-1) to gray (0)
pal_neg <- col_numeric(
  palette = c("#FC6C85", "#E0E0E0"),
  domain = c(-1, .8)
)

# Palette for positive values: from gray (0) to blue (1)
pal_pos <- col_numeric(
  palette = c("#E0E0E0",  "#5DADEC"),
  domain = c(0.8, 1)
)


# Custom palette function that applies the appropriate palette based on value
custom_pal <- function(x) {
  sapply(x, function(val) {
    if (val < 0) {
      pal_neg(val)
    } else {
      pal_pos(val)
    }
  })
}


get_dtu_column <- function(se, cd1, cd2) {
  rd_names <- names(SummarizedExperiment::rowData(se))
  
  direct <- paste0("fitDTUResult_", cd2, "_vs_", cd1) #exp_vs_ctl (cd1==ref)
  reverse <- paste0("fitDTUResult_", cd1, "_vs_", cd2)
  
  if (direct %in% rd_names) {
    return(list(column_name = direct, direction = 1))
  }
  
  if (reverse %in% rd_names) {
    return(list(column_name = reverse, direction = -1))
  }
  
  stop(
    sprintf("No DTU result found for '%s vs %s' or '%s vs %s'", cd2, cd1, cd1, cd2)
  )
}


get_sig_res <- function(se, fdr_threshold, cd1, cd2){
  dtu_column <- get_dtu_column(se, cd1, cd2)
  column <- dtu_column$column_name
  dtu_direction <- dtu_column$direction
  
  sig_res <- rowData(se)[[column]] |>
    tibble::as_tibble() |>
    dplyr::bind_cols(as.data.frame(rowData(se)[,1:4])) |>
    dplyr::filter(empirical_FDR < fdr_threshold) |>
    dplyr::select(gene_id, isoform_id, symbol, estimates, empirical_pval, empirical_FDR) |>
    dplyr::arrange(empirical_pval)
  
  sig_res <-  sig_res %>%
    dplyr::mutate(sign = sign(estimates),
                  dtu_column = column,
                  dtu_direction,
                  cd1 = cd1,
                  cd2 = cd2)
  
  return(sig_res)
  
}


get_x_flat <- function(exons, sig_res){
  sig_exons <- exons[names(exons) %in% sig_res$isoform_id] #get GRangesList only from the DTUs 61 - 35 genes
  #62 transcripts GRangesList - 1 duplicate `ENSG00000198467.13-305f0cb0` 
  #remove duplicate 
  name_counts <- table(names(sig_exons))
  dup_names <- names(name_counts)[name_counts > 1]
  sig_exons <- sig_exons[! duplicated(names(sig_exons))]
  
  
  # set if exons and internal or boundary 
  sig_exons@unlistData$internal <- TRUE
  sig_exons@unlistData$internal[start(sig_exons@partitioning)] <- FALSE
  sig_exons@unlistData$internal[end(sig_exons@partitioning)] <- FALSE
  
  flat_sig_exons <- unlist(sig_exons)
  
  #include coef +/- column from the DTU analysis saturn 
  flat_sig_exons$coef <- sig_res$estimates[match(names(flat_sig_exons), sig_res$isoform_id)]
  flat_sig_exons$sign <- sig_res$sign[match(names(flat_sig_exons), sig_res$isoform_id)]
  
  #include gene name for each transcript name 
  flat_sig_exons$gene <- sig_res$gene_id[match(names(flat_sig_exons), sig_res$isoform_id)]
  flat_sig_exons$isoform <- names(flat_sig_exons)
  
  x_flat <- flat_sig_exons
  return(x_flat)
  
}
get_dtu_df <- function(sig_res){
  # Create a data frame to display Transcript, Gene Symbol, and P-value in the UI.
  significant_transcripts <- sig_res$isoform_id
  dtu_df <- data.frame(
    Transcript = significant_transcripts,
    Symbol = sig_res$symbol,
    P_Value = format(
      sig_res$empirical_pval,
      scientific = TRUE,
      digits     = 2     # one significant digit → “3e-03”
    ),
    stringsAsFactors = FALSE
  )
  return(dtu_df)
  
}

# Functions: Counts, DTU and p-values
# ------------------------------------------------------------------------------

# Calculate the proportion of counts per transcript for a given gene symbol
calc_prop <- function(se, symbol, sig_res) {
  cts <- assay(se, "counts")[mcols(se)$symbol == symbol, ]
  prop <- t(cts) / colSums(cts)
  
  sig_ids <- sig_res %>%
    filter(symbol == symbol) %>%
    pull(isoform_id)
  
  prop <- prop[, colnames(prop) %in% sig_ids, drop = FALSE]
  return(prop)
}

# Calculate mean differences for Differential Transcript Usage (DTU)
calc_mean_diff_DTU <- function(se, gene_symbol, sig_res) {
  cd1 <- sig_res$cd1|>unique()
  cd2 <- sig_res$cd2|>unique()
  
  cts <- assay(se, "counts")[mcols(se)$symbol == gene_symbol, ]
  prop <- t(cts) / colSums(cts)
  
  # Compute difference in mean proportions between two conditions for each transcript
  n <- length(colnames(prop))
  mean_diffs_DTU <- sapply(1:n, function(j) {
    mean(prop[se$condition == cd2, j]) - mean(prop[se$condition == cd1, j])
  })
  names(mean_diffs_DTU) <- rownames(cts)  # Use transcript IDs as names
  
  sig_ids <- sig_res %>%
    filter(symbol == gene_symbol) %>%
    pull(isoform_id)
  
  mean_diffs_DTU <- mean_diffs_DTU[ sig_ids ]
  return(mean_diffs_DTU)
}

# Extract p-values from the Saturn DTU analysis for a given gene symbol
get_pvals <- function(se, gene_symbol, sig_res) {
  # Identify the column in rowData that contains DTU results
  dtu_column <- sig_res$dtu_column |> unique()
  
  # Subset rowData for matching symbol and extract p-values
  pvals <- rowData(se[rowData(se)$symbol == gene_symbol, ])[[dtu_column]]$empirical_pval
  names(pvals) <- rownames(rowData(se[rowData(se)$symbol == gene_symbol, ]))
  sig_ids <- sig_res %>%
    filter(symbol == gene_symbol) %>%
    pull(isoform_id)
  
  pvals <- pvals[ sig_ids ]
  return(pvals)
}

parse_saturnDTU_conditions <- function(se) {
  rd_names <- names(SummarizedExperiment::rowData(se))
  dtu_cols <- grep("^fitDTUResult_", rd_names, value = TRUE)
  
  if (length(dtu_cols) == 0) return(tibble::tibble())
  
  condition_df <- stringr::str_match(dtu_cols, "^fitDTUResult_(.+)_vs_(.+)$")
  
  valid <- complete.cases(condition_df)
  
  tibble::tibble(
    column_name = dtu_cols[valid],
    cd2 = condition_df[valid, 2],
    cd1 = condition_df[valid, 3]
  )
}