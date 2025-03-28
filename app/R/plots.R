

# Function: Plot isoforms per gene using wiggleplotr
# ---------------------------------------------------------------
plot_gene_txs_wiggleplotr <- function(gene) {
  
  # Extract transcript IDs associated with the selected gene
  tx_to_show <- txp %>%
    filter(symbol == gene) %>%
    as_tibble() %>%
    dplyr::pull(txid)
  
  # Extract exon structures for the selected transcripts
  exons_list <- exons_per_txp_GRL[tx_to_show]
  
  # Generate transcript plot
  wiggleplotr::plotTranscripts(exons_list, rescale_introns = FALSE)
}


# Function: Plot isoforms per gene using plotgardener
# ---------------------------------------------------------------
plot_gene_txs <- function(gene_symbol , output_file, mean_diffs_DTU, pvals) {
  
  # Identify the gene of interest by its geneid given its symbol (selectizeInput)
  gene_of_interest <- unique(rowData(se)[rowData(se)$symbol == gene_symbol, "gene_id"])

  
  # Extract genomic coordinates for the gene
  gene_coords <- genes(txdb, filter = list(GENEID = gene_of_interest))
  transcripts <- transcripts(txdb, filter = list(GENEID = gene_of_interest))
  
  # filter to keep only the transcripts in the experiment
  assay_transcripts <-  transcripts |> subset(tx_name %in% names(mean_diffs_DTU)) 
  
  # Extract chromosome, start, and end positions of the gene
  gene_chrom <- as.character(seqnames(gene_coords))
  gene_start <- min(start(gene_coords))
  gene_end <- max(end(gene_coords))
  
  # Define plotting parameters
  par <- pgParams(
    chrom = gene_chrom, 
    chromstart = gene_start, chromend = gene_end,
    assembly = db_assembly, just = c("left", "bottom")
  )
  

  
  #color all txp corresponding to the gene
  df_values <- data.frame(
    transcript = intersect(names(mean_diffs_DTU), names(pvals)),
    mean_diff = mean_diffs_DTU[intersect(names(mean_diffs_DTU), names(pvals))],
    pval = pvals[intersect(names(mean_diffs_DTU), names(pvals))],
    stringsAsFactors = FALSE
  ) %>%
    mutate(score = (1 - pval) * sign(mean_diff),
           color = custom_pal(score))
  
  hilite <- df_values %>% select(transcript, color)
  
  # Dynamically adjust plot dimensions based on transcript count
  plotT_height <- round(length(assay_transcripts) / 3, 3)  # Height scales with transcript count
  page_height <- plotT_height + 1  
  page_width <- 10 
  output_file <-file.path(wd, "app/www/temp.png")
  # Save the plot as a PNG file
  png(output_file, width = page_width * 200, height = page_height * 200, res = 170)
  
  # Create a new plotting page
  pageCreate(
    width = page_width + 1, height = page_height, 
    default.units = "inches", showGuides = FALSE
  )  
  
  # Plot transcripts with dynamic height 
  plotgardener::plotTranscripts(
    params = par, 
    x = 0.5, y = 0, width = page_width, height = plotT_height, 
    just = c("left", "top"), default.units = "inches",
    transcriptHighlights = hilite,
    transcriptFilter = assay_transcripts$tx_name,
    limitLabel = FALSE,
    fontsize = 9,
    labels = "transcript"
  )
  
  # Add genome coordinate labels to the plot
  plotGenomeLabel(
    chrom = gene_chrom, 
    chromstart = gene_start, chromend = gene_end,
    x = 0.5, y = page_height - 0.5,  
    just = c("left", "top"),
    length = page_width
  ) 
  
  dev.off() 
}


# Function: create df_long for ggplot functions given 2 conditions 
# ---------------------------------------------------------------

get_plot_data <- function(prop,
                           cd1 = "ctrl", 
                           cd2 = "exp"
                           ) {
  
  # Subset only samples from the two selected conditions
  subset_indices <- colData(se)$condition %in% c(cd1, cd2)
  prop_subset <- prop[subset_indices, ]
  condition_subset <- se$condition[subset_indices]
  
  # Convert to long format for ggplot
  plot_data <- melt(prop_subset)
  colnames(plot_data) <- c("Sample", "Transcript", "Proportion")
  # add conditions
  sample_metadata <- as.data.frame(colData(se)[, c("sample_id", "condition")])
  plot_data <- merge(plot_data, sample_metadata, by.x = "Sample", by.y = "sample_id") #inner joinis
  plot_data |>
    rename(Condition = condition)
  return(plot_data)
  
}

# Function: Boxplots count comparion KD vs WT (change hardcode)
# ---------------------------------------------------------------

boxplot_count_comparison <- function(prop, cd1 = "ctrl", cd2 = "exp") {
  
  plot_data <- get_plot_data(prop)
  
  # Filter the color palette to only include the selected conditions
  selected_colors <- condition_colors[c(cd1, cd2)]
  
  # Create the boxplot with dynamic color mapping
  ggplot(plot_data, aes(x = Transcript, y = Proportion, fill = condition)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
    theme_classic() +
    labs(x = "Transcript", y = "Proportion", title = "Proportion of Transcripts by Condition") +
    scale_fill_manual(values = selected_colors) +
    theme(axis.text.x = element_text(angle = 0))  # Rotate x-axis labels if needed
}

# Function: spaguetti plot for transcript comparion KD vs WT (change hardcode)
# ---------------------------------------------------------------

line_plot_txp_comparison  <- function(prop, cd1 = "ctrl", cd2 = "exp"){
  
  plot_data <- get_plot_data(prop)
  
  summary_data <- plot_data |>
    group_by(Transcript, condition) |>
    summarise(
      mean_proportion = mean(Proportion),
      sd_proportion = sd(Proportion),
      n = dplyr::n(),
      se_proportion = sd_proportion / sqrt(n)  # Standard error
    ) |>
    # some names are too long - issues to fix
    mutate(Transcript_short = substr(Transcript, 1, 10))
  
  # Filter the color palette to only include the selected conditions
  selected_colors <- condition_colors[c(cd1, cd2)]
  
  p <- ggplot(summary_data, aes(x = condition, y = mean_proportion,
                                group = Transcript, color = Transcript_short)) +
    geom_line() +
    geom_point() + 
    geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                      ymax = mean_proportion + se_proportion), 
                  width = 0.2) +           # Error bars for each point
    
    theme_classic() +
    labs(x = "Condition", y = "Proportion", title = "Proportion of Transcripts by Condition") +
    theme(axis.text.x = element_text(angle = 0))
  ggplotly(p)
}


# Function: barplot of mean difs KD vs WT (change hardcode)
# ---------------------------------------------------------------
barplot_meandifs <- function(mean_diffs_DTU, pvals, cd1 = "ctrl", cd2 = "exp"){
  
  # Extract names
  names_mean_diffs <- names(mean_diffs_DTU)
  names_pvals <- names(pvals)
  
  # Ensure names are aligned properly
  common_transcripts <- intersect(names_mean_diffs, names_pvals)
   
  plot_df <- data.frame(
    Transcript = common_transcripts,
    MeanDiff = mean_diffs_DTU[common_transcripts],
    pval = pvals[common_transcripts],
    stringsAsFactors = FALSE
  )
  
  
  # Define significance threshold
  plot_df <- plot_df %>%
    mutate(Significance = ifelse(pval < 0.05, "Significant", "Not Significant"))
  
  # Compute positions for text (above bars)
  plot_df <- plot_df %>%
    mutate(y_position = MeanDiff + 0.01 * sign(MeanDiff))  # Adjust for visibility
  
  # Create the bar plot with p-value annotations
  p <- ggplot(plot_df, aes(x = reorder(Transcript, MeanDiff), y = MeanDiff, fill = Significance)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Reference line at 0
    geom_text(aes(y = y_position, label = sprintf("p=%.3f", pval)), size = 3, vjust = -0.5) + # Add p-values
    scale_fill_manual(values = c("Significant" = "red", "Not Significant" = "gray")) + # Color scheme
    labs(x = "Transcript", y = "Mean Difference", title = paste("Mean Differences", cd1,"vs", cd2)) +
    # some names are too long - issues to fix
    scale_x_discrete(labels = function(x) substr(x, 1, 10)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45))#, hjust = 1)) # Rotate x labels for readability
  ggplotly(p)
  
}

