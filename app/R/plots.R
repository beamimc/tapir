

# Function: Plot isoforms per gene using wiggleplotr
# ---------------------------------------------------------------
exon_level_plot <- function(gene_symbol) {

  # sig_res |> filter(symbol=="TPM2")|> pull(gene_id)|>unique()
  # #tx from the gene present in the se
  # present_tx <-  names(se)[ rowData(se)$gene_id ==gene ]
  # # get the sig tx 
  sig_tx <- sig_res |> filter(symbol== gene_symbol) |> select(isoform_id)
  p <- wiggleplotr::plotTranscripts(exons[sig_tx$isoform_id], rescale_introns = TRUE)
  return(p)
}


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
  
  hilite <- df_values %>% dplyr::select(transcript, color)
  
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
  prop_subset <- prop[subset_indices, ,  drop = FALSE]
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

get_transcript_colors <- function(transcripts, palette_name = "Set2") {
  transcripts <- sort(unique(transcripts))  # ensure consistent order
  n_colors <- max(3, length(transcripts))
  
  # If Set2 runs out, fall back to hue palette
  if (n_colors <= 8) {
    palette <- RColorBrewer::brewer.pal(n = n_colors, name = palette_name)
  } else {
    palette <- scales::hue_pal()(n_colors)
  }
  
  setNames(palette[seq_along(transcripts)], transcripts)
}


# Function: spaguetti plot for transcript comparion KD vs WT (change hardcode)
# ---------------------------------------------------------------

line_plot_txp_comparison  <- function(prop, txp_colors, cd1 = "ctrl", cd2 = "exp"){
  
  plot_data <- get_plot_data(prop)
  
  summary_data <- plot_data |>
    group_by(Transcript, condition) |>
    summarise(
      mean_proportion = mean(Proportion),
      sd_proportion = sd(Proportion),
      n = dplyr::n(),
      se_proportion = sd_proportion / sqrt(n)  # Standard error
    ) 
  
  
  p <- ggplot(summary_data, aes(x = condition, y = mean_proportion,
                                group = Transcript, color = Transcript)) +
    geom_line() +
    geom_point() + 
    geom_errorbar(aes(ymin = mean_proportion - se_proportion, 
                      ymax = mean_proportion + se_proportion), 
                  width = 0.2) +           # Error bars for each point
    scale_color_manual(values = txp_colors) +
    theme_classic() +
    labs(x = "Condition", y = "Proportion"
         # title = "Proportion of Transcripts by Condition"
         ) +
    
    theme(axis.text.x = element_text(angle = 0))
  ggplotly(p)%>%
    layout(legend = list(
      orientation = "h",      # horizontal keys
      x           = 0.5,      # centered horizontally
      xanchor     = "center",
      y           = 1.05,     # just above the plotting area
      yanchor     = "bottom"
    ))
  # %>% 
  #   layout(margin = list(t = 100))
}


# Function: barplot of mean difs KD vs WT (change hardcode)
# ---------------------------------------------------------------
barplot_meandifs <- function(mean_diffs_DTU, pvals,txp_colors, cd1 = "ctrl", cd2 = "exp") {
  
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

  
  # Compute positions for text (above bars)
  plot_df <- plot_df %>%
    mutate(
      Significance = ifelse(pval < 0.05, "Significant", "Not Significant"),
      y_position = MeanDiff + 0.01 * sign(MeanDiff)
    )
  
  # Plot
  p <- ggplot(plot_df, aes(x = reorder(Transcript, MeanDiff), y = MeanDiff, fill = Transcript)) +
    geom_bar(stat = "identity", width = 0.7) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(y = y_position*1.2, label = sprintf("p=%.0e", pval)), size = 3, vjust = -1) +
    scale_fill_manual(values = txp_colors) +  # ✅ fill not color
    labs(x = "Transcript", y = "Mean Difference") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplotly(p) %>%
    layout(showlegend = FALSE)
    # layout(legend = list(
    #   orientation = "h",
    #   x = 0.5,
    #   xanchor = "center",
    #   y = 1.05,
    #   yanchor = "bottom"
    # ))
}



plot_updownstream_windows <- function(df1, df2,
                                size      = 100,
                                step      = 5,
                                gap       = 8,
                                rect_fill = "grey70",
                                line_size = 1,
                                palette   = c(A="#F84040", 
                                              C="skyblue", 
                                              G="#FFB400", 
                                              U="#06D6A0"),
                                exon_label = "exon") {
  # reshape to long form, keeping a replicate ID
  longify <- function(df, set_label) {
    as.data.frame(df) %>%
      tibble::rownames_to_column("replicate") %>%
      pivot_longer(
        cols = -replicate,
        names_to  = c("window", "nt"),
        names_sep = "_",
        values_to = "value"
      ) %>%
      mutate(set = set_label)
  }
  
  l1   <- longify(df1, "upstream")
  l2   <- longify(df2, "downstream")
  nwin <- length(unique(l1$window))   # e.g. 8
  xmin <- nwin + 0.5
  xmax <- nwin + gap + 0.5
  ymin <- 0.2
  ymax <- 0.4
  
  # assign numeric x positions
  l1 <- l1 %>%
    mutate(x = as.numeric(sub("w","", window)))
  l2 <- l2 %>%
    mutate(x = as.numeric(sub("w","", window)) + nwin + gap)
  
  combined <- bind_rows(l1, l2)
  
  # compute breaks and labels
  window_number <- 1:nwin
  # x positions
  xpos1 <- window_number
  xpos2 <- window_number + nwin + gap
  # labels
  labs1 <- -size + (window_number ) * step
  labs2 <-          (window_number ) * step
  
  all_breaks <- c(xpos1, xpos2)
  all_labels <- c(labs1, labs2)
  
  # now plot
  ggplot(combined,
         aes(x        = x,
             y        = value,
             color    = nt,
             linetype = set,
             group    = interaction(set, nt))) +
    # gap rectangle
    annotate("rect",
             xmin = xmin, xmax = xmax,
             ymin = ymin, ymax = ymax,
             fill = rect_fill) +
    # add text centered in that rectangle
    annotate("text",
             x     = (xmin + xmax)/2,
             y     = (ymin + ymax)/2,
             label = exon_label,
             size  = 3,        # text size
             color = "black",  # choose contrasting color
             fontface = "bold"
    ) +
    # mean lines
    stat_summary(fun      = mean,
                 geom     = "line",
                 size     = line_size) +
    # error bars
    stat_summary(fun.data = mean_se,
                 geom     = "errorbar",
                 width    = 0.2) +
    scale_color_manual(values = palette) +
    scale_linetype_manual(values = c(upstream="solid",
                                     downstream="solid"),
                          guide = FALSE) +
    scale_x_continuous(breaks = all_breaks,
                       labels = all_labels) +
    scale_y_continuous(limits = c(0,0.7), expand = c(0,0)) +
    theme_classic() +
    labs(x = "Relative location (bp)",
         y = "Nucleotide percentage (%)",
         color    = "Nucleotide") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
}

go_plot2 <- function(gene_symbol_list){
  
  # ─── Helper to make an “empty” panel ────────────────────────────────────────
  blank_panel <- function(txt){
    ggplot() +
      annotate("text", x = .5, y = .5, label = txt, size = 6, hjust = .5) +
      theme_void()
  }
  
  # ─── Build one plot per ontology ─────────────────────────────────────────────
  plots <- lapply(c("BP", "MF", "CC"), function(ont) {
    ego <- enrichGO(
      gene       = gene_symbol_list,
      OrgDb      = org.Hs.eg.db,
      keyType    = "SYMBOL",
      ont        = ont,
      pvalueCutoff = 0.05,
      readable   = TRUE
    )
    
    # If there's nothing enriched, return a blank panel
    if (is.null(ego) || nrow(ego@result) == 0) {
      blank_panel(paste0(ont, " — no enriched terms"))
      
      # Otherwise make the cnetplot
    } else {
      cnetplot(
        ego,
        showCategory = 10,
        node_label   = "category"
      ) +
        ggtitle(ont) +
        theme(plot.title = element_text(hjust = 0.5))
    }
  })
  
  # ─── Combine side-by-side ────────────────────────────────────────────────────
  combined <- wrap_plots(plots, ncol = 3)
  combined ##ggplot
  
}

go_plot <- function(gene_symbol_list) {
  # Helper: blank panel when no enriched terms
  blank_panel <- function(txt) {
    ggplot() +
      annotate("text", x = .5, y = .5, label = txt, size = 6, hjust = .5) +
      theme_void()
  }
  
  # Run enrichGO + return barplot per ontology
  plots <- lapply(c("BP", "MF", "CC"), function(ont) {
    ego <- enrichGO(
      gene          = gene_symbol_list,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pvalueCutoff  = 0.05,
      readable      = TRUE
    )
    print(head(ego@result))
    
    
    if (is.null(ego) || nrow(ego@result) == 0) {
      return(blank_panel(paste0(ont, " — no enriched terms")))
    }
    
    return(
      barplot(
        ego,
        showCategory = 10,
        title        = ont,
        font.size    = 10
      )
    )
  })
  
  # Combine the 3 barplots side-by-side
  wrap_plots(plots, ncol = 3)
}

