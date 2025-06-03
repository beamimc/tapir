####################################
# Post DTU Analysis Auxiliary Functions
####################################

# Load required libraries
library(Biostrings)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(plotly)

#######################################################
# get_bp_count: Get Nucleotide Counts
#######################################################
#' Get Base Pair Counts
#'
#' Computes the count (frequency) of each nucleotide (A, U, C, G) for each sequence 
#' in an RNAStringSet object.
#'
#' @param stringSet An RNAStringSet object.
#' @return A matrix with counts of each nucleotide.
get_bp_count <- function(stringSet){ # input is RNAStringSet object
  bp_counts <- oligonucleotideFrequency(stringSet, width = 1)
}

#######################################################
# get_bp_percent: Get Nucleotide Percentages
#######################################################
#' Get Base Pair Percentages
#'
#' Computes the percentage of each nucleotide (A, U, C, G) in each sequence of an
#' RNAStringSet object.
#'
#' @param stringSet An RNAStringSet object.
#' @return A matrix with percentage (proportions) of each nucleotide.
get_bp_percent <- function(stringSet){
  bp_counts <- oligonucleotideFrequency(stringSet, width = 1)
  bp_percent <-  bp_counts / rowSums(bp_counts)
}

#######################################################
# get_dibp_features: Get Di-nucleotide Features
#######################################################
#' Get Di-nucleotide Features
#'
#' Calculates the frequency of di-nucleotides (pairs, e.g. AA, AU, AC, etc.) as a 
#' percentage for each sequence in an RNAStringSet.
#'
#' @param stringSet An RNAStringSet object.
#' @return A matrix with percentages for each di-nucleotide.
get_dibp_features <- function(stringSet){
  di_freq <- oligonucleotideFrequency(stringSet, width = 2)
  di_perc <- di_freq / (width(stringSet) - 1) # equivalent to rowSums per sequence normalization
}

#######################################################
# get_tribp_features: Get Tri-nucleotide Features
#######################################################
#' Get Tri-nucleotide Features
#'
#' Calculates the frequency of tri-nucleotides (triplets, e.g. AAA, AAU, AAC, etc.) as a 
#' percentage for each sequence in an RNAStringSet.
#'
#' @param stringSet An RNAStringSet object.
#' @return A matrix with percentages for each tri-nucleotide.
get_tribp_features <- function(stringSet){
  tri_freq <- oligonucleotideFrequency(stringSet, width = 3)
  tri_perc <- tri_freq / (width(stringSet) - 2)
}

#######################################################
# get_windows: Create Non-Overlapping Windows and Calculate % bp
#######################################################
#' Get Windows of Upstream Regions
#'
#' Divides the upstream sequences (provided as a GRanges object) into non-overlapping windows 
#' and calculates the percentage of each nucleotide per window.
#'
#' @param upstr_ranges A GRanges object containing upstream regions.
#' @param window_width Integer width of each window (default is 10).
#' @param width_upstream Total width of upstream region (default is 100); must match flankupsteam.
#' @return A matrix with nucleotide percentages for each window (features named by window and nucleotide).
get_windows<- function(upstr_ranges, # input is GRanges object 
                       window_width = 10,
                       width_upstream = 100 ## make sure it matches flankupsteam
){
  
  n_windows <- width_upstream / window_width # number of windows generated given the total region length
  windows <- upstr_ranges |> tile_ranges(window_width) # divide into windows
  seq_windows <-  Hsapiens |>  
    getSeq(windows) |> 
    RNAStringSet() 
  # RNAStringSet result has n_windows * (# of exons) rows, e.g. for 158 exons and 10 windows, we get 1580 rows.
  windows_bp <- get_bp_percent(seq_windows) # calculate %bp for all windows
  
  result_list <- vector("list", n_windows)
  for (window_index in 1:n_windows){ 
    indices <- seq(from = window_index, 
                   by = n_windows, 
                   length.out = length(upstr_ranges))
    
    sub_matrix <- windows_bp[indices, ]
    colnames(sub_matrix) <- paste0("w", window_index, "_", colnames(sub_matrix))
    
    result_list[[window_index]] <- sub_matrix
  }
  final_matrix <- do.call(cbind, result_list)
  return(final_matrix)
}

#######################################################
# get_sliding_windows: Create Overlapping Windows and Calculate % bp
#######################################################
#' Get Sliding Windows of Upstream Regions
#'
#' Divides upstream sequences (as a GRanges object) into overlapping windows and 
#' calculates the percentage of each nucleotide per window.
#'
#' @param upstr_ranges A GRanges object containing upstream regions.
#' @param window_width Integer width of each window (default is 10).
#' @param width_upstream Total width of upstream region (default is 100); must match flankupsteam.
#' @param overlap Number of overlapping nucleotides between adjacent windows (default is 5).
#' @return A matrix of features created from overlapping windows.
get_sliding_windows <- function(upstr_ranges, # input is GRanges object 
                                window_width = 10,
                                width_upstream = 100, ## should match flankupsteam - could be improved to avoid hardcoding 
                                overlap = 5 # use overlap = 0 for non-overlapping windows
){
  
  # Calculate step size: if step == window_width, slide_ranges is equivalent to tile_ranges 
  step <-  window_width - overlap 
  n_windows <- ceiling((width_upstream - window_width) / step) + 1
  
  windows <- upstr_ranges |> slide_ranges( width = window_width,  # GRanges
                                           step = step) # divide into windows
  seq_windows <-  Hsapiens |>
    getSeq(windows) |>
    RNAStringSet()  # RNAStringSet
  # The RNAStringSet result has n_windows * (# of exons) rows.
  windows_bp <- get_bp_percent(seq_windows) # calculate %bp for all windows

  # windows_meta <- as.data.frame(windows)[, c("strand", "partition")]
  # 
  # # Combine the matrix and the metadata by column-binding
  # combined <- cbind(windows_bp, windows_meta)
  # 
  # windows_df <- as.data.frame(combined)
  tb <- as_tibble(windows)            # seqnames, start, end, strand, partition, etc.
  windows_meta <- tb[, c("strand","partition")]
  combined     <- cbind(windows_bp, windows_meta)
  windows_df   <- as.data.frame(combined)
  # Group by partition (and strand, if desired) and create a new column "window_label"
  windows_df <- windows_df %>%
    group_by(partition, strand) %>% 
    mutate(window_order = if (dplyr::first(as.character(strand)) == "+") {
      row_number()                # For '+' strand: 1, 2, 3, ...
    } else {
      dplyr::n() - row_number() + 1       # For '-' strand: reverse order
    },
    window_label = paste0("w", window_order)) %>%
    ungroup()
  
  # Create a new data frame with one row per partition and columns for each window-nucleotide combination
  final_df <- windows_df %>%
    # Pivot the nucleotide columns (A, C, G, U) into long format
    pivot_longer(
      cols = c("A", "C", "G", "U"),
      names_to = "Nucleotide",
      values_to = "Value"
    ) %>%
    # Create a new column combining window label and nucleotide (e.g., "w1_A")
    mutate(Window_Nucleotide = paste0(window_label, "_", Nucleotide)) %>%
    # Select only the columns needed for pivoting: partition, the combined label, and value
    select(partition, Window_Nucleotide, Value) %>%
    # Pivot wider so that each Window_Nucleotide becomes its own column
    pivot_wider(
      id_cols = partition,
      names_from = Window_Nucleotide,
      values_from = Value
    ) %>%
    # Remove the partition column
    select(-partition)
  
  final_df
  matrix_df <- as.matrix(final_df)
  return(matrix_df)
}

#######################################################
# create_exon_df: Create DataFrame with Exon Features
#######################################################
#' Create Exon Feature DataFrame
#'
#' Combines nucleotide percentage features, di-nucleotide and tri-nucleotide features, 
#' and sliding window features for exons. It also adds an auxiliary ID and an optional label.
#'
#' @param seq_exons An RNAStringSet object of exon sequences.
#' @param upstr_exons A GRanges object representing upstream exon regions.
#' @param label (Optional) A binary label for machine learning purposes.
#' @return A data frame with all the computed features.
create_exon_df <- function(seq_exons, 
                           upstr_exons,
                           label=NA # binary label for ML
) {
  # Convert the RNAStringSet to a data frame and add the auxiliary ID
  df <- as.data.frame(seq_exons) |>
    dplyr::mutate(aux_id = metadata(seq_exons)$aux_id) |>
    dplyr::mutate(label = label) |>
    dplyr::rename(seq = x)
  
  # Add % base pair features, di-nucleotide, tri-nucleotide, and sliding window features
  df <- cbind(df, get_bp_percent(seq_exons))
  df <- cbind(df, get_dibp_features(seq_exons))
  df <- cbind(df, get_tribp_features(seq_exons))
  df <- cbind(df, get_sliding_windows(upstr_exons,  # sliding windows
                                      window_width = 10,
                                      width_upstream = width_upstream, ## should match flankupsteam 
                                      overlap = 5
  )) 
  
  return(df)
}

#######################################################
# barplot_bppercent: Plot Base Pair Percentage Barplot
#######################################################
#' Bar Plot for Base Pair Percentage
#'
#' Creates a bar plot of the mean percentage of base pairs for each bp combination, 
#' with individual data points overlaid.
#'
#' @param matrix A matrix of bp percentage features (rows represent samples).
#' @return A ggplot2 plot of the bp percentages.
barplot_bppercent <- function(matrix) {
  # Convert the matrix to a data frame and pivot it to long format
  df <- as.data.frame(matrix)
  df$Sample <- rownames(df)
  
  df_long <- pivot_longer(df, cols = -Sample, names_to = "bp", values_to = "Value")
  
  # Calculate the mean value for each bp
  mean_df <- df_long %>% 
    group_by(bp) %>% 
    summarize(Mean = mean(Value)) %>%
    arrange(desc(Mean))
  df_long$bp <- factor(df_long$bp, levels = mean_df$bp)
  
  # Create the ggplot2 plot
  p <- ggplot() +
    # Bar layer for the mean of each bp
    geom_bar(data = mean_df, aes(x = bp, y = Mean), 
             stat = "identity", fill = "darkgray", width = 0.7) +
    # Jitter layer for the individual data points
    geom_jitter(data = df_long, aes(x = bp, y = Value), 
                width = 0.15, size = 2, color = "#2972b6") +
    # Minimal theme and rotate x-axis labels 45 degrees
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "bp combination", 
         y = "%", 
         title = "Percentage of bp in exon's upstream region (150bp)")
  
  # Print the plot
  print(p)
}

#######################################################
# barplot_bppercent2: Grouped Barplot for bp Percentage Comparison
#######################################################
#' Bar Plot for Base Pair Percentage Comparison
#'
#' Creates a grouped bar plot comparing bp percentage features between two groups 
#' (e.g., spliced vs. non-spliced), with jittered individual data points.
#'
#' @param mat1 A matrix of bp percentage features for group1.
#' @param mat2 A matrix of bp percentage features for group2.
#' @param group1_label Label for the first group (default "spliced").
#' @param group2_label Label for the second group (default "non-spliced").
#' @param group1_color Color for group1 (default "#F84040").
#' @param group2_color Color for group2 (default "skyblue").
#' @return A plotly object with the grouped bar plot.
barplot_bppercent2 <- function(mat1, mat2, 
                               group1_label = "spliced", 
                               group2_label = "non-spliced",
                               group1_color = "#F84040", 
                               group2_color = "skyblue") {
  
  # Convert each matrix to a data frame, add sample names and a group identifier
  df1 <- as.data.frame(mat1)
  df1$Sample <- rownames(df1)
  df1$Group <- group1_label
  
  df2 <- as.data.frame(mat2)
  df2$Sample <- rownames(df2)
  df2$Group <- group2_label
  
  # Combine the two data frames
  df_all <- bind_rows(df1, df2)
  
  # Pivot the combined data frame to long format
  df_long <- pivot_longer(df_all, 
                          cols = -c(Sample, Group), 
                          names_to = "bp", 
                          values_to = "Value")
  
  # Calculate the mean for each bp within each group
  mean_df <- df_long %>% 
    group_by(bp, Group) %>% 
    summarize(Mean = mean(Value), .groups = "drop")
  
  # Calculate group1 means and sort bps by them (largest to smallest)
  group1_means <- mean_df %>% 
    filter(Group == group1_label) %>% 
    arrange(desc(Mean))
  
  # Update factor levels based on group1 ordering
  df_long$bp <- factor(df_long$bp, levels = group1_means$bp)
  mean_df$bp <- factor(mean_df$bp, levels = group1_means$bp)
  
  
  # Create the ggplot: grouped bars with jittered data points
  p <- ggplot() +
    # Bar layer: use position_dodge to separate groups
    geom_bar(data = mean_df, 
             aes(x = bp, y = Mean, fill = Group), 
             stat = "identity",
             position = position_dodge(width = 0.8), 
             width = 0.7) +
    # Jitter layer: position_jitterdodge to add jitter while respecting dodge groups
    geom_jitter(data = df_long, 
                aes(x = bp, y = Value, color = Group), 
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
                size = 1) +
    scale_fill_manual(values = setNames(c(group1_color, group2_color), 
                                        c(group1_label, group2_label)))+
    scale_color_manual(values = setNames(c(group1_color, group2_color), 
                                         c('gray', 'gray')))+
    
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "bp combination", 
         y = "%", 
         title = "Percentage of bp in exon's upstream region")
  
  ggplotly(p)
}

#######################################################
# plot_exon_summary: Plot Exon Summary for a Selected Nucleotide
#######################################################
#' Plot Exon Summary
#'
#' Creates a summary plot for exon sliding window bp percentages for a chosen nucleotide.
#' The function processes both spliced and non-spliced exon window matrices, pivots the data to long 
#' format and plots the mean percentage per window along with standard error bars.
#'
#' @param spliced_exons_windows A matrix of sliding window bp percentages for spliced exons.
#' @param nonspliced_exons_windows (Optional) A matrix of sliding window bp percentages for non-spliced exons.
#' @param nucleotide The nucleotide to plot (default is "U").
#' @param custom_colors Named vector of colors for the groups (default: Spliced = "#F84040", Non-Spliced = "skyblue").
#' @param size A numeric value for scaling window labels (default 100).
#' @param step The step size used when creating the windows (default 5).
#' @return A ggplot2 plot showing bp percentage trends across upstream windows.
plot_exon_summary <- function(spliced_exons_windows, 
                              nonspliced_exons_windows = NULL, 
                              nucleotide = "U",
                              custom_colors = c("Downregulated" = "#F84040", 
                                                "Non-Regulated" = "skyblue"),
                              size = 100, 
                              step = 5)
{
  # -------------------------------
  # Process the spliced_exons_windows data
  # -------------------------------
  selected_cols_spliced <- grep(nucleotide, colnames(spliced_exons_windows), value = TRUE)
  if (length(selected_cols_spliced) == 0) {
    stop("No columns found matching the pattern ", nucleotide, " in spliced_exons_windows")
  }
  
  sub_matrix_spliced <- spliced_exons_windows[, selected_cols_spliced, drop = FALSE]
  
  # Convert spliced matrix to long format and label dataset as "Spliced"
  df_spliced <- melt(sub_matrix_spliced, varnames = c("Row", "Window"), value.name = "Value")
  df_spliced$Dataset <- "Downregulated"
  
  # -------------------------------
  # Process the nonspliced_exons_windows data (if provided)
  # -------------------------------
  if (!is.null(nonspliced_exons_windows)) {
    selected_cols_nonspliced <- grep(nucleotide, colnames(nonspliced_exons_windows), value = TRUE)
    if (length(selected_cols_nonspliced) == 0) {
      stop("No columns found matching the pattern ", nucleotide, " in nonspliced_exons_windows")
    }
    
    sub_matrix_nonspliced <- nonspliced_exons_windows[, selected_cols_nonspliced, drop = FALSE]
    
    # Convert nonspliced matrix to long format and label dataset as "Non-Spliced"
    df_nonspliced <- melt(sub_matrix_nonspliced, varnames = c("Row", "Window"), value.name = "Value")
    df_nonspliced$Dataset <- "Non-Regulated"
    
    # Combine both datasets
    long_df <- rbind(df_spliced, df_nonspliced)
    
    # Ensure the Window factor levels come from the spliced dataset (assumed common order)
    long_df$Window <- factor(long_df$Window, levels = selected_cols_spliced)
    
  } else {
    long_df <- df_spliced
    long_df$Window <- factor(long_df$Window, levels = selected_cols_spliced)
  }
  
  ## Create x-axis labels based on window number extraction
  levels_x <- levels(long_df$Window)
  substr <- paste0("_", nucleotide)
  # Remove the "w" and convert the remaining part to numeric
  window_number <- as.numeric(sub(substr,"",sub("w", "", levels_x)))
  
  # Compute the new labels
  x_labels <- -size + (window_number - 1) * step
  # -------------------------------
  # Create the plot
  # -------------------------------
  # - stat_summary (with fun = mean) plots a line of mean values per group (Dataset)
  # - stat_summary (with fun.data = mean_se) adds error bars based on the standard error
  p <- ggplot(long_df, aes(x = Window, y = Value, color = Dataset)) +
    stat_summary(fun = mean, geom = "line", aes(group = Dataset), size = 1) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    scale_color_manual(values = custom_colors) +
    # Modify the x-axis labels:
    scale_x_discrete(labels = x_labels) +
    labs(title = paste(nucleotide, "% per window"),
         x = "Window(upstream position)", y = "%bp") +
    theme_classic()
  
  print(p)
  
  # Optionally, return the long data frame for inspection
  invisible(long_df)
  return(p)
}

plot_window_comparison <- function(upstr_downreg_exons,
                                   upstr_nonreg_exons,
                                   width_upstream=100){
  nucs <- c("G","A","U","C")
  #windoes is #exons x #windows*4 (not summarized by window yet) # 
  # downreg_exons_windows <- get_sliding_windows(upstr_downreg_exons,#slid windows # 
  #                                              window_width = 10, # 
  #                                              width_upstream = width_upstream, ## maches flankupsteam # 
  #                                              overlap = 5) 
  # 
  # nonreg_exons_windows <- get_sliding_windows(upstr_nonreg_exons,#slid windows 
  #                                             window_width = 10, 
  #                                             width_upstream = width_upstream, ## maches flankupsteam 
  #                                             overlap = 5 )
  downreg_exons_windows <- upstr_downreg_exons
  nonreg_exons_windows <- upstr_nonreg_exons
  
  plots <- lapply(nucs, function(nt) {
    plot_exon_summary(downreg_exons_windows, 
                      nonreg_exons_windows, 
                      nucleotide = nt) +
      ggtitle(nt)
  })
  
  # 3) Combine and enforce common y‐axis --------------------------------------
  combined <- wrap_plots(plots, ncol = 4) & 
    coord_cartesian(ylim  = c(0,0.5)) & 
    theme(legend.position = "bottom")
  
  combined
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
