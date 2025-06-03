# Core Libraries ----------------------------------------------------------------
library(shiny)                 # Web application framework
library(ggplot2)               # Plotting
library(grid)                  # Grid-based graphics
library(dplyr)                 # Data manipulation
library(here)                  # Project root management
library(scales)                # Scale functions
library(stringr)               # String manipulation
library(DT)                    # Interactive tables
library(plotly)                # Interactive plotting
library(bslib)                 # Bootstrap themes
library(purrr)                 # Functional programming

# Genomics & Annotation --------------------------------------------------------
library(SummarizedExperiment)  # Genomic data structures
library(GenomicFeatures)       # Process genomic features
library(GenomicRanges)         # Manipulate genomic ranges
library(org.Hs.eg.db)          # Human gene annotation
library(GO.db)                 # Gene Ontology terms
library(biomaRt)               # Query Ensembl data
library(Biostrings)            # DNA/RNA/protein sequences
library(BSgenome)              # Genome data handling
library(BSgenome.Hsapiens.UCSC.hg38)  # Human genome (hg38)
# ─── Packages ────────────────────────────────────────────────────────────────
library(clusterProfiler)
library(enrichplot)
library(patchwork)   # for side-by-side layouts
# Plotting & Visualization -----------------------------------------------------
library(plotgardener)          # Genomic plotting
library(wiggleplotr)           # Wiggle plot visualization
library(plyranges)             # Tidy-style genomic range manipulation
library(reshape2)              # Data reshaping
library(viridis)               # Color palettes
library(RColorBrewer)          # Color palettes
library(reactable)
library(shinycssloaders)

# Load Local Development Package -----------------------------------------------
setwd("/work/users/b/e/beacm/wiggleplotr")
devtools::load_all()

# Set Working Directory and Source Scripts -------------------------------------
setwd(here::here())  
wd <- getwd() 

# Source additional scripts
source(file.path(wd, "app/R/plots.R"))
source(file.path(wd, "app/R/gene_descriptions.R"))
source(file.path(wd, "app/R/exonModules.R"))
source(file.path(wd, "app/R/isoformModules.R"))
source(file.path(wd, "app/R/summaryModules.R"))
source(file.path(wd, "app/R/aux_postDTU.R"))
source(file.path(wd, "app/R/isoform_plots.R"))
source(file.path(wd, "app/R/exon_detection.R"))

# GLOBAL Variables, Palettes & Helper Functions --------------------------------
# ------------------------------------------------------------------------------

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


# Load SummarizedExperiment Object ---------------------------------------------
se <- readRDS(here::here("data", "glinos_saturn_dtu.rds"))

# # Load table Object with satuRn DTU Analysis 
# sig_res <- read.csv(here::here("data", "glinos_saturn_dtu.csv"), 
#                     stringsAsFactors = FALSE)|>
#             as_tibble()

# Load Transcript Database (TxDb) & Prepare for Annotation ---------------------
txdb <- loadDb(here("data","flair_filter_transcripts.sqlite"))

exons <- readRDS(here::here("data", "glinos_exons.rds"))
get_sig_res <- function(fdr_threshold){
  sig_res <- rowData(se)[["fitDTUResult_exp_vs_ctrl"]] |>
    tibble::as_tibble() |>
    dplyr::bind_cols(as.data.frame(rowData(se)[,1:4])) |>
    dplyr::filter(empirical_FDR < fdr_threshold) |>
    dplyr::select(gene_id, isoform_id, symbol, estimates, empirical_pval, empirical_FDR) |>
    dplyr::arrange(empirical_pval)

  sig_res <-  sig_res %>%
    dplyr::mutate(sign = sign(estimates))
  return(sig_res)
  
}

sig_res <- get_sig_res(0.05) #default 0.05 fdr

# # Gene Symbols for UI
# gene_ids <- sort(unique(sig_res$symbol))
# symbol <- gene_ids[1]


#  Assembly for plotgardener
db_assembly <- plotgardener::assembly(
  Genome = "hg38",
  TxDb = txdb,
  OrgDb = org.Hs.eg.db,
  gene.id.column = "GENEID",
  display.column = "GENEID",
  BSgenome = NULL
)


# Condition Colors 
colors <- c(
  "#3BD4D0",  # WT
  "#AF2AAF"   # KD PTBP1
)
unique_conditions <- unique(colData(se)$condition)
condition_colors <- setNames(colors[seq_along(unique_conditions)], unique_conditions)




get_x_flat <- function(sig_res){
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
calc_prop <- function(symbol, sig_res) {
  cts <- assay(se, "counts")[mcols(se)$symbol == symbol, ]
  prop <- t(cts) / colSums(cts)
  
  sig_ids <- sig_res %>%
    filter(symbol == symbol) %>%
    pull(isoform_id)
  
  prop <- prop[, colnames(prop) %in% sig_ids, drop = FALSE]
  return(prop)
}

# Calculate mean differences for Differential Transcript Usage (DTU)
calc_mean_diff_DTU <- function(gene_symbol, sig_res,
                               cd1 = "ctrl", cd2 = "exp") {
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
get_pvals <- function(gene_symbol, sig_res, 
                      cd1 = "ctrl", cd2 = "exp") {
  # Identify the column in rowData that contains DTU results
  saturn_col <- paste0("fitDTUResult_", cd2, "_vs_", cd1)
  
  # Subset rowData for matching symbol and extract p-values
  pvals <- rowData(se[rowData(se)$symbol == gene_symbol, ])[[saturn_col]]$empirical_pval
  names(pvals) <- rownames(rowData(se[rowData(se)$symbol == gene_symbol, ]))
  sig_ids <- sig_res %>%
    filter(symbol == gene_symbol) %>%
    pull(isoform_id)
  
  pvals <- pvals[ sig_ids ]
  return(pvals)
}


