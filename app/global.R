
# Core Libraries ---------------------------------------------------------------
library(shiny)                 # Web application framework
library(ggplot2)               # Plotting
library(grid)                  # Grid-based graphics
library(dplyr)                 # Data manipulation
library(here)                  # Project root management
library(scales)                # Scale functions


# Genomics & Annotation --------------------------------------------------------
library(SummarizedExperiment)  # Genomic data structures
library(GenomicFeatures)       # Process genomic features
library(GenomicRanges)         # Manipulate genomic ranges
library(org.Hs.eg.db)          # Human gene annotation
library(biomaRt)               # Query Ensembl data


# Plotting & Visualization -----------------------------------------------------
library(plotgardener)          # Genomic plotting
library(wiggleplotr)           # Wiggle plot visualization
library(plyranges)             # Tidy-style genomic range manipulation
library(reshape2)              # Data reshaping
library(viridis)               # Color palettes
library(RColorBrewer)          # Color palettes


# Set Working Directory and Source Files -------------------------------------
setwd(here::here())  
wd <- getwd() 

# Source additional scripts
source(file.path(wd, "app/R/plots.R"))
source(file.path(wd, "app/R/gene_descriptions.R"))


# GLOBAL Variables, Palettes & Helper Functions --------------------------------
# ------------------------------------------------------------------------------
# Palette for negative values: from red (-1) to gray (0)
pal_neg <- col_numeric(
  palette = c("#FF2900", "#ffd4cc", "#E0E0E0"),
  domain = c(-1, 0)
)

# Palette for positive values: from gray (0) to blue (1)
pal_pos <- col_numeric(
  palette = c("#E0E0E0", "#cae3ff", "#4FA4FF"),
  domain = c(0, 1)
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
se <- readRDS(file.path(wd, "data/glinos_saturn_dtu.rds"))

# Clean isoform_id: trim extra text and assign as row names
rowData(se)$isoform_id <- rowData(se)$isoform_id %>%
  str_replace("(_ENSG.*?)_ENSG.*", "\\1") %>%
  str_replace("(ENST.*?)_ENSG.*", "\\1")
rownames(se) <- rowData(se)$isoform_id

# get Significant DTU Results
sig_res <- rowData(se)[["fitDTUResult_exp_vs_ctrl"]] |>
  tibble::as_tibble() |>
  bind_cols(as.data.frame(rowData(se)[, 1:4])) |>
  filter(empirical_FDR < 0.1) |>
  select(gene_id, isoform_id, symbol, estimates, empirical_pval, empirical_FDR) |>
  arrange(empirical_pval)

# Clean isoform_id in significant results
sig_res <- sig_res |>
  mutate(
    isoform_id = str_replace(isoform_id, "(_ENSG.*?)_ENSG.*", "\\1"),
    isoform_id = str_replace(isoform_id, "(ENST.*?)_ENSG.*", "\\1")
  )

# Gene Symbols for UI
gene_ids <- sort(unique(sig_res$symbol))
symbol <- gene_ids[1]


# Load Transcript Database (TxDb) & Prepare for Annotation ---------------------
txdb <- loadDb(here("data", "flair_filter_transcripts.sqlite"))

# Group exons by transcript
exons <- exonsBy(txdb, by = "tx")

# Get transcript annotations and convert to tibble
txps <- AnnotationDbi::select(txdb, keys(txdb, "TXID"), c("TXNAME", "GENEID"), "TXID") |>
  tibble::as_tibble() |>
  mutate(TXID = as.character(TXID))

# Check consistency and rename exons by transcript name
length(exons)
all.equal(names(exons), txps$TXID)
names(exons) <- txps$TXNAME

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


# Create a data frame to display Transcript, Gene Symbol, and P-value in the UI.
significant_transcripts <- sig_res$isoform_id
gene_data <- data.frame(
  Transcript = significant_transcripts,
  Symbol = sig_res$symbol,
  P_Value = sig_res$empirical_pval,
  stringsAsFactors = FALSE
)


# Functions: Counts, DTU and p-values
# ------------------------------------------------------------------------------

# Calculate the proportion of counts per transcript for a given gene symbol
calc_prop <- function(symbol) {
  cts <- assay(se, "counts")[mcols(se)$symbol == symbol, ]
  prop <- t(cts) / colSums(cts)
  return(prop)
}

# Calculate mean differences for Differential Transcript Usage (DTU)
calc_mean_diff_DTU <- function(symbol, cd1 = "ctrl", cd2 = "exp") {
  cts <- assay(se, "counts")[mcols(se)$symbol == symbol, ]
  prop <- t(cts) / colSums(cts)
  
  # Compute difference in mean proportions between two conditions for each transcript
  n <- length(colnames(prop))
  mean_diffs_DTU <- sapply(1:n, function(j) {
    mean(prop[se$condition == cd2, j]) - mean(prop[se$condition == cd1, j])
  })
  names(mean_diffs_DTU) <- rownames(cts)  # Use transcript IDs as names
  return(mean_diffs_DTU)
}

# Extract p-values from the Saturn DTU analysis for a given gene symbol
get_pvals <- function(symbol, cd1 = "ctrl", cd2 = "exp") {
  # Identify the column in rowData that contains DTU results
  saturn_col <- paste0("fitDTUResult_", cd2, "_vs_", cd1)
  
  # Subset rowData for matching symbol and extract p-values
  pvals <- rowData(se[rowData(se)$symbol == symbol, ])[[saturn_col]]$pval
  names(pvals) <- rownames(rowData(se[rowData(se)$symbol == symbol, ]))
  return(pvals)
}
