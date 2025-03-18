# Core Libraries
library(shiny)                 # Web applications
library(ggplot2)               # Plotting
library(grid)                  # Grid-based graphics
library(dplyr)                 # Data manipulation
library(here)

# Genomics & Annotation
library(SummarizedExperiment)  # Genomic data structures
library(GenomicFeatures)       # Genomic features processing
library(GenomicRanges)         # Genomic ranges manipulation
library(org.Hs.eg.db)          # Human gene annotation
library(biomaRt)               # Querying Ensembl data

# Plotting & Visualization
library(plotgardener)          # Genomic plotting
library(wiggleplotr)           # Wiggle plot visualization
library(plyranges)             # Tidy-style genomic range manipulation
library(reshape2)
library(viridis)
library(RColorBrewer)

# Set the working directory to the project root
setwd(here::here())  
wd <- getwd() 

# Sources
source(file.path(wd,"app/R/plots.R"))
source(file.path(wd,"app/R/gene_descriptions.R"))


# GLOBAL 
# ---------------------------------------------------------------
# Load required data files
txp_counts_RSE <- readRDS(file.path(wd, "data/glinos_saturn_dtu.rds"))  # RNA-seq transcript count data
txp <- rowRanges(txp_counts_RSE) 
exons_per_txp_GRL <- readRDS(file.path(wd, "data/exons_per_txp_GRL.rds"))  # Exon  per transcript

se <- readRDS(file.path(wd, "/DE_results/glinos_saturn_dtu.rds"))



# gene symbols for selection in selectizeInput (ui) 
gene_ids <- sort(unique(rowData(se)$symbol))

# Load the transcript database (TxDb) for annotations
txdb1 <- loadDb(file.path(wd, "anno/gencode.v47.sqlite"))

# genome assembly metadata for plotgardener
gencode_assembly <- assembly(
  # Genome = "hg38",
  TxDb = txdb,
  # OrgDb = org.Hs.eg.db,
  gene.id.column = "GENEID",
  display.column = "GENEID",
  BSgenome = NULL
)

# Define full color palette for all conditions
colors <- c("#FFC20A",  # WT
            "#0C7BDC",  
            "#FB6B42",  
            "#3BD4D0",
            "#AF2AAF"   
)
unique_conditions <- unique(colData(se)$condition)
condition_colors <- setNames(colors[seq_along(unique_conditions)], unique_conditions)



### significant transcripts
# Identify significant transcripts with pval < 0.05
significant_indices <- rowData(se)$fitDTUResult_KD_vs_WT$pval < 0.05
significant_transcripts <- rownames(se)[significant_indices]

# Extract corresponding gene symbols and p-values
gene_symbols <- rowData(se)$symbol[significant_indices]
p_values <- rowData(se)$fitDTUResult_KD_vs_WT$pval[significant_indices]

# Create a data frame with Transcript, Gene Symbol, and P-value
gene_data <- data.frame(
  Transcript = significant_transcripts,
  Symbol = gene_symbols,
  P_Value = p_values
)
# Function: counts, DTU and pvalues 
# ---------------------------------------------------------------


calc_prop <- function(symbol ){
  cts <- assay(se, "counts")[ mcols(se)$symbol == symbol , ]
  prop <- t(cts) / colSums(cts)
  return(prop)
  
}

calc_mean_diff_DTU <-  function(symbol, cd1 = "WT", cd2 = "KD"){
  cts <- assay(se, "counts")[ mcols(se)$symbol == symbol , ]
  prop <- t(cts) / colSums(cts)
  # boxplot(prop[,2] ~ se$condition)
  
  n <- length(colnames(prop))
  mean_diffs_DTU <- sapply(1:n, function(j) mean(prop[se$condition == cd2, j]) - mean(prop[se$condition == cd1, j]))
  names(mean_diffs_DTU) <- rownames(cts)  # Assuming rownames(cts) contains transcript IDs
  return(mean_diffs_DTU)
}

get_pvals <- function(symbol){
  # hardcoded to KD vs WT  (change)
  pvals <- rowData(se[ mcols(se)$symbol == symbol , ])$fitDTUResult_KD_vs_WT$pval
  names(pvals) <-  rownames(rowData(se[ mcols(se)$symbol == symbol , ])$fitDTUResult_KD_vs_WT)
  return(pvals)
}

