library(shiny)
library(bslib)
library(shinycssloaders)
library(DT)
library(plotly)
library(reactable)
library(dplyr)
library(ggplot2)
library(scales)
library(purrr)
library(stringr)
library(grid)
library(patchwork)

library(SummarizedExperiment)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(GO.db)
library(biomaRt)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

library(clusterProfiler)
library(enrichplot)
library(plyranges)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(here)

# Load Local Development Package -----------------------------------------------
setwd("/work/users/b/e/beacm/wiggleplotr")
devtools::load_all()

# Set Working Directory and Source Scripts -------------------------------------
setwd(here::here())  
wd <- getwd() 

app <- function(se, exons, app_dir = ".") {
  
  #load functions in app/R/
  helpers_env <- new.env()
  r_scripts <- list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(r_scripts, function(f) source(f, local = helpers_env)))
  list2env(as.list(helpers_env), envir = environment())
  
  #load data
  data <- list(se = se, exons = exons)
  se <- data$se 
  
  conditions <- parse_saturnDTU_conditions(see)
  condition_choices_df <- purrr::map_dfr(seq_len(nrow(conditions)), function(i) {
    c1 <- conditions$cd1[i]
    c2 <- conditions$cd2[i]
    tibble::tibble(
      label = c(paste(c1, "vs", c2), paste(c2, "vs", c1)),
      cd1 = c(c1, c2),
      cd2 = c(c2, c1)
    )
  })
  data$condition_choices <-  condition_choices_df
  
  
  source(file.path(app_dir, "ui.R"), local = TRUE)
  ui <- build_ui(condition_choices)
  server_func <- source(file.path(app_dir, "server.R"), local = TRUE)$value
  
  #launch App
  shinyApp(
    ui = ui,
    server = function(input, output, session) {
      server_func(input, output, session, data)
    }
  )
}
