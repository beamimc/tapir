library(shiny)
library(DT)
library(plotly)
library(bslib)

ui <- page_sidebar(
  title   = "Isoform Analysis",      # appears in the top navbar
  theme   = bs_theme(version = 5, bootswatch = "flatly"),   # Bootstrap 5
  
  # ─── Your off-canvas sidebar panel ─────────────────────────────────────────
  sidebar = sidebarPanel(
    width = 12,    
    selectizeInput(
      "global_genes", "Select genes:",
      choices  = gene_ids,
      multiple = TRUE,
      options  = list(create = TRUE)
    )
    # actionButton("select_all", "Select all"),
    # tags$hr(),
    # textInput("view_name", "Save selection as:"),
    # actionButton("save_view", "Save view"),
    # selectInput("load_view", "Load view:", choices = NULL),
    # actionButton("apply_view", "Apply"),
    # actionButton("delete_view", "Delete")
  ),
  
  # ─── Main body: tabset of your modules ──────────────────────────────────────
  tabsetPanel(
    id = "main_tabs",
    tabPanel("DTU Exploration",   value = "dtu",     isoformAnalysisUI("isoform")),
    tabPanel("Exon-Level",        value = "exon",    exonLevelUI("exon")),
    tabPanel("Summary Stats",     value = "summary", summaryStatsUI("summary"))
  )
)