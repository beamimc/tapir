build_ui <- function(conditions) {
  page_sidebar(
    title   = "Isoform Analysis",      
    theme   = bs_theme(version = 5, 
                       bootswatch = "cosmo"),  # Bootstrap 5
  
    sidebar = sidebarPanel(
      width = 12,   
      fluidRow(
        column(
          width = 6,
          radioButtons(
            inputId  = "ref",
            label    = "Reference",
            choices  = conditions$cd1,
            selected = conditions$cd1[1]
          )
        ),
        column(
          width = 6,
          radioButtons(
            inputId  = "condition2",
            label    = "Condition 2",
            choices  = conditions$cd2,
            selected = conditions$cd2[1]
          )
      )
      ),
      br(), br(),br(),
      sliderInput(
        "fdr_threshold", "FDR threshold:",
        min = 0, max = 0.5, value = 0.05, step = 0.01
      ),
      actionButton("apply_fdr", "Apply FDR Filter"),
      br(), br(),br(),
      radioButtons(
        inputId = "exon_filter",
        label = "Select structural change:",
        choices = c("Downregulated exon", "Upregulated exon", "UTR change"),
        selected = "Downregulated exon"
      )
      
      
      # selectizeInput(
      #   "global_genes", "Select genes:",
      #   choices  = gene_ids,
      #   multiple = TRUE,
      #   options  = list(create = TRUE)
      # )
      # actionButton("select_all", "Select all"),
      # tags$hr(),
      # textInput("view_name", "Save selection as:"),
      # actionButton("save_view", "Save view"),
      # selectInput("load_view", "Load view:", choices = NULL),
      # actionButton("apply_view", "Apply"),
      # actionButton("delete_view", "Delete")
    ),
    
    
    tabsetPanel(
      id = "main_tabs",
      tabPanel("DTU Exploration",   value = "dtu",     isoformAnalysisUI("isoform")),
      tabPanel("Exon-Level",        value = "exon",    exonLevelUI("exon")),
      tabPanel("Summary",     value = "summary", summaryStatsUI("summary"))
    )
  )
}