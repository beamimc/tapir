build_ui <- function(condition_choices) {
  
  page_sidebar(
    title   = "SPLain",      
    theme   = bs_theme(version = 5, bootswatch = "cosmo"),
    
    sidebar = sidebarPanel(
      width = 12,
      
      fluidRow(
        column(
          width = 6,
          radioButtons(
            inputId  = "cd1",
            label    = "Reference",
            choices  = c("Loading..." = "loading_placeholder",
                         selected = "ctrl")   # filled in server
          )
        ),
        column(
          width = 6,
          radioButtons(
            inputId  = "cd2",
            label    = "Contrast",
            choices  = c("Select reference first" = "loading_placeholder",
                         selected = "exp")  # updated dynamically based on cd1
          )
        )
      ),
      actionButton("apply_pair", "Apply comparison"),
      textOutput("current_comparison"),
      
      
      br(), br(),
      
      sliderInput("fdr_threshold", "FDR threshold:",
                  min = 0, max = 0.5, value = 0.05, step = 0.01),
      
      actionButton("apply_fdr", "Apply FDR Filter"),
      br(), br(),
      
      radioButtons(
        inputId = "exon_filter",
        label = "Select structural change:",
        choices = c("Downregulated exon", "Upregulated exon", "UTR change"),
        selected = "Downregulated exon"
      )
    ),
    
    tabsetPanel(
      id = "main_tabs",
      tabPanel("DTU Exploration", value = "dtu", isoformAnalysisUI("isoform")),
      tabPanel("Exon-Level", value = "exon", exonLevelUI("exon")),
      tabPanel("Summary", value = "summary", summaryStatsUI("summary"))
    )
  )
}