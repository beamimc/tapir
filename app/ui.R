library(shiny)
library(DT)
library(plotly)

ui <- fluidPage(
  titlePanel("Isoform Analysis Dashboard"),
  
  # Top section (gene selection, description, and GO table)
  fluidRow(
    column(5,
           wellPanel(
             selectizeInput("selected_gene", "Select a Gene:", 
                            choices = NULL, multiple = FALSE),
             
             tags$div(
               style = "border: 1px solid #ddd; padding: 5px; margin-top: 10px;",
               strong("Gene Description:"),
               textOutput("gene_description")
             ),
             
             tags$h4("Significant DTU"),
             DTOutput("gene_table")
           )
    ),
    column(6,
           tags$h4("Gene Ontology Terms"),
           DTOutput("go_table")
    )
  ),
  
  hr(),
  
  # Bottom section (side-by-side plots)
  fluidRow(
    column(6,
           tags$div(
             
             style = "width: 100%; height: 1000px; overflow-y: auto; border: 1px solid #ddd;",
             imageOutput("gene_plot", inline = TRUE)
             )
    ),
    column(3,
           plotlyOutput("barplot", width  = "100%")
    ),
    column(3,
           plotlyOutput("lineplot", width = "100%")
    )
  )
)
