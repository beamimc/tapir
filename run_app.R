# launch_app.R
library(here)
# library(GenomicFeatures)


see <- readRDS(here::here("data", "glinos_saturn_dtu.rds"))
# txdbb <- loadDb(here("data","flair_filter_transcripts.sqlite"))
exonss <- readRDS(here::here("data", "glinos_exons.rds"))


app_dir <-file.path(here::here(), "app")
source(file.path(app_dir, "app.R"))
shiny::runApp(app(se= see, 
                  exons=exonss, 
                  app_dir = app_dir))