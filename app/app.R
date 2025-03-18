# Source the UI and server definitions
source("ui.R")
source("server.R")
source("global.R")


# define app 
shinyApp(ui, server)
