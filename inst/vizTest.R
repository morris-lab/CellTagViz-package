


library(shiny)


source("./shinyModules.R")

ui <- navbarPage(

  theme = shinythemes::shinytheme("cosmo"),
  title = "CellTagViz",

  welcomePanelUI("welcome"),
  plotsPanelUI("plots"),
  dataPanelUI("data"),
  peoplePanelUI("people")

)


server <- function(input, output){

  output$input_out <- renderPrint({

    str(sapply(names(input), function(id){input[[id]]}, simplify = FALSE))

  })
}


shinyApp(ui = ui, server = server)
