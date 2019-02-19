


library(shiny)




source("R/shinyModules.R")

source("R/shinyModulesServer.R")

ui <- navbarPage(

  theme = shinythemes::shinytheme("cosmo"),
  title = "CellTagViz",

  welcomePanelUI("welcome"),
  plotsPanelUI("plots"),
  dataPanelUI("data"),
  peoplePanelUI("people")

)


server <- function(input, output, session){

  userInput <- reactive({

    return(input)

  })

  #callModule(plotOpts, id = "plots", brent = choices)

  callModule(createPlot, id = "plots", plotOptions = userInput())

}




shinyApp(ui = ui, server = server)


