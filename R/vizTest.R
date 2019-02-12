


library(shiny)


ui <- navbarPage(

  theme = shinythemes::shinytheme("cosmo"),

  title = "CellTagViz",

  welcomePanelUI("welcome"),

  plotsPanelUI("plots"),

  dataPanelUI("data"),

  peoplePanelUI("people")

)


server <- function(input, output){


}

shinyApp(ui = ui, server = server)
