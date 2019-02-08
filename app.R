

library(shiny)

source("./linked_scatter.R")

ui <- navbarPage("CellTagViz",
  theme = shinythemes::shinytheme("cosmo"),
  welcomePanel("scatters"),
  plotsPanel("scatters"),
  dataPanel("scatters"),
  peoplePanel("scatters")
)

server <- function(input, output, session) {
  df <- callModule(linkedScatter, "scatters", reactive(mpg),
                   left = reactive(c("cty", "hwy")),
                   right = reactive(c("drv", "hwy"))
  )

  output$summary <- renderText({
    sprintf("%d observation(s) selected", nrow(dplyr::filter(df(), selected_)))
  })
}

shinyApp(ui, server)

