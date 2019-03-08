


#' Title
#'
#' @return Yes
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' yes
#'
#' }

CellTagViz <- function(){

  appPath <- system.file("shinyApps", "vizTest.R", package = "CellTagViz")

  shiny::shinyAppFile(appFile = appPath)

}


#' Title
#' @param vizData sce object
#' @return Yes
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' Yes
#'
#' }

TestApp <- function(vizData){

  # require(shiny)
  #
  # ui <- navbarPage(
  #
  #   theme = shinythemes::shinytheme("yeti"),
  #   title = "CellTagViz",
  #
  #   plotsPanelMinimalUI("plots", inputData = sce),
  #   dataPanelUI("data", inputData = sce)
  #
  # )
  #
  #
  # server <- function(input, output, session){
  #
  #   userInput <- reactive({
  #
  #     return(input)
  #
  #   }) %>% debounce(2000)
  #
  #   callModule(createPlot, id = "plots", plotOptions = userInput(), inputSCE = sce)
  #
  # }
  #
  # app <- shiny::shinyApp(ui = ui, server = server)

  .GlobalEnv$sce <- vizData

  appPath <- system.file("shinyApps", "testApp.R", package = "CellTagViz")

  shiny::shinyAppFile(appPath)

}



