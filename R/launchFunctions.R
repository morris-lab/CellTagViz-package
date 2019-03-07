


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
#'
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

  appPath <- system.file("shinyApps", "testApp.R", package = "CellTagViz")

  ui <- NULL

  server <- NULL

  source(appPath, local = TRUE)

  serverEnv <- environment(server)

  serverEnv$sce <- vizData

  app <- shiny::shinyApp(ui = ui, server = server)

  shiny::runApp(app)

}



