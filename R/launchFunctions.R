


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

TestApp <- function(){

  appPath <- system.file("shinyApps", "testApp.R", package = "CellTagViz")

  shiny::shinyAppFile(appFile = appPath)

}



