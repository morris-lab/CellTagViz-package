

devtools::install_github(

  repo = "babiddy/working.Viz",

  auth_token = "9c4b394c583d882d894e93b69dfc7df1d11b32a9"
)


vizFile <- system.file("vizTest.R", package = "CellTagViz")

shiny::shinyAppFile(appFile = vizFile)
