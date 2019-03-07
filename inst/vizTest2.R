

library(CellTagViz)

library(SingleCellExperiment)

library(shiny)


source("~/GitHub/working.Viz/R/dataTools2.R")

monoclePath <- "~/../Desktop/unsupervised timeline all data.RDS"

seuratPath <- "~/../Desktop/Warner_CCA_after_TSNE.RDS"

seuratPathV3 <- "~/../Desktop/integrated.subset.RDS"

monocleCDS <- readRDS(monoclePath)

seurat <- readRDS(seuratPath)

seuratv3 <- readRDS(seuratPathV3)

dataSets <- list("seuratv2" = seurat, "seuratv3" = seuratv3, "monocle" = monocleCDS)

sce <- combineSCE(dataSets = dataSets)


ui <- navbarPage(

  theme = shinythemes::shinytheme("yeti"),
  title = "CellTagViz",

  plotsPanelMinimalUI("plots", inputData = sce),
  dataPanelUI("data", inputData = sce)

)


server <- function(input, output, session){

  userInput <- reactive({

    return(input)

  }) %>% debounce(1000)

  callModule(createPlot, id = "plots", plotOptions = userInput(), inputSCE = sce)

}




shinyApp(ui = ui, server = server)

