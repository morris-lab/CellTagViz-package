
library(tidyverse)
library(shiny)
library(shinythemes)
library(ggplot2)
library(SingleCellExperiment)
library(Seurat)

library(CellTagViz)


monoclePath <- "~/../Desktop/unsupervised timeline all data.RDS"

#monoclePath <- "../../../Desktop/unsupervised timeline all data.RDS"

#seuratPath <- "../../../Desktop/Warner_CCA_after_TSNE.RDS"

seuratPath <- "~/../Desktop/Warner_CCA_after_TSNE.RDS"

monocleCDS <- readRDS(monoclePath)

seuratObj <- readRDS(seuratPath)


dataList <- list(seurat = seuratObj,
  monocle = monocleCDS)


sce <- CellTagViz::makeVizData(dataSets = dataList)



source("R/shinyModules.R")

source("R/shinyModulesServer.R")

source("R/dataTools.R")

source("R/sceTools.R")

ui <- navbarPage(

  theme = shinythemes::shinytheme("cosmo"),
  title = "CellTagViz",

  welcomePanelUI("welcome"),
  plotsPanelUI("plots", inputData = sce),
  dataPanelUI("data", inputData = sce),
  peoplePanelUI("people")

)


server <- function(input, output, session){

  userInput <- reactive({

    return(input)

  }) %>% debounce(2000)

  callModule(createPlot, id = "plots", plotOptions = userInput(), inputSCE = sce)

}




shinyApp(ui = ui, server = server)


