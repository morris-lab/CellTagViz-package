

library(tidyverse)

library(SingleCellExperiment)

library(Seurat)

dataTools <- "../working.Viz/R/dataTools.R"

source(file = dataTools)

monoclePath <- "../../../Desktop/unsupervised timeline all data.RDS"

seuratPath <- "../../../Desktop/Warner_CCA_after_TSNE.RDS"

monocleCDS <- readRDS(monoclePath)

seuratObj <- readRDS(seuratPath)


dataList <- list(seurat = seuratObj,
                monocle = monocleCDS)


sce <- makeVizData(dataSets = dataList)





plotData <- as_tibble(reducedDim(sce, "cca.aligned.Seurat"))

plotData <- reducedDim(sce, "cca.aligned.Seurat")

cellBCs <- rownames(plotData)

plotData <- as_tibble(plotData, rownames = NA)

plotData$CellBCs <- cellBCs


metaData <- colData(sce)

metaData$CellBCs <- rownames(metaData)


metaData <- as_tibble(metaData, rownames = NA)



ggplot(data = plotData) + geom_point(aes(x = ACC1, y = ACC2))


generatePlotData <- function(sce, metaData, geneChoice, redMethod){

  plotData <- as_tibble(reducedDim(sce, redMethod))

  plotData[, colnames(sce@colData)] <- sce@colData[rownames(plotData), colnames(sce@colData)]

}

allData <- plyr::join(plotData, metaData)

ggplot(data = allData) + geom_point(aes(x = ACC1, y = ACC2, color = as.factor(State.Monocle)))


nba <- ggplot(data = allData) + geom_point(aes(x = ACC1, y = ACC2, color = as.factor(State.Monocle)))

plotly::ggplotly(nba)

plotly::orca(nba, file = "cca.test.svg", format = "svg")




testRun <- function(){

 shiny::shinyAppFile("~/working/working.Viz/R/app.R")

}
