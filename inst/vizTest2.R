

library(CellTagViz)

source("~/GitHub/working.Viz/R/dataTools2.R")

monoclePath <- "~/../Desktop/unsupervised timeline all data.RDS"

seuratPath <- "~/../Desktop/Warner_CCA_after_TSNE.RDS"

seuratPathV3 <- "~/../Desktop/integrated.subset.RDS"

monocleCDS <- readRDS(monoclePath)

seurat <- readRDS(seuratPath)

seuratv3 <- readRDS(seuratPathV3)

sceSeuratV2 <- makeSCESeurat(seuratObj = seurat)

sceSeuratV3 <- makeSCESeuratV3(seuratObj = seuratv3)

sceMonocle <- makeSCEMonocle(monocleObj = monocleCDS)

dim(sceSeuratV2)

dim(sceSeuratV3)

dim(sceMonocle)

SingleCellExperiment::reducedDimNames(sceSeuratV2)

SingleCellExperiment::reducedDimNames(sceSeuratV3)

SingleCellExperiment::reducedDimNames(sceMonocle)

SummarizedExperiment::assayNames(sceSeuratV2)

SummarizedExperiment::assayNames(sceSeuratV3)

SummarizedExperiment::assayNames(sceMonocle)

sceList <- makeVizDataNew(list("SeUratv2" = seurat, "SeuratV3" = seuratv3, "monocle" = monocleCDS))



comboList <- function(x) Reduce(c, x)
