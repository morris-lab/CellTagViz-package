

dim(sceSeuratV2)

dim(sceSeuratV3)

dim(sceMonocle)

sceList <- list(SeuratV2 = sceSeuratV2, SeuratV3 = sceSeuratV3, Monocle = sceMonocle)

cellCounts <- sapply(sceList, ncol, simplify = TRUE)

cellBCs <- colnames(sceList[[which.min(cellCounts)]])

featureList <- rownames(sceList[[which.min(cellCounts)]])

sceList <- sapply(sceList, function(sce){

  sce[featureList, cellBCs]

})

rDat <- sapply(sceList, simplify = TRUE, function(sce){

  rowData(sce)

})

rDat <- Reduce(c, rDat)

cDat <- sapply(sceList, simplify = TRUE, function(sce){

  colData(sce)

})

cDat <- Reduce(c, cDat)

redDat <- sapply(sceList, simplify = TRUE, function(sce){

  reducedDims(sce)

})

redDat <- Reduce(c, redDat)

datAssay <- sapply(sceList, function(sce){

  assays(sce)

})


datAssay <- Reduce(c, datAssay)


finalSCE <- SingleCellExperiment::SingleCellExperiment(colData = cDat, rowData = rDat, assays = datAssay, reducedDims = redDat)


dim(finalSCE)

assayNames(finalSCE)

reducedDimNames(finalSCE)




