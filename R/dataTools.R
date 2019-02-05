
# ==============================================================================
# This file contains functions which generate
# SingleCellExperiment objects.  Currently, these functions
# take in Seurat and Monocle objects and creates an object
# which is compatible with the CellTagViz shiny app. This
# contains expression data, and dimension reduction cell
# embeddings to visualize the data in two dimensions.
# ==============================================================================


makeVizData <- function(dataSets) {

    names(dataSets) <- toupper(names(dataSets))

    seuratObject <- dataSets[["SEURAT"]]

    monocleObject <- dataSets[["MONOCLE"]]

    cellBarcodes <- getCellBarcodes(seurat = seuratObject,
                                   monocle = monocleObject)

    featureList <- getFeatures(seurat = seuratObject,
                              monocle = monocleObject)

    sceObject <- initializeSCE(cellBCs = cellBarcodes,
                              features = featureList)

    sceObject <- addSeuratData(sceObject, seuratObject)

    sceObject <- addMonocleData(sceObject, monocleObject)

    return(sceObject)
}

getCellBarcodes <- function(seurat, monocle) {

    barcodes <- union(seurat@cell.names, colnames(monocle))

}

getFeatures <- function(seurat, monocle) {

    features <- union(rownames(seurat@scale.data), rownames(monocle))

}

initializeSCE <- function(cellBCs, features) {

    cellDF <- data.frame(row.names = cellBCs)

    featureDF <- data.frame(row.names = features)

    sce <- SingleCellExperiment::SingleCellExperiment(colData = cellDF,
                                                      rowData = featureDF)

    return(sce)
}

addExprData <- function(sce, exprMat, datName) {

    exprMat <- Matrix::Matrix(exprMat, sparse = TRUE)

    missingGenes <- rownames(sce)[!rownames(sce) %in% rownames(exprMat)]

    missingCells <- colnames(sce)[!colnames(sce) %in% colnames(exprMat)]

    cellMat <- Matrix::Matrix(nrow = nrow(exprMat), ncol = length(missingCells))

    geneMat <- Matrix::Matrix(nrow = length(missingGenes), ncol = ncol(exprMat))

    exprMat <- cbind(exprMat, cellMat)

    exprMat <- rbind(exprMat, geneMat)

    sce@assays$data[[datName]] <- exprMat

    return(sce)
}


makeDataMatrix <- function(data2add, cellBCs, features) {

    data2add <- Matrix::Matrix(data2add, sparse = TRUE)

    missingBCs <- cellBCs[!cellBCs %in% colnames(data2add)]

    missingFeatures <- features[!features %in% rownames(data2add)]

    if (length(missingBCs) > 0) {

        bcMat <- Matrix::Matrix(nrow = nrow(data2add),
                                ncol = length(missingBCs))

        data2add <- cbind(data2add, bcMat)

    }

    if (length(missingFeatures) > 0) {

        featureMat <- Matrix::Matrix(nrow = length(missingFeatures),
                                     ncol = ncol(data2add))

        data2add <- rbind(data2add, featureMat)

    }

    return(data2add)

}

addSeuratData <- function(sce, seurat) {

    barcodeList <- colnames(sce)

    featureList <- rownames(sce)

    colnames(seurat@meta.data) <- paste0(colnames(seurat@meta.data),
        ".Seurat")

    sce@assays$data[["RawData.Seurat"]] <- makeDataMatrix(dat = seurat@raw.data,
                                                      cellBCs = barcodeList,
                                                     features = featureList)

    sce@assays$data[["ScaleData.Seurat"]] <- makeDataMatrix(dat = seurat@scale.data,
                                                        cellBCs = barcodeList,
                                                       features = featureList)

    sce@colData[rownames(seurat@meta.data), colnames(seurat@meta.data)] <- seurat@meta.data[rownames(seurat@meta.data),
        colnames(seurat@meta.data)]

    cellEmbeddings <- lapply(names(seurat@dr), function(id) {
        seurat@dr[[id]]@cell.embeddings
    })

    names(cellEmbeddings) <- paste0(names(seurat@dr), ".Seurat")

    sce@reducedDims@listData <- cellEmbeddings

    return(sce)

}


addMonocleData <- function(sce, monocleObj) {

    barcodeList <- colnames(sce)

    featureList <- rownames(sce)

    sce@reducedDims@listData[["Pseudotime.Monocle"]] <- Matrix::Matrix(data = t(monocleObj@reducedDimS),
        sparse = TRUE)

    sce@assays$data[["Counts.Monocle"]] <- makeDataMatrix(dat = monocleObj@assayData$exprs,
                                                      cellBCs = barcodeList,
                                                     features = featureList)

    names(monocleObj@phenoData@data) <- paste0(names(monocleObj@phenoData@data),
        ".Monocle")

    sce@colData[rownames(monocleObj@phenoData@data), colnames(monocleObj@phenoData@data)] <- monocleObj@phenoData@data[rownames(monocleObj@phenoData@data),
        colnames(monocleObj@phenoData@data)]

    return(sce)

}
