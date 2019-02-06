
# ==============================================================================
# This file contains functions which generate
# SingleCellExperiment objects.  Currently, these functions
# take in Seurat and Monocle objects and creates an object
# which is compatible with the CellTagViz shiny app. This
# contains expression data, and dimension reduction cell
# embeddings to visualize the data in two dimensions.
# ==============================================================================



#' Creates CellTagViz compatable data sets.
#'
#' \code{makeVizData} currently accepts Seurat and Monocle objects and returns
#' a SingleCellExperiment object for use with CellTagViz.
#'
#' This function accepts Seurat and Monocle objects. The data from these objects
#' is then combined and used to create a SingleCellExperiment object. This object
#' is compatable with CellTagViz and allows simple visualiztion of the analysis
#' performed with Seurat and Monocle.
#'
#' @param dataSets List This is a named list which contains the Seurat and Monocle
#' objects. The names for the list must be Seurat or Monocle corresponding to the
#' object.
#'
#' @return The function returns a SCE object with the data from Monocle and Seurat
#' objects combined.
#'
#' @examples
#'
#' \dontrun{
#'
#' dataSets <- list(seurat = seuratObj, monocle = monocleObj)
#'
#' sce <- makeVizData(dataSets)
#'
#' }
#'
#' @export

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

#' Returns vector of Cell Barcodes.
#'
#' \code{getCellBarcodes} is a helper function used to generate a list of Cell
#' Barcodes from Monocle and Seurat objects.
#'
#' This function accepts a Seurat and Monocle object, finds the union of the
#' cell barcodes in each object, and returns the character vector of the union
#' of the cell barcodes in both the Monocle and Seurat objects.
#'
#' @param seurat Seurat Object
#'
#' @param monocle Monocle Object
#'
#' @return THe function returns the union of the cell barcodes present in the
#' given Monocle and Seurat objects.
#'
#' @examples
#'
#' \dontrun{
#'
#' cellBCs <- getCellBarcodes(seuratObj, monocleObj)
#'
#' }
#'
#' @export

getCellBarcodes <- function(seurat, monocle) {

    barcodes <- union(seurat@cell.names, colnames(monocle))

}

#' Returns character vector of Feature Names.
#'
#' \code{getFeatures} is a helper function used to generate a list of the feature
#' names from Monocle and Seurat objects.
#'
#' This function accepts a Seurat and Monocle object, finds the union of the
#' feature names in each object, and returns the character vector of the union
#' of the feature names present in both the Monocle and Seurat objects.
#'
#' @param seurat Seurat Object
#'
#' @param monocle Monocle Object
#'
#' @return THe function returns the union of the feature names present in the
#' given Monocle and Seurat objects.
#'
#' @examples
#'
#' \dontrun{
#'
#' featureList <- getFeatures(seuratObj, monocleObj)
#'
#' }
#'
#' @export


getFeatures <- function(seurat, monocle) {

    features <- union(rownames(seurat@scale.data), rownames(monocle))

}

#' Initializes SCE object for Monocle and Seurat Data.
#'
#' \code{initializeSCE} is a helper function used to initialize an SCE object
#' with the correct cell barcodes and feature list.
#'
#' This function takes a list of cell barcodes and feature names and creates
#' an empty SCE object with the given cell barcodes and features. This empty
#' SCE object is then utilized to combine and add Monocle and Seurat data.
#'
#' @param cellBCs Character Vector of Cell Barcodes. Typically the output from
#' the function \code{getCellBarcodes}.
#'
#' @param features Character Vector of Feature Names. Typically the output
#' from the function \code{getFeatures}
#'
#' @return The function returns an empty SCE object with the given cell barcodes
#' and feature names. The SCE object can then be used to add Monocle and Seurat
#' data to.
#'
#' @examples
#'
#' \dontrun{
#'
#' cellBCs <- getCellBarcodes(seuratObj, monocleObj)
#'
#' featureList <- getFeatures(seuratObj, monocleObj)
#'
#' sce <- initializeSCE(cellBCs, featureList)
#'
#' }
#'
#' @export

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

        data2add <- base::cbind(data2add, bcMat)

    }

    if (length(missingFeatures) > 0) {

        featureMat <- Matrix::Matrix(nrow = length(missingFeatures),
                                     ncol = ncol(data2add))

        data2add <- base::rbind(data2add, featureMat)

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
