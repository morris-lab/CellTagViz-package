
#' Add Seurat object data to an SCE object
#'
#' \code{addSeuratData} is a function used to parse Seurat objects and add the
#' data contained in seurat objects to a SCE object.
#'
#' This is the main function used to add data to a SCE object from a seurat object.
#' This function identifies the cell barcodes and features listed in the seurat object.
#' The raw and scaled expression data from the seurat object are added to the assays slot
#' of the SCE object. These data are processed to make sure all of the cells
#' and features from both the monocle object and seurat object are represented.
#' Furthermore, the meta data from the seurat object is added to the colData of
#' the SCE object. Each variable name in the seurat meta data is given the
#' suffix '.Seurat' to identify its origins once stored in the SCE object.
#' The cell embeddings for each dimension reduction mehtod used in seurat are
#' stored in the SCE object. Once again the names for each set of cell embeddings
#' is suffixed with '.Seurat' to identify the origins of the embeddings.
#'
#' @param sce Single Cell Experiment Object used to store data
#'
#' @param seurat Seurat Object from which data will be transferred from.
#'
#' @return The function returns a SingleCellExperiment object which stores select
#' data from the given seurat object.
#'
#' @examples
#'
#' \dontrun{
#'
#' sce <- SingleCellExperiment()
#'
#' sceData <- addSeuratData(sce, seuratObj)
#' }
#'
#' @export

makeSCESeurat <- function(seuratObj) {

  colDat <- seuratObj@meta.data

  colnames(colDat) <- paste0(colnames(colDat), ".SeuratV2")

  rowDat <- data.frame(rownames(seuratObj@raw.data))

  rownames(rowDat) <- rowDat[[1]]

  colnames(rowDat) <- "GeneShortName.SeuratV2"

  assayList <- list("RawData.SeuratV2" = seuratObj@raw.data, "ScaleData.SeuratV2" = seuratObj@scale.data, "Data.SeuratV2" = seuratObj@data)

  cellEmbeddings <- lapply(names(seuratObj@dr), function(id) {
    seuratObj@dr[[id]]@cell.embeddings
  })

  names(cellEmbeddings) <- paste0(names(seuratObj@dr), ".SeuratV2")

  cellEmbeddings <- as(cellEmbeddings, "SimpleList")

  sce <- SingleCellExperiment::SingleCellExperiment(colData = colDat, rowData = rowDat, assays = assayList)

  SingleCellExperiment::reducedDims(sce) <- cellEmbeddings

  return(sce)
}


#' Title
#'
#' @param seuratObj
#'
#' @return
#' @export
#'
#' @examples
#'

makeSCESeuratV3 <- function(seuratObj){

  colDat <- seuratObj@meta.data

  colnames(colDat) <- paste0(colnames(colDat), ".SeuratV3")

  rowDat <- data.frame(rownames(seuratObj@assays$RNA@counts))

  rownames(rowDat) <- rowDat[[1]]

  colnames(rowDat) <- "GeneShortName.SeuratV3"

  exprList <- list()

  slots <- c("counts", "data", "scale.data")

  for(i in names(seuratObj@assays)){

    for(j in slots){

      id <- paste0(i, ".", j, ".SeuratV3")

      temp <- GetAssayData(seuratObj, assay = i, slot = j)

      if(length(temp) > 0){

        if(nrow(temp) < nrow(rowDat)){

          datAssay <- addMissingFeatures(temp, rowDat[[1]])

          datAssay <- datAssay[rowDat[[1]], ]

          exprList[[id]] <- datAssay

        } else(

          exprList[[id]] <- temp[rowDat[[1]], ]

        )
      }
    }
  }

  redMethods <- names(seuratObj@reductions)

  names(redMethods) <- paste0(redMethods, ".SeuratV3")

  cellEmbeddings <- lapply(redMethods, function(id){

    Seurat::Embeddings(seuratObj, reduction = id)

  })

  cellEmbeddings <- as(cellEmbeddings, "SimpleList")

  sce <- SingleCellExperiment::SingleCellExperiment(colData = colDat, rowData = rowDat, assays = exprList)

  SingleCellExperiment::reducedDims(sce) <- cellEmbeddings

  return(sce)

}















#' Title
#'
#' @param monocleObj
#'
#' @return
#' @export
#'
#' @examples
#'

makeSCEMonocle <- function(monocleObj){

  featureDat <- Biobase::fData(monocleObj)

  geneCol <- grep(pattern = "gene_short_name", x = colnames(featureDat), fixed = TRUE)

  colnames(featureDat)[geneCol] <- "GeneShortName"

  phenoDat <- Biobase::pData(monocleObj)

  pseudoDims <- t(monocleObj@reducedDimS)

  monCounts <- Biobase::exprs(monocleObj)

  colnames(featureDat) <- paste0(colnames(featureDat), ".Monocle")

  colnames(phenoDat) <- paste0(colnames(phenoDat), ".Monocle")

  monSCE <- SingleCellExperiment::SingleCellExperiment(colData = phenoDat,
    rowData = featureDat)

  SingleCellExperiment::reducedDims(monSCE) <- S4Vectors::SimpleList("Pseudotime.Monocle" = pseudoDims)

  SummarizedExperiment::assay(monSCE, "Monocle.Counts") <- monCounts

  return(monSCE)

}






#' Title
#'
#' @param dataMat
#' @param features
#'
#' @return
#' @export
#'
#' @examples
#'

addMissingFeatures <- function(dataMat, features){

  missingFeatures <- features[! features %in% rownames(dataMat)]

  featureMat <- Matrix::Matrix(data = 0, nrow = length(missingFeatures), ncol = ncol(dataMat))

  rownames(featureMat) <- missingFeatures

  colnames(featureMat) <- colnames(dataMat)

  dataMat <- rbind(dataMat, featureMat)

  return(dataMat)

}


#' Title
#'
#' @param dataSets
#'
#' @return
#' @export
#'
#' @examples
#'

makeVizDataNew <- function(dataSets) {

  names(dataSets) <- toupper(names(dataSets))

  sceList <- sapply(names(dataSets), FUN = function(x){createSCE(objName = x, scData = dataSets[[x]])}, USE.NAMES = TRUE, simplify = FALSE)

  dims <- sapply(sceList, ncol)

  small <- names(dims)[which.min(dims)]

  cellBCs <- colnames(sceList[[small]])

  sceList <- lapply(sceList, FUN = function(x){x[, cellBCs]})

  rowDat <- sapply(sceList, FUN = function(x){

     SingleCellExperiment::rowData(x)

    })

  rowDat <- S4Vectors::cbind.DataFrame(unlist(rowDat))

  colDat <- sapply(sceList, FUN = function(x){

    SingleCellExperiment::colData(x)[cellBCs, , drop = FALSE]

  })

  colDat <- S4Vectors::cbind.DataFrame(unlist(colDat))

  metaNames <- sapply(colnames(colDat), FUN = function(x){

    x <- unlist(strsplit(x, split = ".", fixed = TRUE))

    x <- paste0(x[-1], collapse = ".")

  })

  colnames(colDat) <- metaNames

  featureNames <- sapply(colnames(rowDat), function(x){

    x <- unlist(strsplit(x, split = ".", fixed = TRUE))

    x <- paste0(x[-1], collapse = ".")

  })


  sceSeuratV2 <- sceSeuratV2[rownames(rowDat), rownames(colDat)]

  sceSeuratV3 <- sceSeuratV3[rownames(rowDat), rownames(colDat)]

  sceMonocle <- sceMonocle[rownames(rowDat), rownames(colDat)]


  SummarizedExperiment::rowData(sceSeuratV2) <- rowDat

  SummarizedExperiment::rowData(sceSeuratV3) <- rowDat

  SummarizedExperiment::rowData(sceMonocle) <- rowDat

  SummarizedExperiment::colData(sceSeuratV2) <- colDat

  SummarizedExperiment::colData(sceSeuratV3) <- colDat

  SummarizedExperiment::colData(sceMonocle) <- colDat


  comboSCE <- cbind(sceSeuratV2, sceSeuratV3, sceMonocle)






  reductionDat <- sapply(sceList, FUN = function(x){

    SingleCellExperiment::reducedDims(x)

  })

  reductionDat <- Reduce(c, reductionDat)

  comboSCE <- SingleCellExperiment::SingleCellExperiment(colDat = colDat, rowDat = rowDat)

  SingleCellExperiment::reducedDims(comboSCE) <- reductionDat


  #
  # seuratObject <- dataSets[["SEURAT2"]]
  #
  # monocleObject <- dataSets[["MONOCLE"]]
  #
  # seuratObjectV3 <- dataSets[["SEURAT3"]]
  #
  # cellBarcodes <-
  #   getCellBarcodes(seurat = seuratObject, monocle = monocleObject)
  #
  # featureList <-
  #   getFeatures(seurat = seuratObject, monocle = monocleObject)
  #
  # sceObject <-
  #   initializeSCE(cellBCs = cellBarcodes, features = featureList)
  #
  # sceObject <- addSeuratData(sceObject, seuratObject)
  #
  # sceObject <- addMonocleData(sceObject, monocleObject)

  return(sceList)
}


createSCE <- function(scData, objName){

  switch(EXPR = objName, "SEURATV2" = makeSCESeurat(scData), "SEURATV3" = makeSCESeuratV3(scData), "MONOCLE" = makeSCEMonocle(scData))

}




combineSCE <- function(dataSets){

  names(dataSets) <- toupper(names(dataSets))

  sceList <- sapply(names(dataSets), FUN = function(x){createSCE(objName = x, scData = dataSets[[x]])}, USE.NAMES = TRUE, simplify = FALSE)

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


}
