
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
#'
#' @param seuratObj Seurat Object from which data will be transferred from.
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

  cellEmbeddings <- methods::as(cellEmbeddings, "SimpleList")

  sce <- SingleCellExperiment::SingleCellExperiment(colData = colDat, rowData = rowDat, assays = assayList)

  SingleCellExperiment::reducedDims(sce) <- cellEmbeddings

  return(sce)
}


#' Function to create SingleCellExperiment Object.
#'
#' This function takes a Seurat V3 object as input and creates a
#' SingleCellExperiment Object using the data stored in the seurat object.
#'
#' @param seuratObj SeuratV3 object
#'
#' @return A SingleCellExperiment Object
#'
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

      temp <- Seurat::GetAssayData(seuratObj, assay = i, slot = j)

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

  cellEmbeddings <- methods::as(cellEmbeddings, "SimpleList")

  sce <- SingleCellExperiment::SingleCellExperiment(colData = colDat, rowData = rowDat, assays = exprList)

  SingleCellExperiment::reducedDims(sce) <- cellEmbeddings

  return(sce)

}















#' Function to create SingleCellExperiment object using a monocle CDS.
#'
#' This function takes a monocle CDS object as input and converts the
#' CellDataSet object into a SingleCellExperiment object.
#'
#' @param monocleObj CellDataSet Object with single-cell Data
#'
#' @return A SingleCellExperiment object
#'
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






#' Function to add missing features to a matrix.
#'
#' @param dataMat Matrix
#' @param features List of all features expected
#'
#' @return Matrix with missing features added
#'

addMissingFeatures <- function(dataMat, features){

  missingFeatures <- features[! features %in% rownames(dataMat)]

  featureMat <- Matrix::Matrix(data = 0, nrow = length(missingFeatures), ncol = ncol(dataMat))

  rownames(featureMat) <- missingFeatures

  colnames(featureMat) <- colnames(dataMat)

  dataMat <- rbind(dataMat, featureMat)

  return(dataMat)

}


#' Switch function used to create SingleCellExperiment Objects
#'
#' This function takes in a data object from monocle, seuratv2, or seuratv3 as
#' input as well as the name of the given object. The function used to convert
#' the given object into a SingleCellExperiment object is determined using the
#' given name of the data object.
#'
#' @param scData Object containing Single-cell data
#' @param objName Unique name for object
#'
#' @return A SingleCellExperiment Object
#'

createSCE <- function(scData, objName){

  switch(EXPR = objName, "SEURATV2" = makeSCESeurat(scData), "SEURATV3" = makeSCESeuratV3(scData), "MONOCLE" = makeSCEMonocle(scData))

}




#' Function used to combine multiple SingleCellExperiment objects.
#'
#' This function takes a named list of data objects from monocle, seuratv2, or
#' seuratv3. It is important that the list is named as the names are used to
#' determine how to handle each object.
#'
#' The function accepts a list of single-cell data objects. This list is then
#' converted into a list of SingleCellExperiment objects.
#'
#' @param dataSets Named List
#'
#' @return SingleCellExperiment
#' @export
#'

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

    SingleCellExperiment::rowData(sce)

  })

  rDat <- Reduce(c, rDat)

  cDat <- sapply(sceList, simplify = TRUE, function(sce){

    SingleCellExperiment::colData(sce)

  })

  cDat <- Reduce(c, cDat)

  redDat <- sapply(sceList, simplify = TRUE, function(sce){

    SingleCellExperiment::reducedDims(sce)

  })

  redDat <- Reduce(c, redDat)

  datAssay <- sapply(sceList, function(sce){

    SummarizedExperiment::assays(sce)

  })


  datAssay <- Reduce(c, datAssay)


  finalSCE <- SingleCellExperiment::SingleCellExperiment(colData = cDat, rowData = rDat, assays = datAssay, reducedDims = redDat)


}


#' @importClassesFrom S4Vectors SimpleList
#'

NULL
