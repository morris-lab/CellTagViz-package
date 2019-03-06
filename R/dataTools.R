
# ==============================================================================
# This file contains functions which generate SingleCellExperiment objects.
# Currently, these functions take in Seurat and Monocle objects and creates an
# object which is compatible with the CellTagViz shiny app. This contains
# expression data, and dimension reduction cell embeddings to visualize the data
# in two dimensions.
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
#' }
#'
#' @export

makeVizData <- function(dataSets) {

  names(dataSets) <- toupper(names(dataSets))

  seuratObject <- dataSets[["SEURAT2"]]

  monocleObject <- dataSets[["MONOCLE"]]

  seuratObjectV3 <- dataSets[["SEURAT3"]]

  cellBarcodes <-
    getCellBarcodes(seurat = seuratObject, monocle = monocleObject)

  featureList <-
    getFeatures(seurat = seuratObject, monocle = monocleObject)

  sceObject <-
    initializeSCE(cellBCs = cellBarcodes, features = featureList)

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
#' }
#'
#' @export

initializeSCE <- function(cellBCs, features) {
  cellDF <- data.frame(row.names = cellBCs)

  featureDF <- data.frame(row.names = features)

  sce <-
    SingleCellExperiment::SingleCellExperiment(
      colData = cellDF,
      rowData = featureDF
    )

  return(sce)
}


#' Function to add expression data to SCE object.
#'
#' \code{addExprData} is a helper function used to add expression matrices to an
#' SCE object.
#'
#' This function accepts a SCE object, an expression matrix, and a name for the
#' expression matrix. The expression matrix is typically Gene count data stored
#' in either a monocle or seurat object. This expression matrix is first
#' converted into a sparse matrix. Then cells and features that do not appear
#' in the matrix are added with 0 counts as values. Cells and Features may be
#' missing due to them being present in only one data set (Seurat or Monocle).
#' Once the missing cell barcodes and features have been added the sparse
#' expression matrix is then stored in the SCE object under the given name
#' (\code{datName}).
#'
#' @param sce Single Cell Experiment Object
#'
#' @param exprMat Matrix or object that can be coerced to a matrix containing
#' Feature count data.
#'
#' @param datName String Containing a unique ID for the expression matrix (\code{exprMat})
#'
#' @return The function returns a SCE object with the given expression matrix
#' stored under the name \code{datName} in the assays slot of the SCE object.
#'
#' @examples
#'
#' \dontrun{
#'
#' sce <- SingleCellExperiment::SingleCellExperiment()
#'
#' counts <- matrix(rnorm(100), 10, 10)
#'
#' sce <- addExprData(sce, counts, "Counts")
#' }
#'
addExprData <- function(sce, exprMat, datName) {
  exprMat <- Matrix::Matrix(exprMat, sparse = TRUE)

  missingGenes <-
    rownames(sce)[!rownames(sce) %in% rownames(exprMat)]

  missingCells <-
    colnames(sce)[!colnames(sce) %in% colnames(exprMat)]

  cellMat <-
    Matrix::Matrix(nrow = nrow(exprMat), ncol = length(missingCells))

  geneMat <-
    Matrix::Matrix(nrow = length(missingGenes), ncol = ncol(exprMat))

  exprMat <- cbind(exprMat, cellMat)

  exprMat <- rbind(exprMat, geneMat)

  sce@assays$data[[datName]] <- exprMat

  return(sce)
}

#' Creates sparse matrices of feature count data
#'
#' \code{makeDataMatrix} accepts a matrix like object of feature count data
#' and returns a sparse matrix of the feature counts.
#'
#' This is a helper function for adding feature counts to SCE objects.
#' A matrix of feature counts typically from a Seurat or Monocle object is
#' accepted as an argument. This matrix is then coerced into a sparse matrix and
#' missing cell barcodes and/or features are identified. The missing cell
#' barcodes and features are then added to the matrix with 0 values. The
#' completed matrix is then returned and can be stored as data in the assays
#' slot of an SCE object.
#'
#' @param data2add Matrix-like object of feature counts
#'
#' @param cellBCs Vector of all Cell Barcodes in the SCE
#'
#' @param features Vector of all features in the SCE
#'
#' @return The function returns a sparse matrix of feature counts to be stored
#' in a SCE object
#'
#' @examples
#'
#' barcodes <- c("BC1", "BC2")
#'
#' genes <- c("Apoa1", "Mettl7a1")
#'
#' exprMat <- matrix(rnorm(100), 3, 3)
#'
#' newMat <- makeDataMatrix(exprMat, barcodes, genes)
#' @export

makeDataMatrix <- function(data2add, cellBCs, features) {

  data2add <- Matrix::Matrix(data2add, sparse = TRUE)

  missingBCs <- cellBCs[!cellBCs %in% colnames(data2add)]

  missingFeatures <- features[!features %in% rownames(data2add)]

  if (length(missingBCs) > 0) {
    bcMat <-
      Matrix::Matrix(nrow = nrow(data2add), ncol = length(missingBCs))

    colnames(bcMat) <- missingBCs

    rownames(bcMat) <- rownames(data2add)

    data2add <- base::cbind(data2add, bcMat)
  }

  if (length(missingFeatures) > 0) {
    featureMat <-
      Matrix::Matrix(
        nrow = length(missingFeatures),
        ncol = ncol(data2add)
      )

    colnames(featureMat) <- colnames(data2add)

    rownames(featureMat) <- rownames(missingFeatures)

    data2add <- base::rbind(data2add, featureMat)
  }

  data2add <- as(data2add, "CsparseMatrix")

  return(data2add)
}

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

addSeuratData <- function(sce, seurat) {
  barcodeList <- colnames(sce)

  featureList <- rownames(sce)

  colnames(seurat@meta.data) <-
    paste0(colnames(seurat@meta.data), ".Seurat")

  assay(sce, "RawData.Seurat") <-
    makeDataMatrix(
      data2add = seurat@raw.data,
      cellBCs = barcodeList,
      features = featureList
    )

  assay(sce, "ScaleData.Seurat") <-
    makeDataMatrix(
      data2add = seurat@scale.data,
      cellBCs = barcodeList,
      features = featureList
    )

  sce@colData[rownames(seurat@meta.data), colnames(seurat@meta.data)] <-
    seurat@meta.data[rownames(seurat@meta.data), colnames(seurat@meta.data)]

  cellEmbeddings <- lapply(names(seurat@dr), function(id) {
    seurat@dr[[id]]@cell.embeddings
  })

  names(cellEmbeddings) <- paste0(names(seurat@dr), ".Seurat")

  reducedDims(sce) <- cellEmbeddings

  return(sce)
}


#' Add monocle data to a SCE object
#'
#' \code{addMonocleData} is a function used to transfer data from a monocle
#' object to a SCE object.
#'
#' This function is the main function used to transfer data from a given
#' monocle object into a given SCE object. The expression count data from the
#' monocle object is added to the assays slot of the SCE object. Currently the
#' only cell embeddings transferred from the monocle object are the 'Pseudotime'
#' embeddings from the reducedDimS slot of the monocle object. If other dimension
#' reduction methods were used they will not be included using this function.
#' The meta data generated by monocle is added to the colData slot of the SCE
#' object. Monocle also generates feature level meta data which is currently not
#' transferred using this function. All names are suffixed with '.Monocle' before
#' being stored in the SCE object.
#'
#' @param sce Single Cell Experiment in which data will be stored
#'
#' @param monocleObj Monocle CellDataSet Contains data to transfer to SCE object
#'
#' @return This function returns a SCE object containing the data from the given
#' monocle object.
#'
#' @examples
#'
#' \dontrun{
#'
#' sce <- SingleCellExperiment()
#'
#' sceData <- addMonocleData(sce, monocleObj)
#' }
#'
#' @export

addMonocleData <- function(sce, monocleObj) {
  barcodeList <- colnames(sce)

  featureList <- rownames(sce)


  pseudoDims <- t(monocleObj@reducedDimS)

  missingBCs <- barcodeList[!barcodeList %in% rownames(pseudoDims)]

  missingDat <-
    matrix(
      data = 0,
      nrow = length(missingBCs),
      ncol = 2
    )

  rownames(missingDat) <- missingBCs

  pseudoAdd <- rbind(pseudoDims, missingDat)

  colnames(pseudoAdd) <- c("Dim.1", "Dim.2")

  reducedDim(sce, type = "Pseudotime.Monocle") <- pseudoAdd

  assay(sce, "Counts.Monocle") <-
    makeDataMatrix(
      data2add = monocleObj@assayData$exprs,
      cellBCs = barcodeList,
      features = featureList
    )

  names(monocleObj@phenoData@data) <-
    paste0(names(monocleObj@phenoData@data), ".Monocle")

  sce@colData[rownames(monocleObj@phenoData@data), colnames(monocleObj@phenoData@data)] <-
    monocleObj@phenoData@data[
      rownames(monocleObj@phenoData@data),
      colnames(monocleObj@phenoData@data)
    ]

  return(sce)
}


addMonocleDataNew <- function(monocleObj){

  featureDat <- Biobase::fData(monocleObj)

  phenoDat <- Biobase::pData(monocleObj)

  pseudoDims <- t(monocleObj@reducedDimS)

  monCounts <- Biobase::exprs(monocleObj)

  colnames(featureDat) <- paste0(colnames(featureDat), ".Monocle")

  colnames(phenoDat) <- paste0(colnames(phenoDat), ".Monocle")

  monSCE <- SingleCellExperiment::SingleCellExperiment(colData = phenoDat,
                                             rowData = featureDat,
                                             reducedDims = S4Vectors::SimpleList(Pseudotime.Monocle = pseudoDims))

  SummarizedExperiment::assay(monSCE, "Monocle.Counts") <- monCounts

  return(monSCE)

}



#' @importClassesFrom S4Vectors SimpleList
#'
#'



addseuratV3 <- function(seuratObj){

  colDat <- seuratObj@meta.data

  rowDat <- data.frame(rownames(seuratObj@assays$RNA@counts))

  rownames(rowDat) <- rowDat[[1]]

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

addMissingFeatures <- function(dataMat, features){

  missingFeatures <- features[! features %in% rownames(dataMat)]

  featureMat <- Matrix::Matrix(data = 0, nrow = length(missingFeatures), ncol = ncol(dataMat))

  rownames(featureMat) <- missingFeatures

  colnames(featureMat) <- colnames(dataMat)

  dataMat <- rbind(dataMat, featureMat)

  return(dataMat)

}
