#' Access the cell embeddings of SCE object.
#'
#' \code{getEmbeddings} returns the coordinates of cells after dimension
#' reduction.
#'
#' This is a function used to subset SCE objects. The output is used to generate
#' a data frame which is compatible with ggplot. The cell embeddings are then
#' combined with metadata to generate visualizations with ggplot2.
#'
#' @param sce SingleCellExperiment Object
#'
#' @param redMethod String Name of dimension reduction method to visualize.
#' Must be one of methods listed in \code{names(sce@@reducedDim)}
#'
#' @param cells Character Vector of cells used to subset data. (Optional)
#'
#' @return The function returns the cell embeddings stored in the given slot of
#' the SingleCellExperiment object.
#'
#'
#' @examples
#'
#' sce <- SingleCellExperiment::SingleCellExperiment()
#'
#' sce@@reducedDims@@listData$tSNE <- matrix(data = rnorm(100),
#'                                           nrow = 10,
#'                                           ncol = 10)
#'
#' cellEmbeddings <- CellTagViz:::getEmbeddings(sce, "tSNE")
#'

getEmbeddings <- function(sce, redMethod, cells = FALSE){

  embeddings <- SingleCellExperiment::reducedDim(sce, redMethod)

  if(cells){

    embeddings <- embeddings[cells, ]
  }

  return(embeddings)
}


#' Function used to access meta data (colData) of an SCE object.
#'
#' \code{getMetaData} returns the meta data variable \emph{varName} of the SCE
#' object \emph{sce} for the given vector of cells, \emph{cells}.
#'
#' @param sce SingleCellExperiment Object
#'
#' @param varName String Name of the variable to return from meta/col data.
#'
#' @param cells Character Vector of Cell Barcodes used to subset data.
#'
#' @return This function returns given meta data (\emph{varName}) for the given
#' cells (\emph{cells}) of a SCE object (\emph{sce}).
#'
#' @examples
#'
#' metaData <- sample(letters, 15)
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(colData = list(Letters = metaData))
#'
#' letters <- CellTagViz:::getMetaData(sce, "Letters")
#'

getMetaData <- function(sce, varName, cells = FALSE){

  metaData <- SingleCellExperiment::colData(sce)[varName]

  if(cells){

    metaData <- metaData[cells, ]
  }

  return(metaData)

}


#' Generates ggplot2 compatible data frames.
#'
#' \code{makePlotData} creates data frames for use with ggplot2 visualizations.
#' This is accomplished by calling two functions \code{getEmbeddings} and
#' \code{getMetaData}. The objects returned by both of these functions are
#' then merged into one single data frame. This combined data frame contains the
#' coordinates of each cell from the given dimension reduction method along with
#' the given meta data which can then be used for grouping and coloring purposes.
#'
#' @param sce SingleCellExperiment Object
#'
#' @param redMethod String Name of the dimension reduction method to visualize.
#'
#' @param metaVar String Name of the column to return from meta/col Data.
#'
#' @param cells Character Vector of Cell Barcodes used to subset data.
#'
#' @return This function returns a data frame which can be used to plot and color
#' cells from an SCE object using ggplot2.
#'
#' @examples
#'
#' metaData <- sample(letters, 15)
#'
#' tSNE.embeddings <- matrix(data = rnorm(100),
#'                           nrow = 25,
#'                           ncol = 15)
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(colData = list(Letters = metaData))
#'
#' sce@@reducedDims@@listData$tSNE <- tSNE.embeddings
#'
#' plotData <- CellTagViz:::makePlotData(sce, "tSNE", "Letters")
#'

makePlotData <- function(sce, redMethod, metaVar, cells = FALSE){

  embeddings <- getEmbeddings(sce = sce, redMethod = redMethod, cells = cells)

  embeddings <- methods::as(embeddings[,1:2], "DataFrame")

  metaData <- getMetaData(sce = sce, varName = metaVar, cells = cells)

  plotData <- merge(embeddings, metaData, by = "row.names")

  plotData <- as.data.frame(plotData)

  return(plotData)

}
