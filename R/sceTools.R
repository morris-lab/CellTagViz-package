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
#' \dontrun{
#'
#' sce <- SingleCellExperiment::SingleCellExperiment()
#'
#' reducedDims(sce, "tSNE") <- matrix(
#'   data = rnorm(100),
#'   nrow = 10,
#'   ncol = 10
#' )
#'
#' cellEmbeddings <- getEmbeddings(sce, "tSNE")
#' }
#'
getEmbeddings <- function(sce, redMethod, cells = FALSE) {
  embeddings <- SingleCellExperiment::reducedDim(sce, redMethod)

  if (cells) {
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
#' \dontrun{
#'
#' metaData <- sample(letters, 15)
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(colData = list(
#'   Letters = metaData
#' ))
#'
#' letters <- getMetaData(sce, "Letters")
#' }
#'
getMetaData <- function(sce, varName, cells = FALSE) {
  metaData <- SingleCellExperiment::colData(sce)[varName]

  if (cells) {
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
#' the given meta data which can then be used for grouping and coloring.
#'
#' @param sce SingleCellExperiment Object
#'
#' @param redMethod String Name of the dimension reduction method to visualize.
#'
#' @param metaVar String Name of the column to return from meta/col Data.
#'
#' @param feature String Name of a feature to return from assay data.
#'
#' @param cells Character Vector of Cell Barcodes used to subset data.
#'
#' @return This function returns a data frame which can be used to plot and
#' color cells from an SCE object using ggplot2.
#'
#' @examples
#'
#' \dontrun{
#'
#' metaData <- sample(letters, 15)
#'
#' tSNE.embeddings <- matrix(
#'   data = rnorm(250),
#'   nrow = 25,
#'   ncol = 10
#' )
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(colData = list(
#'   Letters = metaData
#' ))
#'
#' reducedDim(sce, "tSNE") <- tSNE.embeddings
#'
#' plotData <- makePlotData(sce, "tSNE", "Letters")
#' }
#'

makePlotData <-
  function(sce,
    redMethod,
    metaVar,
    feature = FALSE,
    cells = FALSE) {
    if (!shiny::isTruthy(redMethod)) {
      plotData <- SingleCellExperiment::colData(sce)

      plotData <- as.data.frame(plotData)

      if (shiny::isTruthy(feature)) {
        plotData <-
          addFeatureExpr(
            plotData = plotData,
            feature = feature,
            redMethod = "seurat",
            inputSCE = sce
          )

        return(plotData)
      }

      return(plotData)
    }

    embeddings <-
      getEmbeddings(sce = sce,
        redMethod = redMethod,
        cells = cells)

    embeddings <- methods::as(embeddings[, 1:2], "DataFrame")

    colnames(embeddings) <- c("Dim.1", "Dim.2")

    metaData <-
      getMetaData(sce = sce,
        varName = metaVar,
        cells = cells)

    plotData <- merge(embeddings, metaData, by = "row.names")

    plotData <- as.data.frame(plotData)

    rownames(plotData) <- plotData$Row.names

    if (shiny::isTruthy(feature)) {
      plotData <-
        addFeatureExpr(
          plotData = plotData,
          feature = feature,
          redMethod = redMethod,
          inputSCE = sce
        )
    }

    return(plotData)
  }


#' Adds gene expression to plot data.
#'
#' \code{addFeatureExpr} is a function used to add the expression values
#' for a chosen gene to be added to the data frame which contains all of the
#' data needed to construct the plot.
#'
#' @param plotData Data Frame which contains the cell embeddings and meta data
#' for the current plot.
#'
#' @param feature String of the gene name chosen by the user. Expression values
#' for this gene will be added to the plotData data frame.
#'
#' @param redMethod String The current dimension reduction method being
#' visualized. This value is used to determine which gene expression data to use.
#'
#' @param inputSCE SingleCellExperiment object which contains data being used.
#'
#' @return The function returns the given data frame with a column added which
#' contains the expression values for the user chosen gene.
#'
#' @examples
#'
#' \dontrun{
#'
#'
#' plotData <- makePlotData(blah, blah, blah)
#'
#' plotData <- addFeatureExpr(foo, bar, foo)
#'
#' }
#'
addFeatureExpr <- function(plotData, feature, redMethod, inputSCE) {
  if (grepl(pattern = "Seurat", redMethod, ignore.case = TRUE)) {
    exprData <-
      SummarizedExperiment::assay(inputSCE, "ScaleData.Seurat")[feature, ]
  } else if (grepl(pattern = "Monocle", redMethod, ignore.case = TRUE)) {
    exprData <-
      SummarizedExperiment::assay(inputSCE, "Counts.Monocle")[feature, ]
  }

  exprData <- as.data.frame(exprData)

  plotData[[feature]] <- 0

  plotData[rownames(exprData), feature] <-
    exprData[rownames(exprData), 1]

  return(plotData)
}
