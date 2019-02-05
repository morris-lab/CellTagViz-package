
#'Access the cell embeddings of SCE object.
#'
#'\code{getEmbeddings} returns the coordinates of cells after dimension
#'reduction.
#'
#'This is a function used to subset SCE objects. The output is used to generate
#'a data frame which is compatible with ggplot. The cell embeddings are then
#'combined with metadata to generate visualizations with ggplot2.
#'
#'@param sce SingleCellExperiment Object
#'
#'@param redMethod String Name of dimension reduction method to visualize.
#'Must be one of methods listed in \cod{names(sce@@reducedDim)}
#'
#'@param cells Character List of cells used to subset data. (Optional)
#'
#'@return The function returns the cell embeddings stored in the given slot of
#'the SingleCellExperiment object.
#'
#'
#'@examples
#'
#'cellEmbeddings <- getEmbeddings(sce, "cca.Seurat")
#'
#'
#'

getEmbeddings <- function(sce, redMethod, cells = FALSE){

  embeddings <- reducedDim(sce, redMethod)

  if(cells){

    embeddings <- embeddings[cells, ]
  }

  return(embeddings)
}

getMetaData <- function(sce, varName, cells = FALSE){

  metaData <- colData(sce)[varName]

  if(cells){

    metaData <- metaData[cells, ]
  }

  return(metaData)

}

makePlotData <- function(sce, redMethod, metaVar, cells = FALSE){

  embeddings <- getEmbeddings(sce = sce, redMethod = redMethod, cells = cells)

  embeddings <- as(embeddings[,1:2], "DataFrame")

  metaData <- getMetaData(sce = sce, varName = metaVar, cells = cells)

  plotData <- merge(embeddings, metaData, by = "row.names")

  plotData <- as.data.frame(plotData)

  return(plotData)

}
