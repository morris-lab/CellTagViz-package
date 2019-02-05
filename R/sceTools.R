

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
