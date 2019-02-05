

getEmbeddings <- function(sce, redMethod, cells = NULL){
  
  embeddings <- as(reducedDim(sce, redMethod), rownames = NA)
  
  if(cells){
    
    embeddings <- embeddings[cells, ]
  }
  
  return(embeddings)
}

getMetaData <- function(sce, varName, cells = NULL){
  
  metaData <- colData(sce)[varName]
  
  metaData <- met
  
  if(cells){
    
    metaData <- metaData[cells, ]
  }
  
}
