

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

