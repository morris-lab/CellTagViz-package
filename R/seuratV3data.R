





seuratV3 <- function(seuratObj){

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

           emptyMat <- Matrix::Matrix(data = 0, nrow = nrow(rowDat), ncol = nrow(colDat), sparse = TRUE)

           colnames(emptyMat) <- rownames(colDat)

           rownames(emptyMat) <- rowDat[[1]]

           emptyMat[rownames(temp), colnames(temp)] <- temp[rownames(temp), colnames(temp)]

           exprList[[id]] <- emptyMat

         } else(

          exprList[[id]] <- temp

        )
      }
    }
  }

  redMethods <- names(seuratObj@reductions)

  names(redMethods) <- paste0(redMethods, ".Seurat")

  cellEmbeddings <- lapply(redMethods, function(id){

    Seurat::Embeddings(seuratObj, reduction = id)

  })

  cellEmbeddings <- as(cellEmbeddings, "SimpleList")

  sce <- SingleCellExperiment::SingleCellExperiment(colData = colDat, rowData = rowDat, assays = exprList, reducedDims = cellEmbeddings)

}

