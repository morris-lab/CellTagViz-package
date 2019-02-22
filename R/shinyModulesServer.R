

################################################################################
# This file contains functions which define the shiny server modules.
################################################################################

#
# plotOpts <- function(input, output, session, brent){
#
#   output$inputOut <- renderPrint({
#
#     str(sapply(names(brent), function(id){brent[[id]]}, simplify = FALSE))
#
#   })
#
#   output$testPlot <- renderPlot({
#
#     ggplot2::quickplot(x = cars$speed, y = cars$dist)
#
#   })
#
# }
#






createPlot <- function(input, output, session, plotOptions){



  output$inputOut <- renderPrint({

    plotOpts <- reactiveValuesToList(plotOptions)

    str(plotOpts)

  })

  output$testPlot <- renderPlot({

    plotOpts <- reactiveValuesToList(plotOptions)

    switch(plotOpts$plots,

      "tSNE" = plotTsne(sce = sce, metaVar = plotOpts$`TSNE-META`, factor = plotOpts$`TSNE-isFactor`, contour = plotOpts$`TSNE-addContour`, feature = plotOpts$`TSNE-GENE`),

      "Network" = plotNetwork(sce = cars),

      "Pseudotime" = plotPseudotime(sce = sce, metaVar = plotOpts$`PSEUDO-META`, factor = plotOpts$`PSEUDO-isFactor`, contour = plotOpts$`PSEUDO-addContour`, feature = plotOpts$`PSEUDO-GENE`),

      "Stacked Bar Charts" = plotStackedBar(sce = sce, colorVar = plotOpts$`BARCHART-GROUP-META`, groupVar = plotOpts$`BARCHART-META`, horizontal = plotOpts$`BARCHART-plotFlip`),

      "Scatter Plots" = plotScatter(sce = sce, metaVar = plotOpts$`SCATTER-META`, feature = plotOpts$`SCATTER-GENE`, factor = plotOpts$`SCATTER-isFactor`),

      "Meta Data" = plotMeta(sce = cars)

      )
  })
}



easyPlot <- function(...){

  p <- ggplot2::quickplot(x = cars$speed, y = cars$dist, ...)

  return(p)

}


plotTsne <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE){

  plotData <- CellTagViz:::makePlotData(sce = sce, redMethod = "tsne.Seurat", metaVar = metaVar, feature = feature)

  if(factor){

    plotData[[metaVar]] <- as.factor(plotData[[metaVar]])

  }

  if(isTruthy(feature)){

    b <- ggplot2::ggplot(data = plotData) + geom_point(aes_(x = ~tSNE_1, y = ~tSNE_2, color = as.name(feature))) + viridis::scale_color_viridis()

  } else{

    b <- ggplot2::ggplot(data = plotData) + geom_point(aes_(x = ~tSNE_1, y = ~tSNE_2, color = as.name(metaVar)))

  }

  br <- b + labs(title = "tSNE") + theme_classic()

  if(contour){

    datCols <- grep(pattern = "tsne", x = colnames(plotData), ignore.case = TRUE, value = TRUE)

    bre <- br + plotContour(plotData = plotData, metaVar = metaVar, dataCols = datCols)

    return(bre)

  }

  return(br)

}


plotNetwork <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE){

  b <- ggplot2::ggplot(data = sce) + geom_point(aes(x = sce[[1]], y = sce[[2]]))

  br <- b + labs(title = "Network") + theme_classic()

  return(br)

}


plotPseudotime <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE){

  plotData <- CellTagViz:::makePlotData(sce = sce, redMethod = "Pseudotime.Monocle", metaVar = metaVar, feature = feature)

  if(factor){

    plotData[[metaVar]] <- as.factor(plotData[[metaVar]])

  }

  if(isTruthy(feature)){

    b <- ggplot2::ggplot(data = plotData) + geom_point(aes_(x = ~Component.1, y = ~Component.2, color = as.name(feature))) + viridis::scale_color_viridis()

  } else{

    b <- ggplot2::ggplot(data = plotData) + geom_point(aes_(x = ~Component.1, y = ~Component.2, color = as.name(metaVar)))

  }

  br <- b + labs(title = "Pseudotime") + theme_classic()

  if(contour){

    x <- grep(pattern = "1", x = colnames(plotData), ignore.case = TRUE, value = TRUE)

    y <- grep(pattern = "2", x = colnames(plotData), ignore.case = TRUE, value = TRUE)

    datCols <- c(x, y)

    bre <- br + plotContour(plotData = plotData, metaVar = metaVar, dataCols = datCols)

    return(bre)

  }

  return(br)

}



plotStackedBar <- function(sce, groupVar = FALSE, colorVar = FALSE, clones = FALSE, horizontal = FALSE){

  plotData <- CellTagViz:::makePlotData(sce = sce, redMethod = FALSE)

  b <- ggplot2::ggplot(data = plotData) + geom_bar(aes_(x = as.name(groupVar), fill = as.name(colorVar)), position = "fill", na.rm = TRUE) + scale_y_continuous(labels = scales::percent)

  if(isTruthy(horizontal)){

    br <- b + labs(title = "Stacked Bar Charts") + theme_classic() + coord_flip()

  } else(

    br <- b + labs(title = "Stacked Bar Charts") + theme_classic()

  )

  return(br)

}



plotScatter <- function(sce, feature = "Apoa1", metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE){

  plotData <- CellTagViz:::makePlotData(sce = sce, feature = feature, metaVar = metaVar, redMethod = FALSE)

  if(shiny::isTruthy(factor)){

    plotData[[metaVar]] <- as.factor(plotData[[metaVar]])

    b <- ggplot2::ggplot(data = plotData) + geom_violin(aes_(x = as.name(metaVar), y = as.name(feature), fill = as.name(metaVar)), na.rm = TRUE)

  } else(

    b <- ggplot2::ggplot(data = plotData) + geom_point(aes_(x = as.name(metaVar), y = as.name(feature)), na.rm = TRUE)
  )

  br <- b + labs(title = "Scatter Plots") + theme_classic()

  return(br)

}



plotMeta <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE){

  b <- ggplot2::ggplot(data = sce) + geom_point(aes(x = sce[[1]], y = sce[[2]]))

  br <- b + labs(title = "Meta Data") + theme_classic()

  return(br)

}


plotContour <- function(plotData, metaVar, dataCols){

  dimX <- dataCols[[1]]

  dimY <- dataCols[[2]]

  contourLayer <- stat_density2d(data = plotData, aes_(x = as.name(dimX), y = as.name(dimY), color = as.name(metaVar)))

}



addExprToPlot <- function(plotData, feature, exprAssay){

  plotData[[feature]] <- assay(sce, exprAssay)[feature, ]

  return(plotData)

}





