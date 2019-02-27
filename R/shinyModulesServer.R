

################################################################################
# This file contains functions which define the shiny server modules. Many of
# these modules are functions which return ggplot2 plots.
################################################################################

#' This is a shiny server function for generating CellTagViz plots.
#'
#' \code{createPlot} is a function called in the server function of the
#' CellTagViz shiny app. The function assigns the reactive user inputs as a
#' variable. These user inputs are then parsed and passed as arguments to the
#' appropriate plotting functions. This function takes in the user input, parses
#' the user input, and passes the input to the correct function depending on
#' which tab is currently being displayed.
#'
#' @param input Default argument required for shiny modules.
#'
#' @param output Default argument required for shiny modules.
#'
#' @param session Default argument required for shiny modules.
#'
#' @param plotOptions List which contains the user input values.
#'
#' @return Returns plots based on user input values and current tab selection.
#'
#'

createPlot <- function(input, output, session, plotOptions) {
  output$inputOut <- shiny::renderPrint({
    plotOpts <- shiny::reactiveValuesToList(plotOptions)

    utils::str(plotOpts)
  })

  output$testPlot <- shiny::renderPlot({
    plotOpts <- shiny::reactiveValuesToList(plotOptions)

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



easyPlot <- function(...) {
  p <- ggplot2::quickplot(x = cars$speed, y = cars$dist, ...)

  return(p)
}


plotTsne <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE) {
  plotData <- makePlotData(sce = sce, redMethod = "tsne.Seurat", metaVar = metaVar, feature = feature)

  if (factor) {
    plotData[[metaVar]] <- as.factor(plotData[[metaVar]])
  }

  if (shiny::isTruthy(feature)) {
    b <- ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(x = ~tSNE_1, y = ~tSNE_2, color = as.name(feature))) + viridis::scale_color_viridis()
  } else {
    b <- ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(x = ~tSNE_1, y = ~tSNE_2, color = as.name(metaVar)))
  }

  br <- b + ggplot2::labs(title = "tSNE") + ggplot2::theme_classic()

  if (contour) {
    datCols <- grep(pattern = "tsne", x = colnames(plotData), ignore.case = TRUE, value = TRUE)

    bre <- br + plotContour(plotData = plotData, metaVar = metaVar, dataCols = datCols)

    return(bre)
  }

  return(br)
}


plotNetwork <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE) {
  b <- ggplot2::ggplot(data = sce) + ggplot2::geom_point(ggplot2::aes(x = sce[[1]], y = sce[[2]]))

  br <- b + ggplot2::labs(title = "Network") + ggplot2::theme_classic()

  return(br)
}


plotPseudotime <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE) {
  plotData <- makePlotData(sce = sce, redMethod = "Pseudotime.Monocle", metaVar = metaVar, feature = feature)

  if (factor) {
    plotData[[metaVar]] <- as.factor(plotData[[metaVar]])
  }

  if (shiny::isTruthy(feature)) {
    b <- ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(x = ~Component.1, y = ~Component.2, color = as.name(feature))) + viridis::scale_color_viridis()
  } else {
    b <- ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(x = ~Component.1, y = ~Component.2, color = as.name(metaVar)))
  }

  br <- b + ggplot2::labs(title = "Pseudotime") + ggplot2::theme_classic()

  if (contour) {
    x <- grep(pattern = "1", x = colnames(plotData), ignore.case = TRUE, value = TRUE)

    y <- grep(pattern = "2", x = colnames(plotData), ignore.case = TRUE, value = TRUE)

    datCols <- c(x, y)

    bre <- br + plotContour(plotData = plotData, metaVar = metaVar, dataCols = datCols)

    return(bre)
  }

  return(br)
}



plotStackedBar <- function(sce, groupVar = FALSE, colorVar = FALSE, clones = FALSE, horizontal = FALSE) {
  plotData <- makePlotData(sce = sce, redMethod = FALSE)

  b <- ggplot2::ggplot(data = plotData) + ggplot2::geom_bar(ggplot2::aes_(x = as.name(groupVar), fill = as.name(colorVar)), position = "fill", na.rm = TRUE) + ggplot2::scale_y_continuous(labels = scales::percent)

  if (shiny::isTruthy(horizontal)) {
    br <- b + ggplot2::labs(title = "Stacked Bar Charts") + ggplot2::theme_classic() + ggplot2::coord_flip()
  } else {
    (

      br <- b + ggplot2::labs(title = "Stacked Bar Charts") + ggplot2::theme_classic()

    )
  }

  return(br)
}



plotScatter <- function(sce, feature = "Apoa1", metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE) {
  plotData <- makePlotData(sce = sce, feature = feature, metaVar = metaVar, redMethod = FALSE)

  if (shiny::isTruthy(factor)) {
    plotData[[metaVar]] <- as.factor(plotData[[metaVar]])

    b <- ggplot2::ggplot(data = plotData) + ggplot2::geom_violin(ggplot2::aes_(x = as.name(metaVar), y = as.name(feature), fill = as.name(metaVar)), na.rm = TRUE)
  } else {
    (

      b <- ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(x = as.name(metaVar), y = as.name(feature)), na.rm = TRUE)
    )
  }

  br <- b + ggplot2::labs(title = "Scatter Plots") + ggplot2::theme_classic()

  return(br)
}



plotMeta <- function(sce, feature = FALSE, metaVar = FALSE, clones = FALSE, factor = FALSE, contour = FALSE) {
  b <- ggplot2::ggplot(data = sce) + ggplot2::geom_point(ggplot2::aes(x = sce[[1]], y = sce[[2]]))

  br <- b + ggplot2::labs(title = "Meta Data") + ggplot2::theme_classic()

  return(br)
}


plotContour <- function(plotData, metaVar, dataCols) {
  dimX <- dataCols[[1]]

  dimY <- dataCols[[2]]

  contourLayer <- ggplot2::stat_density2d(data = plotData, ggplot2::aes_(x = as.name(dimX), y = as.name(dimY), color = as.name(metaVar)))
}



addExprToPlot <- function(plotData, feature, exprAssay) {
  plotData[[feature]] <- SummarizedExperiment::assay(sce, exprAssay)[feature, ]

  return(plotData)
}
