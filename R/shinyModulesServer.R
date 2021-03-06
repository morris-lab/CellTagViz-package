



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
#' @param inputSCE Single Cell Experiment Object which contains the data to plot.
#'
#' @return Returns plots based on user input values and current tab selection.
#'
#' @export
#'

createPlot <-
  function(input,
             output,
             session,
             plotOptions,
             inputSCE) {
    output$inputOut <- shiny::renderPrint({
      plotOpts <- shiny::reactiveValuesToList(plotOptions)

      utils::str(plotOpts)
    })

    output$testPlot <- shiny::renderPlot({
      plotOpts <- shiny::reactiveValuesToList(plotOptions)

      switch(
        plotOpts$plots,

        "tSNE" = plotTsne(
          sce = inputSCE,
          metaVar = plotOpts$`TSNE-META`,
          factor = plotOpts$`TSNE-isFactor`,
          contour = plotOpts$`TSNE-addContour`,
          feature = plotOpts$`TSNE-GENE`,
          redMethod = plotOpts$`TSNE-REDUCTION`,
          exprDat = plotOpts$`TSNE-EXPR`
        ),

        "Network" = plotNetwork(sce = datasets::cars),

        "Pseudotime" = plotPseudotime(
          sce = inputSCE,
          metaVar = plotOpts$`PSEUDO-META`,
          factor = plotOpts$`PSEUDO-isFactor`,
          contour = plotOpts$`PSEUDO-addContour`,
          feature = plotOpts$`PSEUDO-GENE`
        ),

        "Stacked Bar Charts" = plotStackedBar(
          sce = inputSCE,
          colorVar = plotOpts$`BARCHART-GROUP-META`,
          groupVar = plotOpts$`BARCHART-META`,
          horizontal = plotOpts$`BARCHART-plotFlip`
        ),

        "Scatter Plots" = plotScatter(
          sce = inputSCE,
          metaVar = plotOpts$`SCATTER-META`,
          feature = plotOpts$`SCATTER-GENE`,
          factor = plotOpts$`SCATTER-isFactor`,
          exprDat = plotOpts$`SCATTER-EXPR`
        ),

        "Meta Data" = plotMeta(sce = datasets::cars)
      )
    })
  }



easyPlot <- function(...) {
  p <- ggplot2::quickplot(x = datasets::cars$speed, y = datasets::cars$dist, ...)

  return(p)
}



#' Parses shiny user input and plots dimension reduction visualiztions.
#'
#' This function generates 2d visualizations of single-cell RNA data. The
#' function parses the user input and plots the visualization according to the
#' given user inputs.
#'
#' @param sce Yes
#'
#' @param feature Yes
#'
#' @param metaVar Yes
#'
#' @param clones Yes
#'
#' @param factor Yes
#'
#' @param contour Yes
#'
#' @param redMethod Yes
#'
#' @param exprDat Yes
#'
#' @return This function returns a ggplot2 object
#'


plotTsne <-
  function(sce,
             feature = FALSE,
             metaVar = FALSE,
             clones = FALSE,
             factor = FALSE,
             contour = FALSE,
             redMethod = FALSE,
             exprDat = FALSE) {
    plotData <-
      makePlotData(
        sce = sce,
        redMethod = redMethod,
        metaVar = metaVar,
        feature = feature,
        exprDat = exprDat
      )

    if (factor) {
      plotData[[metaVar]] <- as.factor(plotData[[metaVar]])
    }

    if (shiny::isTruthy(feature)) {
      b <-
        ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(
          x = ~Dim.1,
          y = ~Dim.2,
          color = as.name(feature)
        )) + viridis::scale_color_viridis()
    } else {
      b <-
        ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(
          x = ~Dim.1,
          y = ~Dim.2,
          color = as.name(metaVar)
        ))
    }

    br <-
      b + ggplot2::labs(title = "tSNE") + ggplot2::theme_classic()

    if (contour) {
      datCols <-
        grep(
          pattern = "Dim",
          x = colnames(plotData),
          ignore.case = TRUE,
          value = TRUE
        )

      bre <-
        br + plotContour(
          plotData = plotData,
          metaVar = metaVar,
          dataCols = datCols
        )

      return(bre)
    }

    return(br)
  }

# ==============================================================================
# ==============================================================================

plotNetwork <-
  function(sce,
             feature = FALSE,
             metaVar = FALSE,
             clones = FALSE,
             factor = FALSE,
             contour = FALSE) {
    b <-
      ggplot2::ggplot(data = sce) + ggplot2::geom_point(ggplot2::aes(x = sce[[1]], y = sce[[2]]))

    br <-
      b + ggplot2::labs(title = "Network") + ggplot2::theme_classic()

    return(br)
  }

# ==============================================================================
# ==============================================================================

plotPseudotime <-
  function(sce,
             feature = FALSE,
             metaVar = FALSE,
             clones = FALSE,
             factor = FALSE,
             contour = FALSE) {
    plotData <-
      makePlotData(
        sce = sce,
        redMethod = "Pseudotime.Monocle",
        metaVar = metaVar,
        feature = feature
      )

    if (factor) {
      plotData[[metaVar]] <- as.factor(plotData[[metaVar]])
    }

    if (shiny::isTruthy(feature)) {
      b <-
        ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(
          x = ~Component.1,
          y = ~Component.2,
          color = as.name(feature)
        )) + viridis::scale_color_viridis()
    } else {
      b <-
        ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(
          x = ~Component.1,
          y = ~Component.2,
          color = as.name(metaVar)
        ))
    }

    br <-
      b + ggplot2::labs(title = "Pseudotime") + ggplot2::theme_classic()

    if (contour) {
      x <-
        grep(
          pattern = "1",
          x = colnames(plotData),
          ignore.case = TRUE,
          value = TRUE
        )

      y <-
        grep(
          pattern = "2",
          x = colnames(plotData),
          ignore.case = TRUE,
          value = TRUE
        )

      datCols <- c(x, y)

      bre <-
        br + plotContour(
          plotData = plotData,
          metaVar = metaVar,
          dataCols = datCols
        )

      return(bre)
    }

    return(br)
  }

# ==============================================================================
# ==============================================================================

plotStackedBar <-
  function(sce,
             groupVar = FALSE,
             colorVar = FALSE,
             clones = FALSE,
             horizontal = FALSE) {
    plotData <- makePlotData(sce = sce, redMethod = FALSE)

    b <-
      ggplot2::ggplot(data = plotData) + ggplot2::geom_bar(ggplot2::aes_(x = as.name(groupVar), fill = as.name(colorVar)),
        position = "fill",
        na.rm = TRUE
      ) + ggplot2::scale_y_continuous(labels = scales::percent)

    if (shiny::isTruthy(horizontal)) {
      br <-
        b + ggplot2::labs(title = "Stacked Bar Charts") + ggplot2::theme_classic() + ggplot2::coord_flip()
    } else {
      (br <-
        b + ggplot2::labs(title = "Stacked Bar Charts") + ggplot2::theme_classic())
    }

    return(br)
  }

# ==============================================================================
# ==============================================================================

plotScatter <-
  function(sce,
             feature = "Apoa1",
             metaVar = FALSE,
             clones = FALSE,
             factor = FALSE,
             contour = FALSE,
             exprDat = FALSE) {
    plotData <-
      makePlotData(
        sce = sce,
        feature = feature,
        metaVar = metaVar,
        redMethod = FALSE,
        exprDat = exprDat
      )

    if (shiny::isTruthy(factor)) {
      plotData[[metaVar]] <- as.factor(plotData[[metaVar]])

      b <-
        ggplot2::ggplot(data = plotData) + ggplot2::geom_violin(ggplot2::aes_(
          x = as.name(metaVar),
          y = as.name(feature),
          fill = as.name(metaVar)
        ),
        na.rm = TRUE
        )
    } else {
      (
        b <-
          ggplot2::ggplot(data = plotData) + ggplot2::geom_point(ggplot2::aes_(
            x = as.name(metaVar), y = as.name(feature)
          ), na.rm = TRUE)
      )
    }

    br <-
      b + ggplot2::labs(title = "Scatter Plots") + ggplot2::theme_classic()

    return(br)
  }

# ==============================================================================
# ==============================================================================

plotMeta <-
  function(sce,
             feature = FALSE,
             metaVar = FALSE,
             clones = FALSE,
             factor = FALSE,
             contour = FALSE) {
    b <-
      ggplot2::ggplot(data = sce) + ggplot2::geom_point(ggplot2::aes(x = sce[[1]], y = sce[[2]]))

    br <-
      b + ggplot2::labs(title = "Meta Data") + ggplot2::theme_classic()

    return(br)
  }

# ==============================================================================
# ==============================================================================

plotContour <- function(plotData, metaVar, dataCols) {
  dimX <- dataCols[[1]]

  dimY <- dataCols[[2]]

  contourLayer <-
    ggplot2::stat_density2d(
      data = plotData,
      ggplot2::aes_(
        x = as.name(dimX),
        y = as.name(dimY),
        color = as.name(metaVar)
      )
    )
}
