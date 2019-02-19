

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

      "tSNE" = easyPlot(main = "tSNE"),

      "Network" = easyPlot(main = "Network"),

      "Pseudotime" = easyPlot(main = "Pseudotime"),

      "Stacked Bar Charts" = easyPlot(main = "Stacked Bar Charts"),

      "Scatter Plots" = easyPlot(main = "Scatter Plots"),

      "Meta Data" = easyPlot(main = "Meta Data")

      )
  })
}



easyPlot <- function(...){

  p <- ggplot2::quickplot(x = cars$speed, y = cars$dist, ...)

  return(p)

}




