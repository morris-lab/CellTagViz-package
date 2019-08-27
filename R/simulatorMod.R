
# This function defines the UI elements for a module.
# For the calculator only a few elements are needed.
# When creating UI shiny modules it is important to place the following command at the beginning of the function:
# ns <- NS(id)
# Including this command enables this UI module to be reused multiple times in the same shiny app.
# Each use of the UI module should be passed a unique ID.
# It this unique ID which enables us to differentiate the inputs from the multiple instances of a UI module.
# The id will be prefixed to the inputID for each input of the module.
# Such as calc-libComp/calc-moi or test-libComp/test-moi


calculatorUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::uiOutput(outputId = "calc-sidebar")
}


# This function defines the server module for this calculator.
# Server modules contain any functions necessary to perform the modules task.
# This can be data manipulation, calculations, or plotting.

calculatorServer <- function(input, output, session, userInputs) {

  # Perform simulation in python after the user clicks the Calculate button.
  # output$inputOut <- renderPrint({
  #
  #   reactiveValuesToList(userInputs)
  #
  # })

  # output$calcTabVal <- renderPrint({reactiveValuesToList(userInputs)})

  output$plot1 <- shiny::renderPlot({
    dupsMean <- mean(simResults())

    plotTitle <- paste0("Mean Fraction Duplicated: ", dupsMean)

    ggplot2::ggplot() + ggplot2::aes(simResults()) + ggplot2::geom_histogram() + ggplot2::xlab("Fraction Duplicates") + ggplot2::ylab("Count (Simulations)") + ggplot2::ggtitle(plotTitle) + ggplot2::theme_classic()
  })


  output$plot3 <- shiny::renderPlot({
    plotTitle <- paste0("Mean Clone Size Above Threshold: ", simResults()[[1]])

    ggplot2::ggplot() + ggplot2::aes(simResults()[[2]]) + ggplot2::geom_histogram() + ggplot2::ggtitle(plotTitle) + ggplot2::theme_classic() + ggplot2::xlab("Clone Size (cells)") + ggplot2::ylab("Count (Simulations)")
  })

  output$plot4 <- shiny::renderPlot({
    plotTitle <- paste0("Mean size of all clones: ", simResults()[[1]], "; Mean size of clones above threshold: ", simResults()[[2]])

    ggplot2::ggplot() + ggplot2::aes(simResults()[[3]]) + ggplot2::geom_histogram() + ggplot2::ggtitle(plotTitle) + ggplot2::theme_classic() + ggplot2::xlab("Clone Size (cells)") + ggplot2::ylab("Count (Simulations)")
  })

  output$troubleshootingtext <- shiny::renderUI({
    shiny::HTML(ifelse(userInputs$issues == 1,
      shiny::includeMarkdown("./issues/issue1.Rmd"),
      ifelse(userInputs$issues == 2,
        shiny::includeMarkdown("./issues/issue2.Rmd"),
        ifelse(userInputs$issues == 3,
          shiny::includeMarkdown("./issues/issue3.Rmd"),
          ifelse(userInputs$issues == 4,
            shiny::includeMarkdown("./issues/issue4.Rmd"),
            ifelse(userInputs$issues == 5,
              shiny::includeMarkdown("./issues/issue5.Rmd"),
              ifelse(userInputs$issues == 6,
                shiny::includeMarkdown("./issues/issue6.Rmd"),
                ifelse(userInputs$issues == 7,
                  shiny::includeMarkdown("./issues/issue7.Rmd"),
                  shiny::includeMarkdown("./issues/issue8.Rmd")
                )
              )
            )
          )
        )
      )
    ))
  })

  # output$paperlink <- renderUI({shiny::h4(shiny::tagList("For more information, please see our Nature Protocols paper ", a("here.", href="https://www.google.com/")))})

  output$plot1header <- shiny::renderUI({
    shiny::h3("Calculate Theoretical Duplicates")
  })

  output$plot3header <- shiny::renderUI({
    shiny::h3("Calculate Mean Clone Size in the Entire Cell Population Before Sequencing")
  })

  output$plot4header <- shiny::renderUI({
    shiny::h3("Calculate Sequenced Clone Size Distribution After a Single Experiment")
  })


  simResults <- shiny::eventReactive(
    userInputs$go, {
      switch(
        userInputs$calcTab,

        "one" = runDuplicatesExperiment(
          n_experiments = as.character(userInputs$numExperiments),
          N = as.character(userInputs$popSize),
          L = as.character(userInputs$libComp),
          MOI = as.character(userInputs$moi)
        ),

        "three" = meanSequencedCloneSize(
          n_simulations = userInputs$numExperiments,
          N = userInputs$popSize,
          L = userInputs$libComp,
          MOI = userInputs$moi,
          division_rate = userInputs$division_rate,
          passage_rate = userInputs$passage_rate,
          passage_fraction = userInputs$passage_fraction,
          sequence_time = userInputs$sequence_time,
          n_cells_sequenced = userInputs$n_cells_sequenced,
          threshold_clone_size = userInputs$threshold_clone_size
        ),

        "four" = cloneSizeSingleExperiment(
          N = userInputs$popSize,
          L = userInputs$libComp,
          MOI = userInputs$moi,
          division_rate = userInputs$division_rate,
          passage_rate = userInputs$passage_rate,
          passage_fraction = userInputs$passage_fraction,
          sequence_time = userInputs$sequence_time,
          n_cells_sequenced = userInputs$n_cells_sequenced,
          threshold_clone_size = userInputs$threshold_clone_size
        )
      )
    }
  )

  output$sidebar <- shiny::renderUI({
    switch(
      userInputs$calcTab,

      "one" = shiny::tagList(
        shiny::numericInput("libComp", "Library complexity", 200),

        shiny::numericInput("popSize", "Starting pop", 5000),

        shiny::numericInput("moi", "MOI", 4),

        shiny::numericInput("numExperiments", "Number of sims", 100),

        shiny::actionButton("go", "Calculate")
      ),

      "three" = shiny::tagList(
        shiny::numericInput("libComp", "Library complexity", 200),

        shiny::numericInput("popSize", "Starting pop", 5000),

        shiny::numericInput("moi", "MOI", 4),

        shiny::numericInput("numExperiments", "Number of sims", 100),

        shiny::numericInput("division_rate", "Division rate", 12),

        shiny::numericInput("passage_rate", "Passage rate", 30),

        shiny::sliderInput("passage_fraction", "Passage fraction", min = 0.01, max = 1.0, value = 0.8, step = 0.01),

        shiny::numericInput("sequence_time", "Sequence time", 168),

        shiny::numericInput("n_cells_sequenced", "Number of Cells Sequenced", 10000),

        shiny::numericInput("threshold_clone_size", "Threshold for Clone Size", 3),

        shiny::actionButton("go", "Calculate")
      ),

      "four" = shiny::tagList(
        shiny::numericInput("libComp", "Library complexity", 200),

        shiny::numericInput("popSize", "Starting pop", 5000),

        shiny::numericInput("moi", "MOI", 4),

        shiny::numericInput("division_rate", "Division rate", 12),

        shiny::numericInput("passage_rate", "Passage rate", 30),

        shiny::sliderInput("passage_fraction", "Passage fraction", min = 0.01, max = 1.0, value = 0.8, step = 0.01),

        shiny::numericInput("sequence_time", "Sequence time", 168),

        shiny::numericInput("n_cells_sequenced", "Number of Cells Sequenced", 10000),

        shiny::numericInput("threshold_clone_size", "Threshold for Clone Size", 3),

        shiny::actionButton("go", "Calculate")
      ),

      "five" = shiny::selectInput("issues", shiny::h4("Select one"),
        choices = list(
          "Where do I start?" = 1,
          "How does the CellTag Simulator work?" = 8,
          "My library complexity is low" = 2,
          "My cells were transduced with CellTags at a low efficiency" = 3,
          "A low percentage of profiled cells pass the CellTag expression threshold for cell tracking" = 4,
          "I have too many clonally-related cells" = 5,
          "I have too few clones, or small clones" = 6,
          "I found lineage collisions" = 7
        ),
        selected = 1
      )
    )
  })
}


# Defining this panel allows us to easily add a new page to our shiny app by calling this function.

calculatorPanel <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    title = "Calculator",
    # value = "calc",
    shiny::sidebarPanel(
      calculatorUI(id = id)
    ),

    shiny::mainPanel(
      shiny::tabsetPanel(
        type = "pills",
        id = "calcTab",
        selected = "one",

        shiny::tabPanel(
          title = "Duplicates",
          value = "one",
          shiny::uiOutput("calc-plot1header"),
          shiny::hr(),
          # verbatimTextOutput("calc-calcTabVal"),
          shinycssloaders::withSpinner(type = 7, shiny::plotOutput("calc-plot1"))
        ),

        shiny::tabPanel(
          title = "Clone Size Pre-Sequencing",
          value = "three",
          shiny::uiOutput("calc-plot3header"),
          shiny::hr(),
          shinycssloaders::withSpinner(type = 7, shiny::plotOutput("calc-plot3"))
        ),

        shiny::tabPanel(
          title = "Clone Size Distribution",
          value = "four",
          shiny::uiOutput("calc-plot4header"),
          shiny::hr(),
          shinycssloaders::withSpinner(type = 7, shiny::plotOutput("calc-plot4"))
        ),

        shiny::tabPanel(
          title = "Troubleshooting",
          value = "five",
          shiny::uiOutput("calc-troubleshootingtext"),
          shiny::hr()
          # ,imageOutput(ns("image")),
          # uiOutput("calc-paperlink")
        )
      )
    )
  )
}

# The UI and Server below create a stand alone app for our calculator.
# This is mainly for development and testing purposes.
# Ensuring the modularized app works correctly as a standalone tool will help with implementation with other shiny apps.
# If the app is working correctly by itself it should work when incorporated with another app.
# To include this app with another we just need to include a few commands to our existing app to implement our new calculator.

appUI <- shiny::navbarPage(
  title = "CellTag Simulator",

  theme = shinythemes::shinytheme(theme = "united"),

  # Include this command in the UI function of the existing app.
  calculatorPanel("calc")
)

appServer <- function(input, output, session) {

  # Include this command in the server function of the existing app.
  # userInputs may need to be changed to ensure the correct inputs are being passed.
  # This will mainly depend on the app being added too.
  shiny::callModule(
    module = calculatorServer,
    id = "calc",
    userInputs = input
  )
}

# Finally we need to run a few more commands in the global environment of the app.
# If you have a two file shiny app. These commands can be included in a new file named global.R
# This file should be in the same directory as the server.R and ui.R files.
# The commands in the global.R file are run to initialize the global environment for the shiny app.
# If our test app above were turned into a two file app. The global file would look like this.
# I will mark the commands we will probably need to place in the global.R file
# of the celltag.org shiny app.


# library(shiny)
# library(shinycssloaders) //This
# library(ggplot2)
# library(reticulate) //This
#
# source("modularCalcSimp.R) //This
#
# reticulate::source_python("duplicates.py") //This
#
