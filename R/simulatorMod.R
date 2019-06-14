
#This function defines the UI elements for a module.
#For the calculator only a few elements are needed.
#When creating UI shiny modules it is important to place the following command at the beginning of the function:
#ns <- NS(id)
#Including this command enables this UI module to be reused multiple times in the same shiny app.
#Each use of the UI module should be passed a unique ID.
#It this unique ID which enables us to differentiate the inputs from the multiple instances of a UI module.
#The id will be prefixed to the inputID for each input of the module.
# Such as calc-libComp/calc-moi or test-libComp/test-moi


calculatorUI <- function(id){

	ns <- shiny::NS(id)
	
	uiOutput(outputId = "calc-sidebar")
	
}


#This function defines the server module for this calculator. 
#Server modules contain any functions necessary to perform the modules task.
#This can be data manipulation, calculations, or plotting.

calculatorServer <- function(input, output, session, userInputs){
  
  #Perform simulation in python after the user clicks the Calculate button.
	 # output$inputOut <- renderPrint({
	 #  
	 #   reactiveValuesToList(userInputs)
	 #  
	 # })
	
	#output$calcTabVal <- renderPrint({reactiveValuesToList(userInputs)})
	
	output$plot1 <- renderPlot({
	  
	  dupsMean <- mean(simResults())
	  
	  plotTitle <- paste0("Mean Fraction Duplicated: ", dupsMean)
	  
	  ggplot() + aes(simResults()) + geom_histogram() + xlab("Fraction Duplicates") + ylab("Count (Simulations)") + ggtitle(plotTitle) + theme_classic()
	  
	})
	
	
	output$plot3 <- renderPlot({
	  
	  plotTitle <- paste0("Mean Clone Size Above Threshold: ", simResults()[[1]])
	  
	  ggplot() + aes(simResults()[[2]]) + geom_histogram() + ggtitle(plotTitle) + theme_classic() + xlab("Clone Size (cells)") + ylab("Count (Simulations)")
	  
	})
	
	output$plot4 <- renderPlot({
	  
	  plotTitle <- paste0("Mean size of all clones: ", simResults()[[1]], "; Mean size of clones above threshold: ", simResults()[[2]])
	  
	  ggplot() + aes(simResults()[[3]]) + geom_histogram() + ggtitle(plotTitle) + theme_classic() + xlab("Clone Size (cells)") + ylab("Count (Simulations)")
	  
  })
	
	output$troubleshootingtext <- renderUI({HTML(ifelse(userInputs$issues == 1,
	                                                    includeMarkdown("./issues/issue1.Rmd"),
	                                                    ifelse(userInputs$issues == 2,
	                                                           includeMarkdown("./issues/issue2.Rmd"),
	                                                           ifelse(userInputs$issues == 3,
	                                                                  includeMarkdown("./issues/issue3.Rmd"),
	                                                                  ifelse(userInputs$issues == 4,
	                                                                         includeMarkdown("./issues/issue4.Rmd"),
	                                                                         ifelse(userInputs$issues == 5,
	                                                                                includeMarkdown("./issues/issue5.Rmd"),
	                                                                                ifelse(userInputs$issues == 6,
	                                                                                       includeMarkdown("./issues/issue6.Rmd"),
	                                                                                       ifelse(userInputs$issues == 7,
	                                                                                              includeMarkdown("./issues/issue7.Rmd"),
	                                                                                              includeMarkdown("./issues/issue8.Rmd")
	                                                                                              ))))))))})
	
	#output$paperlink <- renderUI({h4(tagList("For more information, please see our Nature Protocols paper ", a("here.", href="https://www.google.com/")))})
	
	output$plot1header <- renderUI({h3("Calculate Theoretical Duplicates")})
	
	output$plot3header <- renderUI({h3("Calculate Mean Clone Size in the Entire Cell Population Before Sequencing")})
	
	output$plot4header <- renderUI({h3("Calculate Sequenced Clone Size Distribution After a Single Experiment")})

	
	simResults <- eventReactive(
	  userInputs$go,
	  
	  {switch(
	      userInputs$calcTab,
	      
	      "one" = runDuplicatesExperiment(
	        n_experiments = as.character(userInputs$numExperiments),
	        N = as.character(userInputs$popSize),
	        L = as.character(userInputs$libComp), 
	        MOI = as.character(userInputs$moi)),
	      
	      "three" =  meanSequencedCloneSize(
	        n_simulations = userInputs$numExperiments,
	        N = userInputs$popSize,
	        L = userInputs$libComp,
	        MOI = userInputs$moi,
	        division_rate = userInputs$division_rate,
	        passage_rate = userInputs$passage_rate,
	        passage_fraction = userInputs$passage_fraction,
	        sequence_time = userInputs$sequence_time,
	        n_cells_sequenced = userInputs$n_cells_sequenced,
	        threshold_clone_size = userInputs$threshold_clone_size),
	      
	      "four" = cloneSizeSingleExperiment(
	        N = userInputs$popSize,
	        L = userInputs$libComp,
	        MOI = userInputs$moi,
	        division_rate = userInputs$division_rate,
	        passage_rate = userInputs$passage_rate,
	        passage_fraction = userInputs$passage_fraction,
	        sequence_time = userInputs$sequence_time,
	        n_cells_sequenced = userInputs$n_cells_sequenced,
	        threshold_clone_size = userInputs$threshold_clone_size)
	      
	      )
	    })
	
	output$sidebar <- renderUI({
	  
	  switch(
	    userInputs$calcTab,
	    
	    "one" = tagList(
	      
	      numericInput("libComp", "Library complexity", 200),
	      
	      numericInput("popSize", "Starting pop", 5000),
	      
	      numericInput("moi", "MOI", 4),
	      
	      numericInput("numExperiments", "Number of sims", 100),
	      
	      actionButton("go", "Calculate")
	    ),
	    
	    "three" = tagList(
	      
	      numericInput("libComp", "Library complexity", 200),
	      
	      numericInput("popSize", "Starting pop", 5000),
	      
	      numericInput("moi", "MOI", 4),
	      
	      numericInput("numExperiments", "Number of sims", 100),
	      
	      numericInput("division_rate", "Division rate", 12),
	      
	      numericInput("passage_rate", "Passage rate", 30),
	      
	      sliderInput("passage_fraction", "Passage fraction", min = 0.01, max = 1.0, value = 0.8, step = 0.01),
	      
	      numericInput("sequence_time", "Sequence time", 168),
	      
	      numericInput("n_cells_sequenced", "Number of Cells Sequenced", 10000),
	      
	      numericInput("threshold_clone_size", "Threshold for Clone Size", 3),
	      
	      actionButton("go", "Calculate")
	      
	    ),
	    
	    "four" = tagList(
	      
	      numericInput("libComp", "Library complexity", 200),
	      
	      numericInput("popSize", "Starting pop", 5000),
	      
	      numericInput("moi", "MOI", 4),
	      
	      numericInput("division_rate", "Division rate", 12),
	      
	      numericInput("passage_rate", "Passage rate", 30),
	      
	      sliderInput("passage_fraction", "Passage fraction", min = 0.01, max = 1.0, value = 0.8, step = 0.01),
	      
	      numericInput("sequence_time", "Sequence time", 168),
	      
	      numericInput("n_cells_sequenced", "Number of Cells Sequenced", 10000),
	      
	      numericInput("threshold_clone_size", "Threshold for Clone Size", 3),
	      
	      actionButton("go", "Calculate")
	      
	    ),
	    
	    "five" = selectInput("issues", h4("Select one"), choices = list("Where do I start?" = 1,
	                                                                    "How does the CellTag Simulator work?" = 8,
	                                                                    "My library complexity is low" = 2,
	                                                                    "My cells were transduced with CellTags at a low efficiency" = 3,
	                                                                    "A low percentage of profiled cells pass the CellTag expression threshold for cell tracking" = 4,
	                                                                    "I have too many clonally-related cells" = 5,
	                                                                    "I have too few clones, or small clones" = 6,
	                                                                    "I found lineage collisions" = 7
	                                                             ),
	                  selected = 1)
	    
	  )
	  
	})
	
	
}


#Defining this panel allows us to easily add a new page to our shiny app by calling this function. 

calculatorPanel <- function(id){
  
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    
    title = "Calculator",
    #value = "calc",
    shiny::sidebarPanel(
      
      calculatorUI(id = id)
      
    ),
    
    shiny::mainPanel(
      
      tabsetPanel(
        type = "pills",
        id = "calcTab",
        selected = "one",
        
        tabPanel(
          title = "Duplicates", 
          value = "one",
          uiOutput("calc-plot1header"),
          hr(),
          #verbatimTextOutput("calc-calcTabVal"),
          shinycssloaders::withSpinner(type = 7, plotOutput("calc-plot1"))
        ),
        
        tabPanel(
          title = "Clone Size Pre-Sequencing",
          value = "three",
          uiOutput("calc-plot3header"),
          hr(),
          shinycssloaders::withSpinner(type = 7, plotOutput("calc-plot3"))
        ),
        
        tabPanel(
          title = "Clone Size Distribution",
          value = "four",
          uiOutput("calc-plot4header"),
          hr(),
          shinycssloaders::withSpinner(type = 7, plotOutput("calc-plot4"))
        ),
        
        tabPanel(
          title = "Troubleshooting",
          value = "five",
          uiOutput("calc-troubleshootingtext"),
          hr()
          #,imageOutput(ns("image")),
          #uiOutput("calc-paperlink")
        )
        
      )
      
    )
  )
  
}

#The UI and Server below create a stand alone app for our calculator.
#This is mainly for development and testing purposes. 
#Ensuring the modularized app works correctly as a standalone tool will help with implementation with other shiny apps.
#If the app is working correctly by itself it should work when incorporated with another app. 
#To include this app with another we just need to include a few commands to our existing app to implement our new calculator. 

appUI <- navbarPage(
  
  title = "CellTag Simulator",
  
  theme = shinythemes::shinytheme(theme = "united"),
  
  #Include this command in the UI function of the existing app.
  calculatorPanel("calc")
  
)

appServer <- function(input, output, session){
  
  #Include this command in the server function of the existing app.
  #userInputs may need to be changed to ensure the correct inputs are being passed.
  #This will mainly depend on the app being added too. 
  shiny::callModule(
    
    module = calculatorServer,
    id = "calc",
    userInputs = input
  )
}

#Finally we need to run a few more commands in the global environment of the app. 
#If you have a two file shiny app. These commands can be included in a new file named global.R
#This file should be in the same directory as the server.R and ui.R files. 
#The commands in the global.R file are run to initialize the global environment for the shiny app.
#If our test app above were turned into a two file app. The global file would look like this. 
#I will mark the commands we will probably need to place in the global.R file
#of the celltag.org shiny app.


#library(shiny)
#library(shinycssloaders) //This
#library(ggplot2)
#library(reticulate) //This
#
#source("modularCalcSimp.R) //This
#
#reticulate::source_python("duplicates.py") //This
#