

welcomePanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "Welcome!",

    mainPanel(

      includeHTML("../markdown/WELCOME.html")

    )


  )

}


plotsPanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "Plots",

    sidebarLayout(

      position = "left",

      fluid = TRUE,

      sidebarPanel = sidebarPanel(

        plotPanelSideBar(id)

      ),

      mainPanel(

        tabsetPanel(

          type = "pills",

          tabPanel("tSNE"),

          tabPanel("Pseudotime"),

          tabPanel("Network")

      )

    )

  )
  )
}



dataPanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "Data",

    mainPanel()

  )

}


peoplePanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "People",

    mainPanel(

      includeHTML("../markdown/PEOPLE.html")

    )

  )

}

geneChoiceUI <- function(id){

  ns <- NS(id)

  tagList(

    selectizeInput(ns("GENE"), "Choose a Gene", letters)

  )

}

metaChoiceUI <- function(id){

  ns <- NS(id)

  tagList(selectizeInput(ns("META"), "Choose a Variable", LETTERS)

    )

}

plotPanelSideBar <- function(id){

  ns <- NS(id)


  conditionalPanel( F,

    selectizeInput(ns("META"), "Choose a Variable", LETTERS),

    selectizeInput(ns("GENE"), "Choose a Gene", letters)

  )

}
