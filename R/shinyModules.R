

welcomePanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "Welcome!",

    mainPanel(

      includeHTML("markdown/WELCOME.html")

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

          id = "plotPanel",

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

      includeHTML("markdown/PEOPLE.html")

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

  tagList(

    tsneSideBarUI(id),

    networkSideBarUI(id),

    pseudotimeSideBarUI(id)

  )

}


tsneSideBarUI <- function(id){

  conditionalPanel(

    condition = "input.plotPanel == 'tSNE'",

    geneChoiceUI(id),

    metaChoiceUI(id),

    helpText("tSNE Panel")

  )

}


networkSideBarUI <- function(id){

  conditionalPanel(

    condition = "input.plotPanel == 'Network'",

    geneChoiceUI(id),

    metaChoiceUI(id),

    helpText("Network Panel")

  )

}


pseudotimeSideBarUI <- function(id){

  conditionalPanel(

    condition = "input.plotPanel == 'Pseudotime'",

    geneChoiceUI(id),

    metaChoiceUI(id),

    helpText("Pseudotime Panel")

  )

}

