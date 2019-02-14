
# ==============================================================================
# This file contains function, or modules, used to construct the UI for the
# CellTagViz website. These modules will be used to also construct the UI for
# visualization of non-CellTagged data as well.
# ==============================================================================


welcomePanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "Welcome!",

    mainPanel(

      includeHTML(path = "./markdown/WELCOME.html")

    )
  )
}


# ==============================================================================
# ==============================================================================


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
          tabPanel("Network"),
          tabPanel("Pseudotime"),
          tabPanel("Stacked Bar Charts"),
          tabPanel("Scatter Plots"),
          tabPanel("Meta Data")

        ),

        verbatimTextOutput("input_out")

      )
    )
  )
}


# ==============================================================================
# ==============================================================================


dataPanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "Data",

    mainPanel()

  )
}


# ==============================================================================
# ==============================================================================


peoplePanelUI <- function(id){

  ns <- NS(id)

  tabPanel(

    title = "People",

    mainPanel(

      includeHTML("./markdown/PEOPLE.html")

    )
  )
}


# ==============================================================================
# ==============================================================================


geneChoiceUI <- function(id){

  ns <- NS(id)

  selectizeInput(

      inputId = "GENE",
        label = "Choose a Gene",
      choices = letters,
    #selectize = FALSE,
     selected = "Apoa1"

  )
}


# ==============================================================================
# ==============================================================================


metaChoiceUI <- function(id){

  ns <- NS(id)

  selectizeInput(

      inputId = "META",
        label = "Choose a Variable",
      choices = LETTERS,
    #selectize = FALSE,
     selected = "State.Monocle"

  )
}


# ==============================================================================
# ==============================================================================


cloneChoiceUI <- function(id){

  ns <- NS(id)

  selectizeInput(

      inputId = "CLONES",
        label = "Choose a Clone",
      choices = LETTERS,
    #selectize = TRUE,
     selected = NULL

  )
}


# ==============================================================================
# ==============================================================================


contourButtonUI <- function(id){

  ns <- NS(id)

  checkboxInput(

    inputId = "addContour",
      label = "Add Countour Lines",
      value = FALSE

  )
}


# ==============================================================================
# ==============================================================================


factorButtonUI <- function(id){

  ns <- NS(id)

  checkboxInput(

    inputId = "isFactor",
      label = "Is variable a factor?",
      value = FALSE

  )
}


# ==============================================================================
# ==============================================================================


plotDownloadUI <- function(id){

  ns <- NS(id)

  downloadButton(

    outputId = "plotDownload",
       label = "Download Plot"

  )
}


# ==============================================================================
# ==============================================================================


plotPanelSideBar <- function(id){

  ns <- NS(id)

  tagList(

    tsneSideBarUI(id),
    networkSideBarUI(id),
    pseudotimeSideBarUI(id),
    scatterSideBarUI(id),
    stackedSideBarUI(id),
    metaSideBarUI(id),
    plotDownloadUI(id)

  )
}


# ==============================================================================
# ==============================================================================


tsneSideBarUI <- function(id){

  ns <- NS(id)

  conditionalPanel(

    condition = "input.plotPanel == 'tSNE'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    contourButtonUI(id),
    factorButtonUI(id),
    helpText("tSNE Panel")

  )
}


# ==============================================================================
# ==============================================================================


networkSideBarUI <- function(id){

  ns <- NS(id)

  conditionalPanel(

    condition = "input.plotPanel == 'Network'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    cloneChoiceUI(id),
    helpText("Network Panel")

  )
}


# ==============================================================================
# ==============================================================================


pseudotimeSideBarUI <- function(id){

  ns <- NS(id)

  conditionalPanel(

    condition = "input.plotPanel == 'Pseudotime'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    contourButtonUI(id),
    factorButtonUI(id),
    helpText("Pseudotime Panel")

  )
}


# ==============================================================================
# ==============================================================================


stackedSideBarUI <- function(id){

  ns <- NS(id)

  conditionalPanel(

    condition = "input.plotPanel == 'Stacked Bar Charts'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    helpText("Stacked Bar Panel")

  )
}


# ==============================================================================
# ==============================================================================


scatterSideBarUI <- function(id){

  ns <- NS(id)

  conditionalPanel(

    condition = "input.plotPanel == 'Scatter Plots'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    factorButtonUI(id),
    helpText("Scatter Chart Panel")

  )
}


# ==============================================================================
# ==============================================================================


metaSideBarUI <- function(id){

  ns <- NS(id)

  conditionalPanel(

    condition = "input.plotPanel == 'Meta Data'",
    geneChoiceUI(id),
    helpText("Meta Data Panel")

  )
}



