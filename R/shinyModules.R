
# ==============================================================================
# This file contains functions, or modules, used to construct the UI for the
# CellTagViz website. These modules will be used to also construct the UI for
# visualization of non-CellTagged data as well. Initially these modules were
# developed to construct the CellTagViz app, but they can be repurposed as well.
# ==============================================================================


#' Module for Welcome page of CellTag Viz
#'
#' This is a UI module for use with the CellTagViz shiny app.
#' The function includes an HTML file constructed using RStudio and R Markdown
#' files. This module is useful for manipulating and changing the welcome page
#' for the CellTagViz app. This function can also be replaced to change the
#' Welcome page. In the past there has been issues with including HTML files
#' as there were conflicts with other functionalities. But if there are no
#' issues including HTML files I feel they are probably the best to use as the
#' HTML files themselves can be rendered as stand alone documents.
#'
#' @param id String passed to module to identify the namespace of the module.
#'

welcomePanelUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tabPanel(

    title = "Welcome!",

    shiny::mainPanel(

      shiny::includeHTML(path = "./markdown/WELCOME.html")

    )
  )
}


# ==============================================================================
# ==============================================================================


#' Module for plots panel of CellTagViz
#'
#' This is a shiny UI module that defines the layout of the plots panel of the
#' navbar menu. For the sidebar of this panel layout the module calls the
#' function \code{plotPanelSideBar} to define the different sidebar inputs for
#' each tab. This module can be modified to include more/different plots. In the
#' main panel it is possible to print the current user inputs. This can be
#' helpful when modifying the shiny app. The function to do this has been
#' commented out. The plots are defined as tabs in a tabset in the main panel.
#' New plots can easily be added this way.
#'
#' @param id String which defines the namespace of the module.
#'


plotsPanelUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tabPanel(

    title = "Plots",

    shiny::sidebarLayout(

          position = "left",
             fluid = TRUE,
      sidebarPanel = shiny::sidebarPanel(

        plotPanelSideBar(id)

      ),

    shiny::mainPanel(

        shiny::tabsetPanel(

            id = "plotPanel",
          type = "pills",

          shiny::tabPanel("tSNE"),
          shiny::tabPanel("Network"),
          shiny::tabPanel("Pseudotime"),
          shiny::tabPanel("Stacked Bar Charts"),
          shiny::tabPanel("Scatter Plots"),
          shiny::tabPanel("Meta Data")

        ),

        shiny:: verbatimTextOutput("input_out")

      )
    )
  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines the Data panel.
#'
#' This module is used to define the UI of the data panel for CellTagViz.
#' This panel will include the option to download the data used by the shiny
#' app. This is also a good place to include links to thinkgs like GEO or the
#' protocols.io webpage.
#'
#' @param id String which defines the namespace for the module
#'

dataPanelUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tabPanel(

    title = "Data",

    shiny::mainPanel()

  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines the people panel.
#'
#' This module is used to define the UI of the people panel for CellTagViz.
#' This panel includes the pictures of the Morris lab members. This is also a
#' good place to include links to things like our lab webpage and twitter/github/
#' social media profiles. This panel utilizes HTML files created using RStudio
#' and R markdown documents. Making these HTML files stand alone will mean
#' the photos do not need to be included with the package separately.
#'
#' @param id String which defines the namespace for the module
#'


peoplePanelUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tabPanel(

    title = "People",

    shiny::mainPanel(

      shiny::includeHTML("./markdown/PEOPLE.html")

    )
  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines feature selection choices for user input.
#'
#' This module is used to define the UI of the feature selection user input for
#' CellTagViz. This module is a drop down selection box of the list of features
#' present in the data being visualized. These features can then be visualized.
#' This is a simple module which can be used in many situations.
#'
#' @param id String which defines the namespace for the module
#'


geneChoiceUI <- function(id){

  ns <- shiny::NS(id)

  shiny::selectizeInput(

      inputId = "GENE",
        label = "Choose a Gene",
      choices = letters,
    #selectize = FALSE,
     selected = "Apoa1"

  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines meta data variable  selection choices for user input.
#'
#' This module is used to define the UI of the meta data selection user input
#' for CellTagViz. This module is a drop down selection box of the list of
#' meta data variables which can be visualized. This is a simple module which
#' can be used in many situations.
#'
#' @param id String which defines the namespace for the module
#'


metaChoiceUI <- function(id){

  ns <- shiny::NS(id)

  shiny::selectizeInput(

      inputId = "META",
        label = "Choose a Variable",
      choices = LETTERS,
    #selectize = FALSE,
     selected = "State.Monocle"

  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines clone selection choices for user input.
#'
#' This module is used to define the UI of the clone selection user input for
#' CellTagViz. This module is a drop down selection box of the list of clones
#' present in the data being visualized. Using the selection from this input the
#' data will be subset in order to only visualize the cells in the selected
#' clones. This is a simple module which can be used in many situations.
#'
#' @param id String which defines the namespace for the module
#'


cloneChoiceUI <- function(id){

  ns <- shiny::NS(id)

  shiny::selectizeInput(

      inputId = "CLONES",
        label = "Choose a Clone",
      choices = LETTERS,
    #selectize = TRUE,
     selected = NULL

  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines a checkbox to add contour lines.
#'
#' This module is used to define the UI of the contour option checkbox for
#' CellTagViz. This module is a checkbox which add contour lines to the current
#' plot when ticked. This is a simple module which can be
#' used for many visualizations.
#'
#' @param id String which defines the namespace for the module
#'


contourButtonUI <- function(id){

  ns <- shiny::NS(id)

  shiny::checkboxInput(

    inputId = "addContour",
      label = "Add Countour Lines",
      value = FALSE

  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines a checkbox to consider meta data variables as factors.
#'
#' This module is used to define the UI of the factor option checkbox for
#' CellTagViz. This module is a checkbox which when ticked will consider the
#' given meta data variable as a factor. In turn the variable is considered
#' categorical and treated as such. This is a simple module which can be
#' used for many visualizations.
#'
#' @param id String which defines the namespace for the module
#'


factorButtonUI <- function(id){

  ns <- shiny::NS(id)

  shiny::checkboxInput(

    inputId = "isFactor",
      label = "Is meta data categorical?",
      value = FALSE

  )
}


# ==============================================================================
# ==============================================================================


plotDownloadUI <- function(id){

  ns <- shiny::NS(id)

  shiny::downloadButton(

    outputId = "plotDownload",
       label = "Download Plot"

  )
}


# ==============================================================================
# ==============================================================================


plotPanelSideBar <- function(id){

  ns <- shiny::NS(id)

  shiny::tagList(

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

  ns <- shiny::NS(id)

  shiny::conditionalPanel(

    condition = "input.plotPanel == 'tSNE'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    contourButtonUI(id),
    factorButtonUI(id),
    shiny::helpText("tSNE Panel")

  )
}


# ==============================================================================
# ==============================================================================


networkSideBarUI <- function(id){

  ns <- shiny::NS(id)

  shiny::conditionalPanel(

    condition = "input.plotPanel == 'Network'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    cloneChoiceUI(id),
    shiny::helpText("Network Panel")

  )
}


# ==============================================================================
# ==============================================================================


pseudotimeSideBarUI <- function(id){

  ns <- shiny::NS(id)

  shiny::conditionalPanel(

    condition = "input.plotPanel == 'Pseudotime'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    contourButtonUI(id),
    factorButtonUI(id),
    shiny::helpText("Pseudotime Panel")

  )
}


# ==============================================================================
# ==============================================================================


stackedSideBarUI <- function(id){

  ns <- shiny::NS(id)

  shiny::conditionalPanel(

    condition = "input.plotPanel == 'Stacked Bar Charts'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    shiny::helpText("Stacked Bar Panel")

  )
}


# ==============================================================================
# ==============================================================================


scatterSideBarUI <- function(id){

  ns <- shiny::NS(id)

  shiny::conditionalPanel(

    condition = "input.plotPanel == 'Scatter Plots'",
    geneChoiceUI(id),
    metaChoiceUI(id),
    factorButtonUI(id),
    shiny::helpText("Scatter Chart Panel")

  )
}


# ==============================================================================
# ==============================================================================


metaSideBarUI <- function(id){

  ns <- shiny::NS(id)

  shiny::conditionalPanel(

    condition = "input.plotPanel == 'Meta Data'",
    geneChoiceUI(id),
    shiny::helpText("Meta Data Panel")

  )
}



