





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
#' @return Returns a shiny UI element
#'

welcomePanelUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    title = "Welcome!",
    shiny::mainPanel(shiny::includeHTML(path = "inst/markdown/WELCOME.html"))
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
#' @param inputData SingleCellExperiment object which contains current data.
#'
#' @return Returns a shiny UI element
#'


plotsPanelUI <- function(id, inputData) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    title = "Plots",

    shiny::sidebarLayout(
      position = "left",
      fluid = TRUE,
      sidebarPanel = shiny::sidebarPanel(plotPanelSideBar(id, inputSCE = inputData)),

      shiny::mainPanel(
        shiny::tabsetPanel(
          id = id,
          type = "pills",

          shiny::tabPanel("tSNE"),
          shiny::tabPanel("Network"),
          shiny::tabPanel("Pseudotime"),
          shiny::tabPanel("Stacked Bar Charts"),
          shiny::tabPanel("Scatter Plots"),
          shiny::tabPanel("Meta Data")
        ),

        shiny::tagList(
          shiny::plotOutput(ns("testPlot")),
          shiny::verbatimTextOutput(ns("inputOut"))
        )
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
#' @param inputData SingleCellExperiment object which contains current data.
#'
#' @return Returns a shiny UI element
#'

dataPanelUI <- function(id, inputData) {
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
#' good place to include links to things like our lab webpage and twitter/
#' github/ social media profiles. This panel utilizes HTML files created using
#' RStudio and R markdown documents. Making these HTML files stand alone will
#' mean the photos do not need to be included with the package separately.
#'
#' @param id String which defines the namespace for the module
#'
#' @return Returns a shiny UI element
#'


peoplePanelUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    title = "People",

    shiny::mainPanel(shiny::includeHTML("inst/markdown/PEOPLE.html"))
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
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a shiny UI element
#'


featureChoiceUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::selectizeInput(
    inputId = ns("GENE"),
    label = "Choose a Gene",
    choices = c(Choose = "", rownames(inputSCE))
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
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a shiny UI element
#'


metaChoiceUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  if (id == "BARCHART-GROUP") {
    return(
      shiny::selectizeInput(
        inputId = ns("META"),
        label = "Choose a Variable",
        choices = names(SingleCellExperiment::colData(inputSCE)),
        selected = "res.0.6.Seurat"
      )
    )
  } else {
    return(
      shiny::selectizeInput(
        inputId = ns("META"),
        label = "Choose a Variable",
        choices = names(SingleCellExperiment::colData(inputSCE)),
        selected = "State.Monocle"
      )
    )
  }
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
#' @return Returns a shiny UI element
#'


cloneChoiceUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::selectizeInput(
    inputId = ns("CLONES"),
    label = "Choose a Clone",
    choices = LETTERS,
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
#' @return Returns a shiny UI element
#'


contourButtonUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::checkboxInput(
    inputId = ns("addContour"),
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
#' @return Returns a shiny UI element
#'


factorButtonUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::checkboxInput(
    inputId = ns("isFactor"),
    label = "Is meta data categorical?",
    value = FALSE
  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines a button to download the current plot as a PDF.
#'
#' This module is used to define the UI of the plot download button for
#' CellTagViz. This module is a button which when clicked will download and save
#' the current plot for the user. Slight differences in the downloaded plot and
#' the plot produced by CellTagViz might occur. This is due to the CellTagViz
#' interactive plot being rendered via plotly while the downloaded plot is
#' rendered with ggplot2. This is a simple module which can be used to download
#' many visualizations.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a shiny UI element
#'


plotDownloadUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::downloadButton(
    outputId = ns("plotDownload"),
    label = "Download Plot"
  )
}


# ==============================================================================
# ==============================================================================


#' Module which defines the sidebar of the Plot Panel.
#'
#' This module is used to define the UI of the Plot Panel sidebar for
#' CellTagViz. This module is a taglist which contains calls to multiple other
#' modules. The function calls this module makes define conditional sidebar
#' layouts for each of the tabs in the tabset on the Plot Panel. By defining
#' multiple sidebar layouts with the function calls allows each tab to have a
#' unique sidebar. This module also calls the function \code{plotDownloadUI}
#' this was done in this module as I could not produce proper formatting when
#' placing the module call inside the modules which define the conditional
#' sidebars. Furthermore, because the option to download plots is available for
#' each of the plots the download button is not neccesarrily conditional. This
#' makes it possible to call the module \code{plotDownload} in this module and
#' generate the desired layout. Sidebars for new plot tabs can be easily added
#' by creating a new module defining the sidebar and placing the function call
#' in this module.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a shiny UI element
#'


plotPanelSideBar <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::tagList(
    tsneSideBarUI("TSNE", inputSCE = inputSCE),
    networkSideBarUI("NETWORK", inputSCE = inputSCE),
    pseudotimeSideBarUI("PSEUDO", inputSCE = inputSCE),
    scatterSideBarUI("SCATTER", inputSCE = inputSCE),
    stackedSideBarUI("BARCHART", inputSCE = inputSCE),
    metaSideBarUI("METAVAR", inputSCE = inputSCE),
    plotDownloadUI("DOWNLOAD", inputSCE = inputSCE)
  )
}


# ==============================================================================
# ==============================================================================


#' This module defines the UI element for the tSNE panel Sidebar.
#'
#' This module defines a conditional panel which is displayed when the tSNE plot
#' panel is displayed. It makes calls to modules which allow the user to choose
#' a meta data variable to color the plot, choose a feature to visualize the
#' expression, add contour lines to the plot, and ensure categorical variables
#' are treated as such. The module also defines help text which is rendered in
#' the panel. This is used to ensure the correct sidebars are displayed on the
#' correct panels and is only used during testing. Importantly this panel is
#' set to display only in the correct circumstances. This sidebar will only be
#' displayed when \code{input$plotPanel} is equal to "tSNE". This module can be
#' easily modified to define new sidebar panels for new plots.
#'
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a Shiny UI element
#'
#'



tsneSideBarUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::conditionalPanel(
    condition = "input.plots == 'tSNE'",
    reductionChoiceUI(id, inputSCE = inputSCE),
    featureChoiceUI(id, inputSCE = inputSCE),
    metaChoiceUI(id, inputSCE = inputSCE),
    contourButtonUI(id),
    factorButtonUI(id),
    shiny::helpText("tSNE Panel")
  )
}


# ==============================================================================
# ==============================================================================


#' This module defines the UI element for the Network panel Sidebar.
#'
#' This module defines a conditional panel which is displayed when the Network
#' plot panel is displayed. It makes calls to modules which allow the user to
#' choose a meta data variable to color the plot, choose a feature to visualize
#' the expression, and choose clones to visualize. The module also defines help
#' text which is rendered in the panel. This is used to ensure the correct
#' sidebars are displayed on the correct panels and is only used during testing.
#' Importantly this panel is set to display only in the correct circumstances.
#' This sidebar will only be displayed when \code{input$plotPanel} is equal to
#' "Network". This module can be easily modified to define new sidebar panels
#' for new plots.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a Shiny UI element
#'


networkSideBarUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::conditionalPanel(
    condition = "input.plots == 'Network'",
    featureChoiceUI(id, inputSCE = inputSCE),
    metaChoiceUI(id, inputSCE = inputSCE),
    cloneChoiceUI(id),
    shiny::helpText("Network Panel")
  )
}


# ==============================================================================
# ==============================================================================


#' This module defines the UI element for the Pseudotime plot panel Sidebar.
#'
#' This module defines a conditional panel which is displayed when the
#' Pseudotime plot panel is displayed. It makes calls to modules which allow the
#' user to choose a meta data variable to color the plot, choose a feature to
#' visualize the expression, add contour lines to the plot, and ensure
#' categorical variables are treated as such. The module also defines help text
#' which is rendered in the panel. This is used to ensure the correct sidebars
#' are displayed on the correct panels and is only used during testing.
#' Importantly this panel is set to display only in the correct circumstances.
#' This sidebar will only be displayed when \code{input$plotPanel} is equal to
#' "Pseudotime". This module can be easily modified to define new sidebar panels
#' for new plots.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a Shiny UI element
#'


pseudotimeSideBarUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::conditionalPanel(
    condition = "input.plots == 'Pseudotime'",
    featureChoiceUI(id, inputSCE = inputSCE),
    metaChoiceUI(id, inputSCE = inputSCE),
    contourButtonUI(id),
    factorButtonUI(id),
    shiny::helpText("Pseudotime Panel")
  )
}


# ==============================================================================
# ==============================================================================


#' This module defines the UI element for the Stacked Bar Chart panel Sidebar.
#'
#' This module defines a conditional panel which is displayed when the Stacked
#' bar chart plot panel is displayed. It makes calls to modules which allow the
#' user to choose a meta data variable to color the plot, choose a feature to
#' visualize the expression. The module also defines help text which is rendered
#' in the panel. This is used to ensure the correct sidebars are displayed on
#' the correct panels and is only used during testing. Importantly this panel is
#' set to display only in the correct circumstances. This sidebar will only be
#' displayed when \code{input$plotPanel} is equal to "Stacked Bar Charts". This
#' module can be easily modified to define new sidebar panels for new plots.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a Shiny UI element
#'


stackedSideBarUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::conditionalPanel(
    condition = "input.plots == 'Stacked Bar Charts'",
    metaChoiceUI("BARCHART-GROUP", inputSCE = inputSCE),
    metaChoiceUI(id, inputSCE = inputSCE),
    horizontalButtonUI(id),
    shiny::helpText("Stacked Bar Panel")
  )
}


# ==============================================================================
# ==============================================================================


#' This module defines the UI element for the Scatter Chart panel Sidebar.
#'
#' This module defines a conditional panel which is displayed when the Scatter
#' chart plot panel is displayed. It makes calls to modules which allow the
#' user to choose a meta data variable to color the plot, choose a feature to
#' visualize the expression, and ensure categorical variables are treated as
#' such. The module also defines help text which is rendered in the panel. This
#' is used to ensure the correct sidebars are displayed on the correct panels
#' and is only used during testing. Importantly this panel is set to display
#' only in the correct circumstances. This sidebar will only be displayed when
#' \code{input$plotPanel} is equal to "Scatter Plots". This module can be
#' easily modified to define new sidebar panels for new plots.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a Shiny UI element
#'


scatterSideBarUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::conditionalPanel(
    condition = "input.plots == 'Scatter Plots'",
    featureChoiceUI(id, inputSCE = inputSCE),
    metaChoiceUI(id, inputSCE = inputSCE),
    factorButtonUI(id),
    shiny::helpText("Scatter Chart Panel")
  )
}


# ==============================================================================
# ==============================================================================


#' This module defines the UI element for the Meta Data panel Sidebar.
#'
#' This module defines a conditional panel which is displayed when the Meta Data
#' panel is displayed. It makes calls to modules which allow the user to choose
#' a feature to visualize the expression of. The module also defines help text
#' which is rendered in the panel. This is used to ensure the correct sidebars
#' are displayed on the correct panels and is only used during testing.
#' Importantly this panel is set to display only in the correct circumstances.
#' This sidebar will only be displayed when \code{input$plotPanel} is equal to
#' "Meta Data". This module can be easily modified to define new sidebar panels
#' for new plots.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a Shiny UI element
#'


metaSideBarUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::conditionalPanel(
    condition = "input.plots == 'Meta Data'",
    featureChoiceUI(id, inputSCE = inputSCE),
    shiny::helpText("Meta Data Panel")
  )
}


#' Module which defines a checkbox to make stacked bar charts horizontal.
#'
#' This module is used to define the UI of the factor option checkbox for
#' CellTagViz. This module is a checkbox which when ticked will consider the
#' given meta data variable as a factor. In turn the variable is considered
#' categorical and treated as such. This is a simple module which can be
#' used for many visualizations.
#'
#' @param id String which defines the namespace for the module
#'
#' @return Returns a shiny UI element
#'


horizontalButtonUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::checkboxInput(
    inputId = ns("plotFlip"),
    label = "Horizontal Plot",
    value = FALSE
  )
}


#' Module which defines a selection input for choosing a dimension reduction
#' method.
#'
#' This module is used to define the UI of the dimension reduction selection
#' input for CellTagViz. This module is a drop down list which allows ths user
#' to choose which dimension reduction method to visualize. The choices for this
#' list are the names of the datasets returned by \code{reducedDimNames}.
#' This is a simple module which can be used for many visualizations.
#'
#' @param id String which defines the namespace for the module
#'
#' @param inputSCE SingleCellExperiment object which contains current data.
#'
#' @return Returns a shiny UI element
#'


reductionChoiceUI <- function(id, inputSCE) {
  ns <- shiny::NS(id)

  shiny::selectizeInput(
    inputId = ns("REDUCTION"),
    label = "Choose a dimension reduction method",
    choices = SingleCellExperiment::reducedDimNames(inputSCE)
  )
}
