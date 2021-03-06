% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyModules.R
\name{plotsPanelMinimalUI}
\alias{plotsPanelMinimalUI}
\title{Module for plots panel of CellTagViz}
\usage{
plotsPanelMinimalUI(id, inputData)
}
\arguments{
\item{id}{String which defines the namespace of the module.}

\item{inputData}{SingleCellExperiment object which contains current data.}
}
\value{
Returns a shiny UI element
}
\description{
This is a shiny UI module that defines the layout of the plots panel of the
navbar menu. For the sidebar of this panel layout the module calls the
function \code{plotPanelSideBar} to define the different sidebar inputs for
each tab. This module can be modified to include more/different plots. In the
main panel it is possible to print the current user inputs. This can be
helpful when modifying the shiny app. The function to do this has been
commented out. The plots are defined as tabs in a tabset in the main panel.
New plots can easily be added this way.
}
