% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sceTools.R
\name{makePlotData}
\alias{makePlotData}
\title{Generates ggplot2 compatible data frames.}
\usage{
makePlotData(sce, redMethod, metaVar, feature = FALSE, cells = FALSE,
  exprDat = FALSE)
}
\arguments{
\item{sce}{SingleCellExperiment Object}

\item{redMethod}{String Name of the dimension reduction method to visualize.}

\item{metaVar}{String Name of the column to return from meta/col Data.}

\item{feature}{String Name of a feature to return from assay data.}

\item{cells}{Character Vector of Cell Barcodes used to subset data.}

\item{exprDat}{YES}
}
\value{
This function returns a data frame which can be used to plot and
color cells from an SCE object using ggplot2.
}
\description{
\code{makePlotData} creates data frames for use with ggplot2 visualizations.
This is accomplished by calling two functions \code{getEmbeddings} and
\code{getMetaData}. The objects returned by both of these functions are
then merged into one single data frame. This combined data frame contains the
coordinates of each cell from the given dimension reduction method along with
the given meta data which can then be used for grouping and coloring.
}
\examples{

\dontrun{

metaData <- sample(letters, 15)

tSNE.embeddings <- matrix(
  data = rnorm(250),
  nrow = 25,
  ncol = 10
)

sce <- SingleCellExperiment::SingleCellExperiment(colData = list(
  Letters = metaData
))

reducedDim(sce, "tSNE") <- tSNE.embeddings

plotData <- makePlotData(sce, "tSNE", "Letters")
}

}
