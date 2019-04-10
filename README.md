# The Package formerly known as CellTagViz

This repo is the home to an R package which contains functions to visualize
single cell genomic datasets. The package consits mainly of Shiny server and ui
modules to create interactive Shiny apps as well as functions to generate the 
desired visualizations with ggplot2. Currently it is compatible with Seurat V2, 
Seurat V3, and Monocle datasets.


# Installation

The package can be installed from github using the package devtools.

```{r}
devtools::install_github(repo = "morris-lab/CellTagViz-package")

```

# Visualizing Data

Now that the package is installed you can load your Seurat or Monocle data and
visualize it using the package.

```{r}
library(CellTagViz)

#Load your Seurat or Monocle Objects

dataList <- list("SeuratV2" = your_seuratv2_data, 
                 "SeuratV3" = your_seuratv3_data, 
                  "Monocle" = your_monocle_cds)

vizData <- combineSCE(dataList)

TestApp(vizData = vizData)

```
When creating the vizData object the package attempts to preserve all dimension
reduction, expression data, and meta data from each object. If visualizing
dimension reduction data from a Seurat object I suggest using the Scaled Data
expression data from the same object. For visualizing pseudotime from monocle
it is probably best to use the Monocle.Counts expression data. 

