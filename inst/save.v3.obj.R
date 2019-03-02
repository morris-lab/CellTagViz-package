


load("~/../Desktop/sham.sbr.integrated.RData")


saveRDS(object = integrated.subset, file = "~/../Desktop/integrated.subset.RDS")



path <- "~/../Desktop/integrated.subset.RDS"

dat <- readRDS(file = path)

monPath <- "~/../Desktop/unsupervised timeline all data.RDS"

mon <- readRDS(file = monPath)


#monSCE <- exportCDS(mon, export_to = "Scater")

library(S4Vectors)

source("~/GitHub/working.Viz/R/dataTools.R")

monSCE <- addMonocleDataNew(mon)

sce <- addSeuratv3(dat)
