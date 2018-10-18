
source("~/Dropbox/codes/scTyper/RFTyper/RF_functions.R")
library(tidyverse)
library(Seurat)
setwd("~/Documents/RFTyper/test3/")
#Load multiple PBMC datasets

dir='~/Documents/RFTyper/Datasets/'

datas <- c(
"CD14posMonocytes",
"CD19posBCells",
"CD34posCells",
"CD4posHelperTCells",
"CD4posCD25posRegulatoryTCells",
"CD4posCD45RAposCD25negNaiveTcells",
"CD4posCD45ROposMemoryTCells",
"CD56posNaturalKillerCells",
"CD8posCytotoxicTcells",
"CD8posCD45RAposNaiveCytotoxicTCells"
)

for (d in datas){
  print(d)
  exp.data <- Read10X(data.dir = paste(dir,d,"/",sep = ""))
  SeuratWrapper(d, exp.data, Normalize = T, Label = d)
  rm(exp.data)
}
