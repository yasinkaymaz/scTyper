.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")
source("/n/home13/yasinkaymaz/codes/scTyper/RFTyper/RF_functions.R")
library(tidyverse)
library(Seurat)
setwd("/n/home13/yasinkaymaz/codes/test/brainCellMarkers/")
#Load multiple PBMC datasets

dir='/n/home13/yasinkaymaz/codes/test/brainCellMarkers/Datasets/'

datasets <- c(
  "pbmc4K_CR.2.1.0",
  "pbmc8K_CR.2.1.0",
  "pbmc33K_CR.1.1.0",
  "pbmc3K_CR.1.1.0",
  "pbmc68K_CR.1.0.0",
  "pbmc6K_CR.1.1.0"
)

for (d in datasets){
  print(d)
  exp.data <- Read10X(data.dir = paste(dir,d,"/",sep = ""))
  SeuratWrapper(d, exp.data, Normalize = T, Label = d)
  rm(exp.data)
}



# Model Training
load(paste(dir,"/pbmc3k_final.Rda",sep = ""))
sampleExpdata <- as.matrix(pbmc@data)
sampleClassdata <- pbmc@meta.data$ClusterNames_0.6

run.name = "pbmc3K"
prepareDataset(ExpressionData = sampleExpdata, CellLabels = sampleClassdata, do.splitTest = F, run.name = run.name)

rm(sampleExpdata)
load(paste(run.name,".trainingData.postPCA.data",sep=""))

CellTyperTrainer(trainingData.postPCA, run.name)

load(paste(run.name,".RF_model.Robj", sep = ""))
rf.strata <- modelname
rm(modelname)
head(pbmc@meta.data)


# Testing
ob.list <- list(pbmc4K_CR.2.1.0, pbmc8K_CR.2.1.0, pbmc33K_CR.1.1.0, pbmc3K_CR.1.1.0, pbmc68K_CR.1.0.0, pbmc6K_CR.1.1.0)

pbmc122K.integrated <- SeuratCCAmerger(listofObjects = ob.list)

pbmc122K.integrated <- CellTyper(SeuratObject = pbmc122K.integrated, model = rf.strata)

PlotPredictions(SeuratObject = pbmc122K.integrated, model = rf.strata, outputFilename = "PBMC122k.Predictions")


