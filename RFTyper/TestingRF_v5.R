
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


merged <- ScaleData(merged,genes.use = hv.genes)

merged <- RunPCA(merged,
                 pc.genes = hv.genes,
                 do.print = FALSE)

merged <- FindClusters(merged,
                       reduction.type = "pca",
                       dims.use = 1:10,
                       resolution = 1,
                       print.output = FALSE,
                       save.SNN = TRUE,
                       force.recalc = T)
merged <- RunTSNE(merged, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)

head(merged@meta.data)

load("rf.with.strata.RF_model.Robj")
rf <- modelname

merged <- CellTyper(SeuratObject = merged, model = rf)
PlotPredictions(SeuratObject = merged, model = rf, outputFilename = "merged.CD14-CD19-CD34.predictions")

TSNEPlot(merged, do.label = T)



#optimize rf model

load("rf.with.strata.RF_model.Robj")
rf <- modelname

print(rf)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

filter_p = 0.05
bestvote_cutoff = 0.7
badinput_ids <- NULL
round_n=1
badinput_stats <- data.frame()

rfvotes <- as.data.frame(rf$votes)
rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
rfvotes$label <- apply(rf$votes, 1, function(x)  names(x)[which.max(x)] )
rfvotes$inputnames <- rownames(rfvotes)
e <- as.data.frame(table(rfvotes$label))

Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()
# 
# for(i in 1:length(e$Var1)){
#   badinputs <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% select(inputnames) %>% c()
#   badinputs.m <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% summarize(mean=mean(bestvote)) %>% c()
#   badinput_ids <- unique(c(badinput_ids, badinputs$inputnames))
#   classBestscore <- rfvotes %>% as.tibble() %>%  filter(., label == e$Var1[i]) %>% summarize(median=median(bestvote)) %>% c()
#   print(paste("The number of input dropped for",e$Var1[i],"class is", length(badinputs$inputnames),sep = " "))
#   class_badinput_stats <- data.frame(class=e$Var1[i], classInputSize=e$Freq[i], allBestscoreMedian=Bestscore$median, classBestscoreMedian=classBestscore$median,  tossedInput=length(badinputs$inputnames), tossedBestvoteMean=badinputs.m$mean, iteration=round_n)
#   badinput_stats <- rbind(badinput_stats, class_badinput_stats)
# }
# 
# badinput_stats[is.nan(badinput_stats)] <- 0
# toss_n <- badinput_stats %>% as.tibble() %>% filter(., iteration == round_n) %>% summarise(n=sum(tossedInput)) %>% c()
# toss_n <- toss_n$n

#Repeat the process until the Median Best score improves
#currentscore <- badinput_stats %>% as.tibble() %>% filter(., iteration == round_n) %>% select(allBestscoreMedian) %>% c()
#currentscore <- currentscore$allBestscoreMedian[1]

#Defaults

currentscore = Bestscore$median
toss_n = dim(rfvotes)[1]

#While Median Best votes score overall is larger in the new iteration AND the number of inputs need to be tossed is larger %1 of all input data, then continue.
while (Bestscore$median >= currentscore && toss_n > round(0.01*dim(rfvotes)[1]) ){

  print(paste("Round number ",round_n))
  print(paste("Current score is", currentscore,". toss_n is", toss_n, ". Fractions is", round(0.01*dim(rfvotes)[1])))

  for(i in 1:length(e$Var1)){
    badinputs <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% select(inputnames) %>% c()
    badinputs.m <- rfvotes %>% as.tibble() %>% filter(., label == e$Var1[i] & bestvote < bestvote_cutoff) %>% arrange(bestvote) %>% top_n(n=round(filter_p*e$Freq[i]), wt=-bestvote) %>% summarize(mean=mean(bestvote)) %>% c()
    badinput_ids <- unique(c(badinput_ids, badinputs$inputnames))
    classBestscore <- rfvotes %>% as.tibble() %>%  filter(., label == e$Var1[i]) %>% summarize(median=median(bestvote)) %>% c()
    print(paste("The number of input dropped for",e$Var1[i],"class is", length(badinputs$inputnames),sep = " "))
    class_badinput_stats <- data.frame(class=e$Var1[i], classInputSize=e$Freq[i], allBestscoreMedian=Bestscore$median, classBestscoreMedian=classBestscore$median, tossedInput=length(badinputs$inputnames), tossedBestvoteMean=badinputs.m$mean, iteration=round_n)
    badinput_stats <- rbind(badinput_stats, class_badinput_stats)
  }
  
  badinput_stats[is.nan(badinput_stats)] <- 0
  toss_n <- badinput_stats %>% as.tibble() %>% filter(., iteration == round_n) %>% summarise(n=sum(tossedInput)) %>% c()
  toss_n <- toss_n$n
  
  print(badinput_stats)
  
  #filter input using the bad input list generated in the previous iteration
  trainingData.postPCA <- trainingData.postPCA[which(!rownames(trainingData.postPCA) %in% badinput_ids ),]
  
  #run the RF again with the updated training set
  rf <- randomForest(CellType~., data = trainingData.postPCA, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500, sampsize=c(table(trainingData.postPCA$CellType)))
  print(rf)
  #Check again to see if there is room to improve:
  rfvotes <- as.data.frame(rf$votes)
  rfvotes$bestvote <- apply(rfvotes, 1, function(x) max(x))
  rfvotes$label <- apply(rf$votes, 1, function(x)  names(x)[which.max(x)] )
  rfvotes$inputnames <- rownames(rfvotes)
  e <- as.data.frame(table(rfvotes$label))
  Bestscore <- rfvotes %>% as.tibble() %>% summarize(median=median(bestvote)) %>% c()
  #update the round
  round_n = round_n + 1
}

