SeuratWrapper <- function(SeuratObjName, logExpData, Label, NewMeta){
  
  SeuratObjName <- CreateSeuratObject(raw.data = logExpData)
  SeuratObjName <- FindVariableGenes(SeuratObjName, do.plot = F, display.progress = F)
  SeuratObjName <- ScaleData(SeuratObjName)
  SeuratObjName <- AddMetaData(SeuratObjName, NewMeta)
  hv.genes <- head(rownames(SeuratObjName@hvg.info), 1000)
  SeuratObjName <- RunPCA(SeuratObjName, 
                              pc.genes = hv.genes, 
                              do.print = FALSE)
  SeuratObjName <- FindClusters(SeuratObjName, 
                                    reduction.type = "pca", 
                                    dims.use = 1:10, 
                                    resolution = 1, 
                                    print.output = FALSE, 
                                    save.SNN = TRUE,
                                    force.recalc = T)
  SeuratObjName <- RunTSNE(SeuratObjName, dims.use = 1:10, do.fast = TRUE,check_duplicates = FALSE)
  
  pdf(paste(Label,".plots.pdf",sep=""),width=8,height = 8)
  PCAPlot(SeuratObjName, dim.1 = 1, dim.2 = 2)
  PCElbowPlot(SeuratObjName)
  TSNEPlot(SeuratObjName, do.label = TRUE)
  TSNEPlot(SeuratObjName, do.label = TRUE,group.by ="tech")
  dev.off()
  
  return(SeuratObjName)
}



FeaturePrep <- function(SeuratObj, gene.set, ScoreName){
  library(matrixStats)
  # Get mean expression of genes of interest per cell
  #mean.exp <- colMeans(x = SeuratObj@data[gene.set, ], na.rm = TRUE)
  mean.exp <- colMaxs(x = 2^as.matrix(SeuratObj@data[gene.set, ]), na.rm = TRUE)
  
  # Add mean expression values in 'object@meta.data$gene.set.score'
  if (all(names(x = mean.exp) == rownames(x = SeuratObj@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
        "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
    SeuratObj@meta.data[[ScoreName]] <- mean.exp
  }
  return(SeuratObj)
}


selecteGenes.best.loadings <- function(trainingExpData, pcs, num, caption="Highest loading"){
  #Second version:
  #p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom
  #Weighted gene picking depending on PC number: Initial PCs give more genes.
  #For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
  print("Performing PCA...")
  #Do PCA here on input data
  library(matrixStats)
  trainingExpData <- trainingExpData[, apply(trainingExpData, 2, var) != 0]
  trainingExpData <- trainingExpData[, which(colVars(as.matrix(trainingExpData)) > 0.05)]
  
  pcatrain <- prcomp(trainingExpData,center = FALSE,scale=FALSE,rank. = 20)
  save(pcatrain,file="PCA.train.data")
  
  print("Selecting the genes as best features...")
  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  pb <- txtProgressBar(min = 0, max = pcs, style = 3)
  for(i in 1:pcs){
    orderedpcai <- pcatrain$rotation[order(abs(pcatrain$rotation[,i]),decreasing = TRUE),i]
    Gn <- round((pcs-i+1)*(num*2)/(pcs*(pcs+1)))
    TotalGenes = as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai,Gn)),bestgenes=head(orderedpcai,Gn))
    load <- rbind(load, top)
    setTxtProgressBar(pb,i)
    cat("\n",'Picking the best genes from first', i,'PCs is completed.',"\n")
  }

  bestgenes <- unique(load$genes)
  return(bestgenes)
}


prepareDataset <- function(ExpressionData, CellLabels, do.splitTest=FALSE, percent.Testset=0.2, regenerate.data=FALSE, run.name){
  #Required: This will take a Normalized expression data matrix, rows as genes and columns as cells. Example: as.matrix(SeuratObject@data)
  #Required: A list of cell labels. same dimension as colnames(input expression). Example: SeuratObject@meta.data$res.1
  #This will have an option to split data into test and training datasets. Default is, 0.2, 20%.
  #Default percent for test is the 20% of the whole data.
  
  
  #Transpose the matrix cols <--> rows t()
  #Keep the data in matrix form, otherwise randomForest will throw error: 'Error: protect(): protection stack overflow'
  trainingData <- as.data.frame(t(ExpressionData))
  #It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
  names(trainingData) <- make.names(names(trainingData))
  trainingData$CellType <- CellLabels
  #Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
  trainingData$CellType <- factor(trainingData$CellType)
  
  
  tsub <- trainingData[,!(names(trainingData) %in% c("CellType"))]
  tsub <- droplevels(tsub)
  indx <- sapply(tsub, is.factor)
  tsub[indx] <- lapply(tsub[indx], function(x) as.numeric(as.character(x)))
  
  #trainingData <- cbind(tsub, trainingData[,"CellType"])
  #tmat <- as.matrix(trainingData)
  
  gc()
  options("expression" =500000)
  #if do.splitTest=TRUE, separate "percent.Testset"% of matrix as test set. save this test data as a different slot.
  if (do.splitTest == TRUE){
    
    if(missing(percent.Testset)){
      ratio.test=0.2
      }else{
      ratio.test=percent.Testset
      }
    
    #First split the data into training and test sets
    smp_size <- floor((1-ratio.test) * nrow(tsub))
    print(smp_size)
    ## set the seed to make your partition reproductible
    set.seed(123)
    #Randomly sample the cells (rows now)
    train_ind <- base::sample(seq_len(nrow(tsub)), size = smp_size)
    train <- tsub[train_ind, ]
    
    save(trainingData, file=paste(run.name,".trainingData.tmp.Rdata",sep = ""))
    #save(tsub,file="tsub.tmp.Rdata")
    rm(trainingData, tsub)
    #Perform PCA on data (optional: only on training portion) and select genes for features
    pca.genes <- as.character(selecteGenes.best.loadings(train,10,500))
    
    load(paste(run.name,".trainingData.tmp.Rdata",sep = ""))
    trainingData.postPCA <- trainingData[train_ind,c(pca.genes,"CellType")]
    trainingData.postPCA <- droplevels(trainingData.postPCA)
    
    #test <- tsub[-train_ind, pca.genes]
    test <- trainingData[-train_ind,c(pca.genes,"CellType")]
    #save the test dataset
    save(test, file=paste(run.name,".testing.data",sep = ""))
    # Then do PCA on only training set:
    }else{
      train <- tsub
      print("Not splitting data for test set...")
      save(trainingData, file=paste(run.name,".trainingData.tmp.Rdata",sep = ""))
      #save(tsub,file="tsub.tmp.Rdata")
      rm(trainingData, tsub)
      #Perform PCA on data (optional: only on training portion) and select genes for features
      pca.genes <- as.character(selecteGenes.best.loadings(train,10,2000))
      load(paste(run.name,".trainingData.tmp.Rdata",sep = ""))
      trainingData.postPCA <- trainingData[,c(pca.genes,"CellType")]
      trainingData.postPCA <- droplevels(trainingData.postPCA)
      
    }#closes the do.splitTest
  #save(train, file="Training.data")
  
  #Subset the data for selected genes.
  #trainfiltered <- train[, pca.genes]
  #trainingData=input data
  #trainingData.postPCA <- data.frame(trainfiltered, CellType=trainingData[rownames(trainfiltered),"CellType"])
  save(trainingData.postPCA, file=paste(run.name,".trainingData.postPCA.data",sep = ""))
  
}#closes the function

CellTyperTrainer <- function(trainingData, run.name){
  library(randomForest)
  modelname <- paste(run.name, "_rf",sep = "")
  modelname = randomForest(CellType~., data = trainingData, norm.votes = TRUE,importance=TRUE, proximity = TRUE, ntree=500)
  save(modelname, file=paste(run.name,".RF_model.Robj",sep = ""))
  print(rf)
}

CellTyper <- function(testSet, model){
  library(caret)
  library(randomForest)
  pred.test.prob <- predict(model, testSet, type = "prob")
  pred.test.prob <- as.data.frame(pred.test.prob)
  
  pred.test.out <- predict(model, testSet, type="response")
  
  testPred <- pred.test.prob
  testPred$Prior <- testSet$CellType
  
  testPred$Prediction <- pred.test.out
  
  confmat <- confusionMatrix(testPred$Prediction, testPred$Prior)
  confmat$table
  
  print(  
    confusionMatrix(data = testPred$Prediction,  
                    reference = testPred$Prior))
  
  model$confusion[,"class.error"]
  
  varImpPlot(model, sort = T, n.var=20,main="Top 10 - Variable Importance")
  return(testPred)
}