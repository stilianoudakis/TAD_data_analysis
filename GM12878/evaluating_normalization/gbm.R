#GBM

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)


setwd("/home/stilianoudakisc/TAD_data_analysis/comparing_normalization/")

#Chromosome 1
chr1data_f <- readRDS("chr1data_f.rds")

#Full data
#logitdata_f <- readRDS("logitdata_f.rds")

#set number of bootstrap samples
bootsamps = 10

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(chr1data_f$y[which(chr1data_f$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1data_f$y=="No"),
                        length(which(chr1data_f$y=="Yes")),
                        replace = TRUE)
}

#set.seed(123)
#for(j in 1:bootsamps){
#  sampids[,j] <- sample(which(logitdata_f$y==0),
#                        length(which(logitdata_f$y==1)),
#                        replace = TRUE)
#}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

########################################################################################

#With log transform

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
gbmlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1data_f$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              fpr <- matrix(nrow=ceiling((length(which(chr1data_f$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              auc <- numeric(bootsamps),
              varimp <- matrix(nrow=dim(chr1data_f)[2]-1,
                               ncol=bootsamps))
rownames(gbmlst[[4]]) <- colnames(chr1data_f)[-1]

#gbmlst <- list(tpr <- matrix(nrow=ceiling((length(which(logitdata_f$y==1))*2)*.3), 
#                            ncol=bootsamps),
#              fpr <- matrix(nrow=ceiling((length(which(logitdata_f$y==1))*2)*.3), 
#                            ncol=bootsamps),
#              auc <- numeric(bootsamps),
#              varimp <- matrix(nrow=dim(logitdata_f)[2]-1,
#                               ncol=bootsamps))
#rownames(gbmlst[[4]]) <- colnames(logitdata_f)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1data_f[which(chr1data_f$y=="Yes"),],
                           chr1data_f[sampids[,i],])
  
  #data <- rbind.data.frame(logitdata_f[which(logitdata_f$y==1),],
  #                         logitdata_f[sampids[,i],])
  
  
  # Splitting the data
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #GBM
  gbmModel <- train(y ~ ., data=train, 
                    method = "gbm", 
                    metric="ROC", 
                    trControl = fitControl, 
                    verbose=FALSE, 
                    tuneLength=5)
  pred.gbmModel <- as.vector(predict(gbmModel, 
                                     newdata=test, 
                                     type="prob")[,"Yes"])
  gbmlst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,1]
  gbmlst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,2]
  gbmlst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.gbmModel))
  gbmlst[[4]][,i] <- varImp(gbmModel)$importance[,1]
  
}

mean(gbmlst[[3]][i])

saveRDS(gbmlst, "gbmlst.rds")

###################################################################################

#Without log transform

cols <- c(grep("dist",colnames(chr1data_f)))
chr1data_f[,cols] <- apply(chr1data_f[,cols], 2, function(x){2^x})

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
gbmlst_nl <- list(tpr <- matrix(nrow=ceiling((length(which(chr1data_f$y=="Yes"))*2)*.3), 
                             ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1data_f$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1data_f)[2]-1,
                                ncol=bootsamps))
rownames(gbmlst_nl[[4]]) <- colnames(chr1data_f)[-1]

#gbmlst <- list(tpr <- matrix(nrow=ceiling((length(which(logitdata_f$y==1))*2)*.3), 
#                              ncol=bootsamps),
#                fpr <- matrix(nrow=ceiling((length(which(logitdata_f$y==1))*2)*.3), 
#                              ncol=bootsamps),
#                auc <- numeric(bootsamps),
#                varimp <- matrix(nrow=dim(logitdata_f)[2]-1,
#                                 ncol=bootsamps))
#rownames(gbmlst[[4]]) <- colnames(logitdata_f)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1data_f[which(chr1data_f$y=="Yes"),],
                           chr1data_f[sampids[,i],])
  
  #data <- rbind.data.frame(logitdata_f[which(logitdata_f$y==1),],
  #                         logitdata_f[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #GBM
  gbmModel <- train(y ~ ., data=train, 
                    method = "gbm", 
                    metric="ROC", 
                    trControl = fitControl, 
                    verbose=FALSE, 
                    tuneLength=5)
  pred.gbmModel <- as.vector(predict(gbmModel, 
                                     newdata=test, 
                                     type="prob")[,"Yes"])
  gbmlst_nl[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,1]
  gbmlst_nl[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,2]
  gbmlst_nl[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.gbmModel))
  gbmlst_nl[[4]][,i] <- varImp(gbmModel)$importance[,1]
  
}

mean(gbmlst_nl[[3]])

saveRDS(gbmlst_nl, "gbmlst_nl.rds")