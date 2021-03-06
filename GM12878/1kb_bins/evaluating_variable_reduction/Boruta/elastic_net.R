#Elastic Net using reduced data from Boruta Package

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
library(DMwR)
library(gridExtra)
library(ggplot2)
library(leaps)


setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")


# Full Dataset

boruta_chr1_gm12878 <- readRDS("boruta_chr1_gm12878.rds")


#set number of bootstrap samples
bootsamps = 50

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 3,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(boruta_chr1_gm12878$y[which(boruta_chr1_gm12878$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(boruta_chr1_gm12878$y=="No"),
                        length(which(boruta_chr1_gm12878$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(boruta_chr1_gm12878)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(boruta_chr1_gm12878)[-1]

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(boruta_chr1_gm12878[which(boruta_chr1_gm12878$y=="Yes"),],
                           boruta_chr1_gm12878[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=TRUE)
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
}

enetlst.b <- enetlst


##################################################################

# Reduced Dataset


boruta_chr1_gm12878_r_r <- readRDS("boruta_chr1_gm12878_r.rds")


#set number of bootstrap samples
bootsamps = 50

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 3,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(boruta_chr1_gm12878_r$y[which(boruta_chr1_gm12878_r$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(boruta_chr1_gm12878_r$y=="No"),
                        length(which(boruta_chr1_gm12878_r$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878_r$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878_r$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(boruta_chr1_gm12878_r)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(boruta_chr1_gm12878_r)[-1]

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(boruta_chr1_gm12878_r[which(boruta_chr1_gm12878_r$y=="Yes"),],
                           boruta_chr1_gm12878_r[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=TRUE)
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
}

enetlst.b.r <- enetlst

