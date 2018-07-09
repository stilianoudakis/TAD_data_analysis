#Elastic Net using reduced data from variable selection

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

#Forward Data

chr1_gm12878_fwd <- readRDS("chr1_gm12878_fwd.rds")


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
                  nrow=length(chr1_gm12878_fwd$y[which(chr1_gm12878_fwd$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_gm12878_fwd$y=="No"),
                        length(which(chr1_gm12878_fwd$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1_gm12878_fwd)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(chr1_gm12878_fwd)[-1]

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_fwd[which(chr1_gm12878_fwd$y=="Yes"),],
                           chr1_gm12878_fwd[sampids[,i],])
  
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

enetlst.fwd <- enetlst



########################################################################

#Backward Data

chr1_gm12878_bwd <- readRDS("chr1_gm12878_bwd.rds")


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
                  nrow=length(chr1_gm12878_bwd$y[which(chr1_gm12878_bwd$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_gm12878_bwd$y=="No"),
                        length(which(chr1_gm12878_bwd$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_bwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_bwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1_gm12878_bwd)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(chr1_gm12878_bwd)[-1]

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_bwd[which(chr1_gm12878_bwd$y=="Yes"),],
                           chr1_gm12878_bwd[sampids[,i],])
  
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

enetlst.bwd <- enetlst


############################################################################

auc.fwd <- mean(enetlst.fwd[[3]])
round(auc.fwd,3)

auc.bwd <- mean(enetlst.bwd[[3]])
round(auc.bwd,3)

plot(rowMeans(enetlst.fwd[[2]]),
     rowMeans(enetlst.fwd[[1]]), 
     type="l", col="red",
     xlab="1-Specificity",
     ylab="Sensitivity")
lines(rowMeans(enetlst.bwd[[2]]),
      rowMeans(enetlst.bwd[[1]]),
      type="l", col="blue")
abline(a=0, b=1)

varimp.enet.fwd <- as.vector(rowMeans(enetlst.fwd[[4]]))
varimp.enet.df.fwd <- data.frame(Feature=rownames(enetlst.fwd[[4]]),
                               Importance=varimp.enet.fwd)
varimp.enet.df.fwd <- varimp.enet.df.fwd[order(varimp.enet.df.fwd$Importance),]
varimp.enet.df.fwd$Feature <- factor(varimp.enet.df.fwd$Feature,
                                     levels=varimp.enet.df.fwd$Feature)
enetp.fwd <- ggplot(varimp.enet.df.fwd, aes(x=Feature, 
                                       y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()


varimp.enet.bwd <- as.vector(rowMeans(enetlst.bwd[[4]]))
varimp.enet.df.bwd <- data.frame(Feature=rownames(enetlst.bwd[[4]]),
                                 Importance=varimp.enet.bwd)
varimp.enet.df.bwd <- varimp.enet.df.bwd[order(varimp.enet.df.bwd$Importance),]
varimp.enet.df.bwd$Feature <- factor(varimp.enet.df.bwd$Feature,
                                     levels=varimp.enet.df.bwd$Feature)
enetp.bwd <- ggplot(varimp.enet.df.bwd, aes(x=Feature, 
                                            y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()
