#Elastic Net with balanced classes using bootstrap sampling 

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)

setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")

#set number of bootstrap samples
bootsamps = 100

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
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
                  nrow=length(chr1_gm12878_f$y[which(chr1_gm12878_f$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_gm12878_f$y=="No"),
                        length(which(chr1_gm12878_f$y=="Yes")),
                        replace = TRUE)
}



#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

########################################################################################

#With log transform and standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   auc <- numeric(bootsamps),
                   varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                    ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(chr1_gm12878_f)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),],
                           chr1_gm12878_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
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

mean(enetlst[[3]])

saveRDS(enetlst, "enetlst_bs.rds")

