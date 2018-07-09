#Model building for k562 cell line

library(caret)
library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
library(DMwR)
library(gridExtra)
library(ggplot2)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

#classic glm

#full data
k562 <- readRDS("k562.rds")

#chromosome 1
chr1_k562 <- k562[which(k562$CHR=="chr1"),]

#Taking log2 transform of continous data
cols <- c(grep("dist",colnames(chr1_k562)))
chr1_k562[,cols] <- apply(chr1_k562[,cols], 2, function(x){log(x + 1, base=2)})

#Changing binary variables to factors
cols <- c(intersect(grep("score",colnames(chr1_k562), invert = TRUE),
                    grep("dist",colnames(chr1_k562), invert = TRUE)))
chr1_k562[,cols] <- lapply(chr1_k562[,cols], factor)

#Changing levels of response (y) to yes no
levels(chr1_k562$y) <- c("No", "Yes")

#removing factor variables with only 1 factor
chr1_k562 <- chr1_k562[,-which(names(chr1_k562)=="tRNA" | names(chr1_k562)=="CHR")]

#Classic GLM without balanced classes
set.seed(3432)
inTrainingSetglm <- createDataPartition(chr1_k562$y,p=.7,list=FALSE)
trainglm <- chr1_k562[inTrainingSetglm,]
testglm <- chr1_k562[-inTrainingSetglm,]


glmModel <- glm(y ~ ., 
                data = trainglm, 
                family = binomial)
pred.glmModel <- predict(glmModel, newdata=testglm, type="response")
roc.glmModel <- roc(testglm$y, pred.glmModel)
auc.glmModel <- pROC::auc(roc.glmModel)
#Area under the curve: 0.7724

plot(roc.glmModel)


#######################################################################

#Random Forest 

#Chromosome 1
chr1_k562_fwd <- readRDS("chr1_k562_fwd.rds")

#Full data
#logitdata_f <- readRDS("logitdata_f.rds")


#set number of bootstrap samples
bootsamps = 50

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
                  nrow=length(chr1_k562_fwd$y[which(chr1_k562_fwd$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_k562_fwd$y=="No"),
                        length(which(chr1_k562_fwd$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_k562_fwd$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              fpr <- matrix(nrow=ceiling((length(which(chr1_k562_fwd$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              auc <- numeric(bootsamps),
              varimp <- matrix(nrow=dim(chr1_k562_fwd)[2]-1,
                               ncol=bootsamps))
rownames(rflst[[4]]) <- colnames(chr1_k562_fwd)[-1]



for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_k562_fwd[which(chr1_k562_fwd$y=="Yes"),],
                           chr1_k562_fwd[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  
  #Random Forest
  rfModel <- train(y ~ ., data=train, 
                   method = "rf", 
                   metric="ROC", 
                   trControl = fitControl, 
                   verbose=FALSE, 
                   tuneLength=5)
  pred.rfModel <- as.vector(predict(rfModel, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  rflst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.rfModel))
  rflst[[4]][,i] <- varImp(rfModel)$importance[,1]
  
}

k562rflst <- rflst

saveRDS(k562rflst, "k562rflst.rds")



########################################################################

#GBM

#set number of bootstrap samples
bootsamps = 50

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
                  nrow=length(chr1_k562_fwd$y[which(chr1_k562_fwd$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_k562_fwd$y=="No"),
                        length(which(chr1_k562_fwd$y=="Yes")),
                        replace = TRUE)
}



#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
gbmlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_k562_fwd$y=="Yes"))*2)*.3), 
                             ncol=bootsamps),
               fpr <- matrix(nrow=ceiling((length(which(chr1_k562_fwd$y=="Yes"))*2)*.3), 
                             ncol=bootsamps),
               auc <- numeric(bootsamps),
               varimp <- matrix(nrow=dim(chr1_k562_fwd)[2]-1,
                                ncol=bootsamps))
rownames(gbmlst[[4]]) <- colnames(chr1_k562_fwd)[-1]



#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_k562_fwd$y=="No"),
                        length(which(chr1_k562_fwd$y=="Yes")),
                        replace = TRUE)
}


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_k562_fwd[which(chr1_k562_fwd$y=="Yes"),],
                           chr1_k562_fwd[sampids[,i],])
  
  
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

k562gbmlst <- gbmlst

saveRDS(k562gbmlst, "k562gbmlst.rds")



########################################################################

# Elastic Net

#set number of bootstrap samples
bootsamps = 50

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
                  nrow=length(chr1_k562_fwd$y[which(chr1_k562_fwd$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_k562_fwd$y=="No"),
                        length(which(chr1_k562_fwd$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_k562_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1_k562_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1_k562_fwd)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(chr1_k562_fwd)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_k562_fwd[which(chr1_k562_fwd$y=="Yes"),],
                           chr1_k562_fwd[sampids[,i],])
  
  
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
                     standardize=FALSE)
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
}

k562enetlst <- enetlst

saveRDS(k562enetlst, "k562enetlst.rds")