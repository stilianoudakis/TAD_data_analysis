#Elastic net

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
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")


#set number of bootstrap samples
bootsamps = 5

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
enetlst_ls <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   auc <- numeric(bootsamps),
                   varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                    ncol=bootsamps))
rownames(enetlst_ls[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_ls <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_ls) <- c("TN",
                           "FN",
                           "FP",
                           "TP",
                           "Total",
                           "Sensitivity",
                           "Specificity",
                           "Kappa",
                           "Accuracy",
                           "Precision",
                           "FPR",
                           "FNR",
                           "FOR",
                           "NPV",
                           "MCC",
                           "F1")

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
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_ls[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_ls[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_ls[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_ls[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_ls[1,i] <- confMat$table[1,1]
  enetperf_ls[2,i] <- confMat$table[1,2]
  enetperf_ls[3,i] <- confMat$table[2,1]
  enetperf_ls[4,i] <- confMat$table[2,2]
  enetperf_ls[5,i] <- sum(confMat$table)
  enetperf_ls[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_ls[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_ls[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_ls[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_ls[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_ls[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_ls[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_ls[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_ls[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_ls[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_ls[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_ls[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_ls[[3]])

saveRDS(enetlst_ls, "enetlst_ls.rds")
saveRDS(enetperf_ls, "enetperf_ls.rds")

#################################################################################

#With log transform and no standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_lns <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    auc <- numeric(bootsamps),
                    varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_lns[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_lns <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_lns) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1")

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
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_lns[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_lns[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_lns[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_lns[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_lns[1,i] <- confMat$table[1,1]
  enetperf_lns[2,i] <- confMat$table[1,2]
  enetperf_lns[3,i] <- confMat$table[2,1]
  enetperf_lns[4,i] <- confMat$table[2,2]
  enetperf_lns[5,i] <- sum(confMat$table)
  enetperf_lns[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_lns[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_lns[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_lns[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_lns[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_lns[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_lns[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_lns[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_lns[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_lns[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_lns[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_lns[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_lns[[3]])

saveRDS(enetlst_lns, "enetlst_lns.rds")
saveRDS(enetperf_lns, "enetperf_lns.rds")



#################################################################################

#Without log transform and standardization

cols <- c(grep("dist",colnames(chr1_gm12878_f)))
chr1_gm12878_f[,cols] <- apply(chr1_gm12878_f[,cols], 2, function(x){2^x})

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_nls <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    auc <- numeric(bootsamps),
                    varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_nls[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_nls <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_nls) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1")

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
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_nls[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_nls[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_nls[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_nls[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_nls[1,i] <- confMat$table[1,1]
  enetperf_nls[2,i] <- confMat$table[1,2]
  enetperf_nls[3,i] <- confMat$table[2,1]
  enetperf_nls[4,i] <- confMat$table[2,2]
  enetperf_nls[5,i] <- sum(confMat$table)
  enetperf_nls[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_nls[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_nls[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_nls[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_nls[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_nls[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_nls[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_nls[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_nls[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_nls[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_nls[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_nls[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_nls[[3]])

saveRDS(enetlst_nls, "enetlst_nls.rds")
saveRDS(enetperf_nls, "enetperf_nls.rds")


####################################################################

#Without log transform and no standardization

enetlst_nlns <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                   ncol=bootsamps),
                     fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                   ncol=bootsamps),
                     auc <- numeric(bootsamps),
                     varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                      ncol=bootsamps))
rownames(enetlst_nlns[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_nlns <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_nlns) <- c("TN",
                             "FN",
                             "FP",
                             "TP",
                             "Total",
                             "Sensitivity",
                             "Specificity",
                             "Kappa",
                             "Accuracy",
                             "Precision",
                             "FPR",
                             "FNR",
                             "FOR",
                             "NPV",
                             "MCC",
                             "F1")

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
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_nlns[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_nlns[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_nlns[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_nlns[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_nlns[1,i] <- confMat$table[1,1]
  enetperf_nlns[2,i] <- confMat$table[1,2]
  enetperf_nlns[3,i] <- confMat$table[2,1]
  enetperf_nlns[4,i] <- confMat$table[2,2]
  enetperf_nlns[5,i] <- sum(confMat$table)
  enetperf_nlns[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_nlns[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_nlns[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_nlns[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_nlns[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_nlns[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_nlns[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_nlns[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_nlns[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_nlns[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_nlns[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_nlns[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_nlns[[3]])

saveRDS(enetlst_nlns, "enetlst_nlns.rds")
saveRDS(enetperf_nlns, "enetperf_nlns.rds")

