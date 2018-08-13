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
chr1data_f <- readRDS("chr1data_f.rds")

#Full data
#logitdata_f <- readRDS("logitdata_f.rds")

#set number of bootstrap samples
bootsamps = 10

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



#####################################################################


#####################################################################

#plotting performance


test.auc <- data.frame(Normalization=c("Log/Standardization", 
                                       "Log/No Standardization",
                                       "No Log/Standardization",
                                       "No Log/No Standardization"),
                       auc=c(mean(enetlst_ls[[3]]),
                             mean(enetlst_lns[[3]]),
                             mean(enetlst_nls[[3]]),
                             mean(enetlst_nlns[[3]])))

test.auc <- test.auc[order(test.auc$auc, decreasing=TRUE),]

test.auc$Normalization <- factor(test.auc$Normalization, levels=test.auc$Normalization)

test.auc

p<-ggplot(data=test.auc, aes(x=Normalization, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") + ylim(0,1) +
  theme_minimal()
p


plot(rowMeans(enetlst_ls[[2]]),rowMeans(enetlst_ls[[1]]), type="l", col="red",xlab="1-Specificity",ylab="Sensitivity")
lines(rowMeans(enetlst_lns[[2]]),rowMeans(enetlst_lns[[1]]), type="l", col="blue")
lines(rowMeans(enetlst_nls[[2]]),rowMeans(enetlst_nls[[1]]), type="l", col="black")
lines(rowMeans(enetlst_nlns[[2]]),rowMeans(enetlst_nlns[[1]]), type="l", col="green")
legend("bottomright", legend = c("Log/ \n Standardization", 
                              "Log/ \n No Standardization",
                              "No Log/ \n Standardization",
                              "No Log/ \n No Standardization"),
       fill=c("red","blue","black","green"),
       cex=.5)


varimp.enetls <- as.vector(rowMeans(enetlst_ls[[4]]))
varimp.enet.dfls <- data.frame(Feature=rownames(enetlst_ls[[4]]),
                              Importance=varimp.enetls)
varimp.enet.dfls <- varimp.enet.dfls[order(varimp.enet.dfls$Importance),]
numvarenet <- dim(varimp.enet.dfls)[1]
varimp.enet.dfls <- varimp.enet.dfls[(numvarenet-19):numvarenet,]
varimp.enet.dfls$Feature <- factor(varimp.enet.dfls$Feature,levels=varimp.enet.dfls$Feature)
enetp1 <- ggplot(varimp.enet.dfls, aes(x=Feature, 
                                      y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()

varimp.enetnls <- as.vector(rowMeans(enetlst_nls[[4]]))
varimp.enet.dfnls <- data.frame(Feature=rownames(enetlst_nls[[4]]),
                               Importance=varimp.enetnls)
varimp.enet.dfnls <- varimp.enet.dfnls[order(varimp.enet.dfnls$Importance),]
numvarenet <- dim(varimp.enet.dfnls)[1]
varimp.enet.dfnls <- varimp.enet.dfnls[(numvarenet-19):numvarenet,]
varimp.enet.dfnls$Feature <- factor(varimp.enet.dfnls$Feature,levels=varimp.enet.dfnls$Feature)
enetp2 <- ggplot(varimp.enet.dfnls, aes(x=Feature, 
                                       y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="black") +
  coord_flip()

varimp.enetlns <- as.vector(rowMeans(enetlst_lns[[4]]))
varimp.enet.dflns <- data.frame(Feature=rownames(enetlst_lns[[4]]),
                                Importance=varimp.enetlns)
varimp.enet.dflns <- varimp.enet.dflns[order(varimp.enet.dflns$Importance),]
numvarenet <- dim(varimp.enet.dflns)[1]
varimp.enet.dflns <- varimp.enet.dflns[(numvarenet-19):numvarenet,]
varimp.enet.dflns$Feature <- factor(varimp.enet.dflns$Feature,levels=varimp.enet.dflns$Feature)
enetp3 <- ggplot(varimp.enet.dflns, aes(x=Feature, 
                                        y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()

varimp.enetnlns <- as.vector(rowMeans(enetlst_nlns[[4]]))
varimp.enet.dfnlns <- data.frame(Feature=rownames(enetlst_nlns[[4]]),
                                Importance=varimp.enetnlns)
varimp.enet.dfnlns <- varimp.enet.dfnlns[order(varimp.enet.dfnlns$Importance),]
numvarenet <- dim(varimp.enet.dfnlns)[1]
varimp.enet.dfnlns <- varimp.enet.dfnlns[(numvarenet-19):numvarenet,]
varimp.enet.dfnlns$Feature <- factor(varimp.enet.dfnlns$Feature,levels=varimp.enet.dfnlns$Feature)
enetp4 <- ggplot(varimp.enet.dfnlns, aes(x=Feature, 
                                        y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="green") +
  coord_flip()

grid.arrange(enetp1, enetp2, enetp3, enetp4, ncol=2)




