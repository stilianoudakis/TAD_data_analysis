#Evaluating Variable Reduction Techniques using recursive feature selection

#loading packages 

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)
library(leaps)

#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/variable_selection/R.S./")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f <- rbind.data.frame(chr1_gm12878_f[samps,],
                                   chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])


#setting rfe parameters
control <- rfeControl(functions=rfFuncs, method="cv", number=10)#, repeats=5)

trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)

#Setting performance	matrix						
enetperf_rfe <- matrix(nrow = 16, ncol=1)
rownames(enetperf_rfe) <- c("TN",
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

#splitting data
set.seed(7215)
inTrainingSet <- sample(length(chr1_gm12878_f2$y),floor(length(chr1_gm12878_f2$y)*.7))
#inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
train <- chr1_gm12878_f2[inTrainingSet,]
test <- chr1_gm12878_f2[-inTrainingSet,]

rfeModel <- rfe(train[,-1], 
                train[,1], 
                sizes=c(2:50), 
                metric="ROC",
                rfeControl=control,
                trControl = trainctrl)

pred.rfeModel <- as.vector(predict(rfeModel, newdata=test, type="prob")[,"Yes"])

pred.rfeModel2 <- predict(rfeModel,
                          newdata=test,
                          type="raw")

confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
TN = as.numeric(confMat$table[1,1])
FN = as.numeric(confMat$table[1,2])
FP = as.numeric(confMat$table[2,1])
TP = as.numeric(confMat$table[2,2])
enetperf_nlns[1,1] <- confMat$table[1,1]
enetperf_nlns[2,1] <- confMat$table[1,2]
enetperf_nlns[3,1] <- confMat$table[2,1]
enetperf_nlns[4,1] <- confMat$table[2,2]
enetperf_nlns[5,1] <- sum(confMat$table)
enetperf_nlns[6,1] <- as.vector(confMat$byClass["Sensitivity"])
enetperf_nlns[7,1] <- as.vector(confMat$byClass["Specificity"])
enetperf_nlns[8,1] <- as.vector(confMat$overall["Kappa"])
enetperf_nlns[9,1] <- as.vector(confMat$overall["Accuracy"])
enetperf_nlns[10,1] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
enetperf_nlns[11,1] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
enetperf_nlns[12,1] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
enetperf_nlns[13,1] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
enetperf_nlns[14,1] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
#enetperf_nlns[15,1] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
enetperf_nlns[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
enetperf_nlns[16,1] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))									 

roc.rfeModel <- pROC::roc(test$y, pred.rfeModel)
auc.rfeModel <- pROC::auc(pROC::roc(test$y, pred.rfeModel))

saveRDS(rfeModel, "rfeModel.rds")
saveRDS(roc.rfeModel, "roc.rfeModel.rds")
saveRDS(auc.rfeModel, "auc.rfeModel.rds")

saveRDS(enetperf_rfe, "enetperf_rfe.rds")

predictors(rfeModel)

vars <- c("y",predictors(rfeModel))

chr1_gm12878_rfe <- chr1_gm12878_f[,which(names(chr1_gm12878_f) %in% vars)]

saveRDS(chr1_gm12878_rfe, "chr1_gm12878_rfe.rds")

