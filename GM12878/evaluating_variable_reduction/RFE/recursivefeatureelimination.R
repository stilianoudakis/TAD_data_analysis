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

#splitting data
set.seed(7215)
inTrainingSet <- sample(length(chr1_gm12878_f$y),floor(length(chr1_gm12878_f$y)*.7))
#inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
train <- chr1_gm12878_f[inTrainingSet,]
test <- chr1_gm12878_f[-inTrainingSet,]

rfeModel <- rfe(train[,-1], 
                train[,1], 
                sizes=c(2:20), 
                metric="Accuracy",
                rfeControl=control,
                trControl = trainctrl)




pred.rfeModel <- as.vector(predict(rfeModel, newdata=test, type="prob")[,"Yes"])

roc.rfeModel <- pROC::roc(test$y, pred.rfeModel)



saveRDS(rfeModel, "rfeModel.rds")
saveRDS(roc.rfeModel, "roc.rfeModel.rds")

