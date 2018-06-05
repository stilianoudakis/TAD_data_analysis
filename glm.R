#classic glm


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

#chromosome 1
chr1data_f <- readRDS("chr1data_f")

#full data
logitdata_f <- readRDS("logitdata_f")

#Classic GLM without balanced classes
set.seed(3432)
inTrainingSetglm <- createDataPartition(chr1data_f$y,p=.7,list=FALSE)
trainglm <- chr1data_f[inTrainingSetglm,]
testglm <- chr1data_f[-inTrainingSetglm,]

set.seed(3432)
inTrainingSetglm <- createDataPartition(logitdata_f$y,p=.7,list=FALSE)
trainglm <- logitdata_f[inTrainingSetglm,]
testglm <- logitdata_f[-inTrainingSetglm,]

trainglm$y <- as.factor(trainglm$y)
testglm$y <- as.factor(testglm$y)
levels(trainglm$y) <- c("No", "Yes")

glmModel <- glm(y ~ ., 
                data = trainglm, 
                family = binomial)
pred.glmModel <- predict(glmModel, newdata=testglm, type="response")
roc.glmModel <- roc(testglm$y, pred.glmModel)
auc.glmModel <- pROC::auc(roc.glmModel)

