#classic glm


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

#chromosome 1
chr1data_f <- readRDS("chr1data_f.rds")

#full data
#logitdata_f <- readRDS("logitdata_f")

########################################################################################

#With log transform


#Classic GLM without balanced classes
set.seed(3432)
inTrainingSetglm <- createDataPartition(chr1data_f$y,p=.7,list=FALSE)
trainglm <- chr1data_f[inTrainingSetglm,]
testglm <- chr1data_f[-inTrainingSetglm,]

#set.seed(3432)
#inTrainingSetglm <- createDataPartition(logitdata_f$y,p=.7,list=FALSE)
#trainglm <- logitdata_f[inTrainingSetglm,]
#testglm <- logitdata_f[-inTrainingSetglm,]

glmModel <- glm(y ~ ., 
                data = trainglm, 
                family = binomial)
pred.glmModel <- predict(glmModel, newdata=testglm, type="response")
roc.glmModel <- roc(testglm$y, pred.glmModel)
auc.glmModel <- pROC::auc(roc.glmModel)

auc.glmModel

saveRDS(roc.glmModel, "roc.glmModel.rds")


###################################################################################

#Without log transform

cols <- c(grep("dist",colnames(chr1data_f)))
chr1data_f[,cols] <- apply(chr1data_f[,cols], 2, function(x){2^x})

set.seed(3432)
inTrainingSetglm <- createDataPartition(chr1data_f$y,p=.7,list=FALSE)
trainglm <- chr1data_f[inTrainingSetglm,]
testglm <- chr1data_f[-inTrainingSetglm,]

glmModel_nl <- glm(y ~ ., 
                data = trainglm, 
                family = binomial)
pred.glmModel_nl <- predict(glmModel_nl, newdata=testglm, type="response")
roc.glmModel_nl <- roc(testglm$y, pred.glmModel_nl)
auc.glmModel_nl <- pROC::auc(roc.glmModel_nl)

auc.glmModel_nl

saveRDS(roc.glmModel_nl, "roc.glmModel_nl.rds")
