#Evaluating Variable Reduction Techniques using variable selection 

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

#chr1_gm12878_f$A <- as.numeric(chr1_gm12878_f$A)
#chr1_gm12878_f$B <- as.numeric(chr1_gm12878_f$B)

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f <- rbind.data.frame(chr1_gm12878_f[samps,],
                                     chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])

# Performing stepwise selection

#center and scaling data to avoid using intercept term
cols <- names(Filter(is.numeric, chr1_gm12878_f))
chr1_gm12878_f[,cols] <- scale(chr1_gm12878_f[,cols], center = TRUE, scale = TRUE)



#Using cross validation (10 fold)


#backward
k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_gm12878_f), replace=TRUE)
cv.preds.bwd=matrix(NA, nrow=ncol(chr1_gm12878_f)-1,ncol=k)
auc.model.bwd <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_gm12878_f[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_gm12878_f[folds!=j,], family = binomial)
  
  best.fit.bwd = step(glm.full, trace=0)
  
  numpreds <- length(names(best.fit.bwd$coefficients)[-1])
  cv.preds.bwd[(1:numpreds),j] <- names(best.fit.bwd$coefficients)[-1]
  cols <- names(best.fit.bwd$model)
  model <- glm(y ~ . , data = chr1_gm12878_f[folds==j,cols], family = binomial)
  pred.model <- predict(model, newdata=chr1_gm12878_f[folds==j,cols], type="response")
  roc.model <- roc(chr1_gm12878_f[folds==j,"y"], pred.model)
  auc.model.bwd[j] <- pROC::auc(roc.model)
  
}

saveRDS(cv.preds.bwd, "cv.preds.bwd.rds")
saveRDS(auc.model.bwd, "auc.model.bwd.rds")

auc.model.bwd <- readRDS("auc.model.bwd.rds")
cv.preds.bwd <- readRDS("cv.preds.bwd.rds")

vars.bwd <- na.omit(cv.preds.bwd[,which.max(auc.model.bwd)])
vars.bwd[grep("_dist",vars.bwd,invert = TRUE)] <- unlist(lapply(vars.bwd[grep("_dist",vars.bwd,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))

chr1_gm12878_bwd <- chr1_gm12878_f[,which((names(chr1_gm12878_f) %in% vars.bwd) | names(chr1_gm12878_f)=="y")]

dim(chr1_gm12878_bwd)
#247632     23

saveRDS(chr1_gm12878_bwd, "chr1_gm12878_bwd.rds")

################################################################################

#Evaluating performance of reduced data

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
                  nrow=length(chr1_gm12878_bwd$y[which(chr1_gm12878_bwd$y=="Yes")]))


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
enetlst_bwd <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_bwd$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_bwd$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    auc <- numeric(bootsamps),
                    varimp <- matrix(nrow=dim(chr1_gm12878_bwd)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_bwd[[4]]) <- colnames(chr1_gm12878_bwd)[-1]

enetperf_bwd <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_bwd) <- c("TN",
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
  data <- rbind.data.frame(chr1_gm12878_bwd[which(chr1_gm12878_bwd$y=="Yes"),],
                           chr1_gm12878_bwd[sampids[,i],])
  
  
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
  enetlst_bwd[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_bwd[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_bwd[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_bwd[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_bwd[1,i] <- confMat$table[1,1]
  enetperf_bwd[2,i] <- confMat$table[1,2]
  enetperf_bwd[3,i] <- confMat$table[2,1]
  enetperf_bwd[4,i] <- confMat$table[2,2]
  enetperf_bwd[5,i] <- sum(confMat$table)
  enetperf_bwd[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_bwd[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_bwd[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_bwd[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_bwd[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_bwd[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_bwd[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_bwd[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_bwd[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_bwd[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_bwd[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_bwd[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
}

mean(enetlst_bwd[[3]])

saveRDS(enetlst_bwd, "enetlst_bwd.rds")
saveRDS(enetperf_bwd, "enetperf_bwd.rds")