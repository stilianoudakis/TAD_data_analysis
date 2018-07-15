# Variable selection using forward selection for k562 cell line

#loading packages 

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

#Chromosome 1
chr1_k562_f <- readRDS("chr1_k562_f.rds")

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_k562_f$y=="No")
samps <- sample(which(chr1_k562_f$y=="No"),length(which(chr1_k562_f$y=="Yes")))
chr1_k562_f_r <- rbind.data.frame(chr1_k562_f[samps,],
                                chr1_k562_f[which(chr1_k562_f$y=="Yes"),])


# Performing stepwise selection

#center and scaling data to avoid using intercept term
#cols <- names(Filter(is.numeric, chr1_k562_f))
#chr1_k562_f[,cols] <- scale(chr1_k562_f[,cols], center = TRUE, scale = TRUE)



#Using cross validation (10 fold)

#forward
k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_k562_f_r), replace=TRUE)
cv.preds.fwd=matrix(NA, nrow=ncol(chr1_k562_f_r)-1,ncol=k)
auc.model.fwd <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_k562_f_r[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_k562_f_r[folds!=j,], family = binomial)
  
  best.fit.fwd = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="forward",
                      trace=0)
  
  numpreds <- length(names(best.fit.fwd$coefficients)[-1])
  cv.preds.fwd[(1:numpreds),j] <- names(best.fit.fwd$coefficients)[-1]
  cols <- names(best.fit.fwd$model)
  model <- glm(y ~ . , data = chr1_k562_f_r[folds==j,cols], family = binomial)
  pred.model <- predict(model, newdata=chr1_k562_f_r[folds==j,cols], type="response")
  roc.model <- roc(chr1_k562_f_r[folds==j,"y"], pred.model)
  auc.model.fwd[j] <- pROC::auc(roc.model)
  
}

saveRDS(cv.preds.fwd, "cv.preds.fwd.rds")
saveRDS(auc.model.fwd, "auc.model.fwd.rds")

auc.model.fwd <- readRDS("auc.model.fwd.rds")
cv.preds.fwd <- readRDS("cv.preds.fwd.rds")

vars.fwd <- na.omit(cv.preds.fwd[,which.max(auc.model.fwd)])
vars.fwd[grep("_dist",vars.fwd,invert = TRUE)] <- unlist(lapply(vars.fwd[grep("_dist",vars.fwd,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))

chr1_k562_fwd <- chr1_k562_f[,which((names(chr1_k562_f) %in% vars.fwd) | names(chr1_k562_f)=="y")]

dim(chr1_k562_fwd)
#247672     15

saveRDS(chr1_k562_fwd, "chr1_k562_fwd.rds")

#############################################################################

#Backward

#Using cross validation (10 fold)

k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_k562_f_r), replace=TRUE)
cv.preds.bwd=matrix(NA, nrow=ncol(chr1_k562_f_r)-1,ncol=k)
auc.model.bwd <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_k562_f_r[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_k562_f_r[folds!=j,], family = binomial)
  
  best.fit.bwd = step(glm.full, trace=0)
  
  numpreds <- length(names(best.fit.bwd$coefficients)[-1])
  cv.preds.bwd[(1:numpreds),j] <- names(best.fit.bwd$coefficients)[-1]
  cols <- names(best.fit.bwd$model)
  model <- glm(y ~ . , data = chr1_k562_f_r[folds==j,cols], family = binomial)
  pred.model <- predict(model, newdata=chr1_k562_f_r[folds==j,cols], type="response")
  roc.model <- roc(chr1_k562_f_r[folds==j,"y"], pred.model)
  auc.model.bwd[j] <- pROC::auc(roc.model)
  
}


saveRDS(cv.preds.bwd, "cv.preds.fwd.rds")
saveRDS(auc.model.bwd, "auc.model.fwd.rds")

auc.model.bwd <- readRDS("auc.model.bwd.rds")
cv.preds.bwd <- readRDS("cv.preds.bwd.rds")

vars.bwd <- na.omit(cv.preds.bwd[,which.max(auc.model.bwd)])
vars.bwd[grep("_dist",vars.bwd,invert = TRUE)] <- unlist(lapply(vars.bwd[grep("_dist",vars.bwd,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))

chr1_k562_bwd <- chr1_k562_f[,which((names(chr1_k562_f) %in% vars.bwd) | names(chr1_k562_f)=="y")]

dim(chr1_k562_bwd)
#247672     20

saveRDS(chr1_k562_bwd, "chr1_k562_bwd.rds")


############################################################################

#Boruta

boruta_chr1_r <- Boruta(y ~ ., data=chr1_k562_f_r,
                        doTrace=2) #,
#getImp=getImpFerns)
print(boruta_chr1_r)

saveRDS(boruta_chr1_r, "boruta_chr1_r.rds")


feats_r <- getSelectedAttributes(boruta_chr1_r, withTentative = T)

boruta_chr1_k562_r <- chr1_k562_f[,c("y",feats_r)]

dim(boruta_chr1_k562_r)
#247672     39

saveRDS(boruta_chr1_k562_r, "boruta_chr1_k562_r.rds")


setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

boruta_chr1_r <- readRDS("boruta_chr1_r.rds")

par(mar=c(7,4,2,2))
plot(boruta_chr1_r, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta_chr1_r$ImpHistory),function(i)
  boruta_chr1_r$ImpHistory[is.finite(boruta_chr1_r$ImpHistory[,i]),i])
names(lz) <- colnames(boruta_chr1_r$ImpHistory)
Labels <- sort(sapply(lz,median))
Labels <- names(Labels)
Labels[grep("k562_", Labels)] <- gsub("k562_","",Labels[grep("k562_", Labels)])
Labels[which(Labels=="novel_sequence_insertion_dist")] <- "NSI_dist"
Labels[which(Labels=="sequence_alteration_dist")] <- "SA_dist"
Labels[which(Labels=="tandem_duplication_dist")] <- "TD_dist"
Labels[which(Labels=="mobile_element_insertion_dist")] <- "MEI_dist"
axis(side = 1,las=2,labels = Labels,
     at = 1:ncol(boruta_chr1_r$ImpHistory), cex.axis = 0.7)
par(mar=c(4,4,2,2))

#######################################################################

#comparing variable selection techniques

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

#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


# Forward

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

enet.fwd <- enetlst

# Backward

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(chr1_k562_bwd$y[which(chr1_k562_bwd$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_k562_bwd$y=="No"),
                        length(which(chr1_k562_bwd$y=="Yes")),
                        replace = TRUE)
}


#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_k562_bwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1_k562_bwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1_k562_bwd)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(chr1_k562_bwd)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_k562_bwd[which(chr1_k562_bwd$y=="Yes"),],
                           chr1_k562_bwd[sampids[,i],])
  
  
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

enet.bwd <- enetlst


# Boruta

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(boruta_chr1_k562_r$y[which(boruta_chr1_k562_r$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(boruta_chr1_k562_r$y=="No"),
                        length(which(boruta_chr1_k562_r$y=="Yes")),
                        replace = TRUE)
}


#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(boruta_chr1_k562_r$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(boruta_chr1_k562_r$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(boruta_chr1_k562_r)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(boruta_chr1_k562_r)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(boruta_chr1_k562_r[which(boruta_chr1_k562_r$y=="Yes"),],
                           boruta_chr1_k562_r[sampids[,i],])
  
  
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

enet.b.r <- enetlst



# Full

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(chr1_k562_f$y[which(chr1_k562_f$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_k562_f$y=="No"),
                        length(which(chr1_k562_f$y=="Yes")),
                        replace = TRUE)
}


#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_k562_f$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1_k562_f$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1_k562_f)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(chr1_k562_f)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_k562_f[which(chr1_k562_f$y=="Yes"),],
                           chr1_k562_f[sampids[,i],])
  
  
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

enet.full <- enetlst



mean(enet.fwd[[3]])  #0.7686359
mean(enet.bwd[[3]])  #0.7670766
mean(enet.b.r[[3]])  #0.7663008
mean(enet.full[[3]]) #0.7639172


