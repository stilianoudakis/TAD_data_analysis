#Elastic Net using Boruta data

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")


#From full data
boruta_chr1_gm12878 <- readRDS("boruta_chr1_gm12878.rds")

boruta_chr1_gm12878$A <- as.numeric(boruta_chr1_gm12878$A)
boruta_chr1_gm12878$B <- as.numeric(boruta_chr1_gm12878$B)

#set number of bootstrap samples
bootsamps = 50

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
                  nrow=length(boruta_chr1_gm12878$y[which(boruta_chr1_gm12878$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(boruta_chr1_gm12878$y=="No"),
                        length(which(boruta_chr1_gm12878$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(boruta_chr1_gm12878)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(boruta_chr1_gm12878)[-1]


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(boruta_chr1_gm12878[which(boruta_chr1_gm12878$y=="Yes"),],
                           boruta_chr1_gm12878[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
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

enetlst <- readRDS("enetlst.boruta.rds")

#Model performance

auc <- mean(enetlst[[3]])
auc
#0.8036025

plot(rowMeans(enetlst[[2]]),rowMeans(enetlst[[1]]), 
     type="l", 
     col="red",
     xlab="1-Specificity",
     ylab="Sensitivity")


varimp.enet <- as.vector(rowMeans(enetlst[[4]]))
varimp.enet.df <- data.frame(Feature=rownames(enetlst[[4]]),
                               Importance=varimp.enet)
varimp.enet.df <- varimp.enet.df[order(varimp.enet.df$Importance),]
numvarenet <- dim(varimp.enet.df)[1]
varimp.enet.df <- varimp.enet.df[(numvarenet-19):numvarenet,]
varimp.enet.df$Feature <- factor(varimp.enet.df$Feature,levels=varimp.enet.df$Feature)
enetp <- ggplot(varimp.enet.df, aes(x=Feature, 
                                       y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()

#########################################################################

#From reduced data
boruta_chr1_gm12878_r <- readRDS("boruta_chr1_gm12878_r.rds")

boruta_chr1_gm12878_r$A <- as.numeric(boruta_chr1_gm12878_r$A)
boruta_chr1_gm12878_r$B <- as.numeric(boruta_chr1_gm12878_r$B)


#set number of bootstrap samples
bootsamps = 50

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
                  nrow=length(boruta_chr1_gm12878_r$y[which(boruta_chr1_gm12878_r$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(boruta_chr1_gm12878_r$y=="No"),
                        length(which(boruta_chr1_gm12878_r$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878_r$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(boruta_chr1_gm12878_r$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(boruta_chr1_gm12878_r)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(boruta_chr1_gm12878_r)[-1]



for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(boruta_chr1_gm12878_r[which(boruta_chr1_gm12878_r$y=="Yes"),],
                           boruta_chr1_gm12878_r[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
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

enetlst_r <- readRDS("enetlst.boruta_r.rds")

#Model performance

auc <- mean(enetlst_r[[3]])
auc
#0.8033327

plot(rowMeans(enetlst_r[[2]]),rowMeans(enetlst_r[[1]]), 
     type="l", 
     col="red",
     xlab="1-Specificity",
     ylab="Sensitivity")


varimp.enet_r <- as.vector(rowMeans(enetlst_r[[4]]))
varimp.enet.df_r <- data.frame(Feature=rownames(enetlst_r[[4]]),
                             Importance=varimp.enet_r)
varimp.enet.df_r <- varimp.enet.df_r[order(varimp.enet.df_r$Importance),]
numvarenet <- dim(varimp.enet.df_r)[1]
varimp.enet.df_r <- varimp.enet.df_r[(numvarenet-19):numvarenet,]
varimp.enet.df_r$Feature <- factor(varimp.enet.df_r$Feature,levels=varimp.enet.df$Feature)
enetp_r <- ggplot(varimp.enet.df, aes(x=Feature, 
                                    y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()


#Comparing Results

x <- intersect(varimp.enet.df$Feature,varimp.enet.df_r$Feature)

varimp.enet.df$ranking <- rank(-varimp.enet.df$Importance)
varimp.enet.df_r$ranking <- rank(-varimp.enet.df_r$Importance)

#rankings between rf, gbm, svm, and elastic net
commonfeatsdf <- data.frame(Features = x,
                            "Full Data" = varimp.enet.df$ranking[varimp.enet.df$Feature %in% x],
                            "FD Importance" = varimp.enet.df$Importance[varimp.enet.df$Feature %in% x],
                            "Reduced Data" = varimp.enet.df_r[order(match(varimp.enet.df_r$Feature, x)),]$ranking[varimp.enet.df_r[order(match(varimp.enet.df_r$Feature, x)),]$Feature %in% x],
                            "FD Importance" = varimp.enet.df_r[order(match(varimp.enet.df_r$Feature, x)),]$Importance[varimp.enet.df_r[order(match(varimp.enet.df_r$Feature, x)),]$Feature %in% x]
)

commonfeatsdf

#saveRDS(commonfeatsdf, "/home/stilianoudakisc/TAD_data_analysis/output/commonfeats_nosmote")

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/common_feats_nosmote")
datatable(commonfeatsdf)
#dev.off()


