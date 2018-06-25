#Testing SMOTE combinations for Elastic Net

#Using Chromosome 1 data

library(caret)
library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
library(DMwR)
library(gridExtra)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

chr1data_f <- readRDS("chr1data_f.rds")

# Splitting the data
set.seed(5228)
inTrainingSet <- createDataPartition(chr1data_f$y,p=.7,list=FALSE)
train <- chr1data_f[inTrainingSet,]
test <- chr1data_f[-inTrainingSet,]

# Establishing tuning/training parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 3,
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

#Testing eight different combinations of perc.over/perc.under
#100/200, 200/200, 300/200, 400/200, 
#100/300, 200/300, 300/300, 400/300

#########################################################################################

#Elastic Net

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_sm <- list(tpr <- matrix(nrow=length(test$y), 
                              ncol=8),
                fpr <- matrix(nrow=length(test$y), 
                              ncol=8),
                auc <- numeric(8),
                varimp <- matrix(nrow=dim(chr1data_f)[2]-1,
                                 ncol=8))
rownames(enetlst_sm[[4]]) <- colnames(chr1data_f)[-1]


for(i in 1:4){
      set.seed(111)
      train_smote <- SMOTE(y ~ ., 
                         data=train, 
                         perc.over = i*100, 
                         perc.under = 200)
      
      #ENET Model
      eNetModel_sm <- train(y ~ ., data=train_smote, 
                         method = "glmnet", 
                         metric="ROC", 
                         trControl = fitControl, 
                         family="binomial", 
                         tuneLength=5,
                         standardize=TRUE)
      pred.eNetModel <- as.vector(predict(eNetModel_sm, 
                                          newdata=test, 
                                          type="prob")[,"Yes"])
      
      
      enetlst_sm[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
      enetlst_sm[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
      enetlst_sm[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
      enetlst_sm[[4]][,i] <- varImp(eNetModel_sm)$importance[,1]
        
      set.seed(111)
      train_smote <- SMOTE(y ~ ., 
                             data=train, 
                             perc.over = i*100, 
                             perc.under = 300)
        
      #ENET Model
      eNetModel_sm <- train(y ~ ., data=train_smote, 
                              method = "glmnet", 
                              metric="ROC", 
                              trControl = fitControl, 
                              family="binomial", 
                              tuneLength=5,
                              standardize=TRUE)
      pred.eNetModel <- as.vector(predict(eNetModel_sm, 
                                            newdata=test, 
                                            type="prob")[,"Yes"])
        
        
      enetlst_sm[[1]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
      enetlst_sm[[2]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
      enetlst_sm[[3]][i+4] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
      enetlst_sm[[4]][,i+4] <- varImp(eNetModel_sm)$importance[,1]
    
}


#########################################################################################

#Random Forest

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst_sm <- list(tpr <- matrix(nrow=length(test$y), 
                                 ncol=8),
                   fpr <- matrix(nrow=length(test$y), 
                                 ncol=8),
                   auc <- numeric(8),
                   varimp <- matrix(nrow=dim(chr1data_f)[2]-1,
                                    ncol=8))
rownames(rflst_sm[[4]]) <- colnames(chr1data_f)[-1]


for(i in 1:4){
  set.seed(111)
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 200)
  
  #RF Model
  rfModel_sm <- train(y ~ ., data=train_smote, 
                   method = "rf", 
                   metric="ROC", 
                   trControl = fitControl, 
                   verbose=FALSE)
  pred.rfModel <- as.vector(predict(rfModel_sm, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  rflst_sm[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst_sm[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst_sm[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.rfModel))
  rflst_sm[[4]][,i] <- varImp(rfModel_sm)$importance[,1]
  
  set.seed(111)
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 300)
  
  #RF Model
  rfModel_sm <- train(y ~ ., data=train_smote, 
                   method = "rf", 
                   metric="ROC", 
                   trControl = fitControl, 
                   verbose=FALSE)
  pred.rfModel <- as.vector(predict(rfModel_sm, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  rflst_sm[[1]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst_sm[[2]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst_sm[[3]][i+4] <- pROC::auc(pROC::roc(test$y, pred.rfModel))
  rflst_sm[[4]][,i+4] <- varImp(rfModel_sm)$importance[,1]
  
}


#########################################################################################

#GBM Forest

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
gbmlst_sm <- list(tpr <- matrix(nrow=length(test$y), 
                               ncol=8),
                 fpr <- matrix(nrow=length(test$y), 
                               ncol=8),
                 auc <- numeric(8),
                 varimp <- matrix(nrow=dim(chr1data_f)[2]-1,
                                  ncol=8))
rownames(gbmlst_sm[[4]]) <- colnames(chr1data_f)[-1]


for(i in 1:4){
  set.seed(111)
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 200)
  
  #GBM Model
  gbmModel_sm <- train(y ~ ., data=train_smote, 
                      method = "gbm", 
                      metric="ROC", 
                      trControl = fitControl, 
                      verbose=FALSE,
                      tuneLength=5)
  pred.gbmModel <- as.vector(predict(gbmModel_sm, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  gbmlst_sm[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,1]
  gbmlst_sm[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,2]
  gbmlst_sm[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.gbmModel))
  gbmlst_sm[[4]][,i] <- varImp(gbmModel_sm)$importance[,1]
  
  set.seed(111)
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 300)
  
  #GBM Model
  gbmModel_sm <- train(y ~ ., data=train_smote, 
                      method = "gbm", 
                      metric="ROC", 
                      trControl = fitControl,
                      tuneLength=5,
                      verbose=FALSE)
  pred.gbmModel <- as.vector(predict(gbmModel_sm, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  gbmlst_sm[[1]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,1]
  gbmlst_sm[[2]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.gbmModel)[,2]
  gbmlst_sm[[3]][i+4] <- pROC::auc(pROC::roc(test$y, pred.gbmModel))
  gbmlst_sm[[4]][,i+4] <- varImp(gbmModel_sm)$importance[,1]
  
}



#####################################################################

#Plotting performance

#enet
test.auc.enet <- data.frame(Combination=c("100/200","200/200","300/200","400/200",
                                     "100/300","200/300","300/300","400/300"),
                       auc=c(enetlst_sm[[3]][1],enetlst_sm[[3]][2],enetlst_sm[[3]][3],
                             enetlst_sm[[3]][4],enetlst_sm[[3]][5],enetlst_sm[[3]][6],
                             enetlst_sm[[3]][7],enetlst_sm[[3]][8]))

test.auc.enet <- test.auc.enet[order(test.auc.enet$auc, decreasing=TRUE),]

test.auc.enet$Combination <- factor(test.auc.enet$Combination, levels=test.auc.enet$Combination)

test.auc.enet

p<-ggplot(data=test.auc.enet, aes(x=Combination, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") + ylim(0,1)+
  theme_minimal()
p

plot(enetlst_sm[[2]][,1],enetlst_sm[[1]][,1], type="l",col=1)
for(i in 2:8){lines(enetlst_sm[[2]][,i],enetlst_sm[[1]][,i], type="l",col=i)}
legend("bottomright",c("100/200","200/200","300/200",
                       "400/200","100/300","200/300",
                       "300/300","400/300"), fill=1:8)

varimp.enet1 <- as.vector(enetlst_sm[[4]][,1])
varimp.enet.df1 <- data.frame(Feature=rownames(enetlst_sm[[4]]),
                             Importance=varimp.enet1)
varimp.enet.df1 <- varimp.enet.df1[order(varimp.enet.df1$Importance),]
numvarenet <- dim(varimp.enet.df1)[1]
varimp.enet.df1 <- varimp.enet.df1[(numvarenet-19):numvarenet,]
varimp.enet.df1$Feature <- factor(varimp.enet.df1$Feature,levels=varimp.enet.df1$Feature)
enetp1 <- ggplot(varimp.enet.df1, aes(x=Feature, 
                                    y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="green") +
  coord_flip()


#############################################################################

#rf
test.auc.rf <- data.frame(Combination=c("100/200","200/200","300/200","400/200",
                                     "100/300","200/300","300/300","400/300"),
                       auc=c(rflst_sm[[3]][1],rflst_sm[[3]][2],rflst_sm[[3]][3],
                             rflst_sm[[3]][4],rflst_sm[[3]][5],rflst_sm[[3]][6],
                             rflst_sm[[3]][7],rflst_sm[[3]][8]))

test.auc.rf <- test.auc.rf[order(test.auc.rf$auc, decreasing=TRUE),]

test.auc.rf$Combination <- factor(test.auc.rf$Combination, levels=test.auc.rf$Combination)

test.auc.rf

p<-ggplot(data=test.auc.rf, aes(x=Combination, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") + ylim(0,1)+ggtitle("Random Forest")+
  theme_minimal()
p

plot(rflst_sm[[2]][,1],rflst_sm[[1]][,1], type="l",col=1,xlab="1-Specificity", ylab="Sensitivity", main="Elastic Net")
for(i in 2:8){lines(rflst_sm[[2]][,i],rflst_sm[[1]][,i], type="l",col=i)}
legend("bottomright",c("100/200","200/200","300/200",
                       "400/200","100/300","200/300",
                       "300/300","400/300"), fill=1:8)

varimp.rf1 <- as.vector(rflst_sm[[4]][,1])
varimp.rf.df1 <- data.frame(Feature=rownames(rflst_sm[[4]]),
                              Importance=varimp.enet1)
varimp.rf.df1 <- varimp.rf.df1[order(varimp.rf.df1$Importance),]
numvarrf <- dim(varimp.rf.df1)[1]
varimp.rf.df1 <- varimp.rf.df1[(numvarrf-19):numvarrf,]
varimp.rf.df1$Feature <- factor(varimp.rf.df1$Feature,levels=varimp.rf.df1$Feature)
rfp1 <- ggplot(varimp.rf.df1, aes(x=Feature, 
                                      y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()

grid.arrange(rfp1,rfp2,rfp3,rfp4,rfp5,rfp6,rfp7,rfp8,ncol=4,top="Random Forest")

#################################################################################

#gbm
test.auc.gbm <- data.frame(Combination=c("100/200","200/200","300/200","400/200",
                               "100/300","200/300","300/300","400/300"),
                       auc=c(gbmlst_sm[[3]][1],gbmlst_sm[[3]][2],gbmlst_sm[[3]][3],
                             gbmlst_sm[[3]][4],gbmlst_sm[[3]][5],gbmlst_sm[[3]][6],
                             gbmlst_sm[[3]][7],gbmlst_sm[[3]][8]))

test.auc.gbm <- test.auc.gbm[order(test.auc.gbm$auc, decreasing=TRUE),]

test.auc.gbm$Combination <- factor(test.auc.gbm$Combination, levels=test.auc.gbm$Combination)

test.auc.gbm

p<-ggplot(data=test.auc.gbm, aes(x=Combination, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") + ylim(0,1)+ggtitle("GBM")+
  theme_minimal()
p

plot(gbmlst_sm[[2]][,1],gbmlst_sm[[1]][,1], type="l",col=1,xlab="1-Specificity", ylab="Sensitivity", main="GBM")
for(i in 2:8){lines(gbmlst_sm[[2]][,i],gbmlst_sm[[1]][,i], type="l",col=i)}
legend("bottomright",c("100/200","200/200","300/200",
                       "400/200","100/300","200/300",
                       "300/300","400/300"), fill=1:8)

varimp.gbm1 <- as.vector(gbmlst_sm[[4]][,1])
varimp.gbm.df1 <- data.frame(Feature=rownames(gbmlst_sm[[4]]),
                            Importance=varimp.gbm1)
varimp.gbm.df1 <- varimp.gbm.df1[order(varimp.gbm.df1$Importance),]
numvargbm <- dim(varimp.gbm.df1)[1]
varimp.gbm.df1 <- varimp.gbm.df1[(numvargbm-19):numvargbm,]
varimp.gbm.df1$Feature <- factor(varimp.gbm.df1$Feature,levels=varimp.gbm.df1$Feature)
gbmp1 <- ggplot(varimp.gbm.df1, aes(x=Feature, 
                                  y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()

grid.arrange(gbmp1,gbmp2,gbmp3,gbmp4,gbmp5,gbmp6,gbmp7,gbmp8,ncol=4,top="GBM")

