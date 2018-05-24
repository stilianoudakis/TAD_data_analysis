#library(MultiAssayExperiment)
#library(GenomicRanges)
#library(IRanges)
library(caret)
library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
library(DMwR)
library(gridExtra)
library(DT)
library(ggplot2)

#setwd("/home/stilianoudakisc/TAD_data_analysis/Rdata")
setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

logitdata <- readRDS("logitdata2.rds")

logitdata <- logitdata[,-which(colnames(logitdata)=="A" 
                               | colnames(logitdata)=="B"
                               | colnames(logitdata)=="UCNE"
                               | colnames(logitdata)=="UCNE_score" 
                               | colnames(logitdata)=="UCNE_dist"
                               | colnames(logitdata)=="gerp_score")]

chr1data <- logitdata[which(logitdata$CHR=="chr1"),]

cols <- c(grep("dist",colnames(chr1data)))
chr1data[,cols] <- apply(chr1data[,cols], 2, function(x){log(x + 1, base=2)})

nzv <- nearZeroVar(chr1data[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])

#Removing zero variance predictors
chr1data_f <- chr1data[, -which(colnames(chr1data) %in% nzvar)]

#check for linear dependencies
comboinfo <- findLinearCombos(chr1data_f)
chr1data_f <- chr1data_f[,-comboinfo$remove]


#Classic GLM without balanced classes
set.seed(3432)
inTrainingSetglm <- createDataPartition(chr1data_f$y,p=.7,list=FALSE)
trainglm <- chr1data_f[inTrainingSetglm,]
testglm <- chr1data_f[-inTrainingSetglm,]

trainglm$y <- as.factor(trainglm$y)
testglm$y <- as.factor(testglm$y)
levels(trainglm$y) <- c("No", "Yes")

glmModel <- glm(y ~ ., 
                data = trainglm, 
                family = binomial)
pred.glmModel <- predict(glmModel, newdata=testglm, type="response")
roc.glmModel <- roc(testglm$y, pred.glmModel)
auc.glmModel <- pROC::auc(roc.glmModel)


############################################################################

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
                  nrow=length(chr1data_f$y[which(chr1data_f$y==1)]))


#set length of list objects that will be filled in with specificities
#and sensitivities
datapart <- ceil(length(chr1data_f$y[which(chr1data_f$y==1)])*2*.3) + 1

auclst <- list(auc.enet <- numeric(bootsamps),
               auc.rf <- numeric(bootsamps),
               auc.gbm <- numeric(bootsamps))
performlst <- list(varimp.enet <- matrix(nrow=dim(chr1data_f)[2]-1,ncol=bootsamps),
                   varimp.rf <- matrix(nrow=dim(chr1data_f)[2]-1,ncol=bootsamps),
                   varimp.gbm <- matrix(nrow=dim(chr1data_f)[2]-1,ncol=bootsamps))
rownames(performlst[[1]]) <- rownames(performlst[[2]]) <- rownames(performlst[[3]]) <- colnames(chr1data_f)[-1]

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1data_f$y==0),
                        length(which(chr1data_f$y==1)),
                        replace = TRUE)
}


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1data_f[which(chr1data_f$y==1),],
                           chr1data_f[sampids[,i],])
  
  #shuffle the data
  #g <- runif(nrow(data))
  #data <- data[order(g),]
  
  # Splitting the data
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #turning response into factor variable "no"/"yes"
  train$y <- as.factor(train$y)
  test$y <- as.factor(test$y)
  levels(train$y) <- c("No", "Yes")
  levels(test$y) <- c("No", "Yes")
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5)
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  roc.eNetModel <- pROC::roc(test$y, pred.eNetModel)
  auclst[[1]][i] <- pROC::auc(roc.eNetModel)
  #enet varimp
  performlst[[1]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Random Forest
  rfModel <- train(y ~ ., data=train, 
                   method = "rf", 
                   metric="ROC", 
                   trControl = fitControl, 
                   verbose=FALSE, 
                   tuneLength=5)
  pred.rfModel <- as.vector(predict(rfModel, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  roc.rfModel <- pROC::roc(test$y, pred.rfModel)
  auclst[[2]][i] <- pROC::auc(roc.rfModel)
  #rf varimp
  performlst[[2]][,i] <- varImp(rfModel)$importance[,1]
  
  #GBM
  gbmModel <- train(y ~ ., data=train, 
                    method = "gbm", 
                    metric="ROC", 
                    trControl = fitControl, 
                    verbose=FALSE, 
                    tuneLength=5)
  pred.gbmModel <- as.vector(predict(gbmModel, 
                                     newdata=test, 
                                     type="prob")[,"Yes"])
  roc.gbmModel <- pROC::roc(test$y, pred.gbmModel)
  auclst[[3]][i] <- pROC::auc(roc.gbmModel)
  #gbm varimp
  performlst[[3]][,i] <- varImp(gbmModel)$importance[,1]
  
  #SVM
  #svmModel <- train(train[,-1], train$y, data = train,
  #                  method='svmRadial',
  #                  trControl=fitControl,
  #                  metric="ROC",
  #                  verbose=FALSE,
  #                  tuneLength = 5)
  #pred.svmModel <- as.vector(predict(svmModel, 
  #                                   newdata=test, 
  #                                   type="prob")[,"Yes"])
  #roc.svmModel <- pROC::roc(test$y, pred.svmModel)
  #auclst[[4]][i] <- pROC::auc(roc.svmModel)
  #svm varimp
  #performlst[[4]][,i] <- varImp(svmModel)$importance[,1]
}



#plotting performance
auc.glmModel <- 0.7884
test.auc <- data.frame(model=c("GLM","ElasticNet","RForest","GBM"),
                       auc=c(auc.glmModel,
                             mean(auclst[[1]]), 
                             mean(auclst[[2]]), 
                             mean(auclst[[3]])))

test.auc <- test.auc[order(test.auc$auc, decreasing=TRUE),]

test.auc$model <- factor(test.auc$model, levels=test.auc$model)

test.auc

perf_nosmote<-ggplot(data=test.auc, aes(x=model, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal()

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/perf_nosmote")
perf_nosmote
#dev.off()


#Variable Importance Plots

#ENET
varimp.enet <- as.vector(rowMeans(performlst[[1]]))
varimp.enet.df <- data.frame(Feature=rownames(performlst[[1]]),
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
           fill="green") +
  coord_flip()
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/enet_varimp_nosmote")
enetp
#dev.off()


#RF
varimp.rf <- as.vector(rowMeans(performlst[[2]]))
varimp.rf.df <- data.frame(Feature=rownames(performlst[[2]]),
                           Importance=varimp.rf)
varimp.rf.df <- varimp.rf.df[order(varimp.rf.df$Importance),]
numvarrf <- dim(varimp.rf.df)[1]
varimp.rf.df <- varimp.rf.df[(numvarrf-19):numvarrf,]
varimp.rf.df$Feature <- factor(varimp.rf.df$Feature,levels=varimp.rf.df$Feature)
rfp <- ggplot(varimp.rf.df, aes(x=Feature, 
                                y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/rf_varimp_nosmote")
rfp
#dev.off()

#GBM
varimp.gbm <- as.vector(rowMeans(performlst[[3]]))
varimp.gbm.df <- data.frame(Feature=rownames(performlst[[3]]),
                            Importance=varimp.gbm)
varimp.gbm.df <- varimp.gbm.df[order(varimp.gbm.df$Importance),]
numvargbm <- dim(varimp.gbm.df)[1]
varimp.gbm.df <- varimp.gbm.df[(numvargbm-19):numvargbm,]
varimp.gbm.df$Feature <- factor(varimp.gbm.df$Feature,levels=varimp.gbm.df$Feature)
gbmp <- ggplot(varimp.gbm.df, aes(x=Feature, 
                                  y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/gbm_varimp_nosmote")
gbmp
#dev.off()

#SVM
#varimp.svm <- as.vector(rowMeans(performlst[[4]]))
#varimp.svm.df <- data.frame(Feature=rownames(performlst[[4]]),
#                            Importance=varimp.svm)
#varimp.svm.df <- varimp.svm.df[order(varimp.svm.df$Importance),]
#numvarsvm <- dim(varimp.svm.df)[1]
#varimp.svm.df <- varimp.svm.df[(numvarsvm-19):numvarsvm,]
#varimp.svm.df$Feature <- factor(varimp.svm.df$Feature,levels=varimp.svm.df$Feature)
#svmp <- ggplot(varimp.svm.df, aes(x=Feature, 
#                                  y=Importance)) +
#  xlab("Predictors") +
#  ylab("Importance") +
#  #ggtitle("Importance Plot for Gradient Boosting Machine") +
#  geom_bar(stat="identity", 
#           width=.5, 
#           position="dodge",
#           fill="red") +
#  coord_flip()
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/varimps_nosmote")
#svmp
#dev.off()

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/svm_varimp_nosmote")
grid.arrange(enetp,rfp,gbmp, ncol=3)
#dev.off()



#Comparing Results
#finding common features between the models

x <- intersect(varimp.rf.df$Feature,varimp.gbm.df$Feature)
z <- intersect(x,varimp.enet.df$Feature)
#z <- intersect(y,varimp.svm.df$Feature)
#z


varimp.rf.df$ranking <- rank(-varimp.rf.df$Importance)
varimp.gbm.df$ranking <- rank(-varimp.gbm.df$Importance)
varimp.enet.df$ranking <- rank(-varimp.enet.df$Importance)
#varimp.svm.df$ranking <- rank(-varimp.svm.df$Importance)


#rankings between rf, gbm, svm, and elastic net
commonfeatsdf <- data.frame(Features = z,
                            RandomForest = varimp.rf.df$ranking[varimp.rf.df$Feature %in% z],
                            RFImp = varimp.rf.df$Importance[varimp.rf.df$Feature %in% z],
                            GBM = varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$ranking[varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$Feature %in% z],
                            GBMImp = varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$Importance[varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$Feature %in% z],
                            ElasticNet = varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$ranking[varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$Feature %in% z],
                            ENetImp = varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$Importance[varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$Feature %in% z]
                            )

commonfeatsdf

#saveRDS(commonfeatsdf, "/home/stilianoudakisc/TAD_data_analysis/output/commonfeats_nosmote")

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/common_feats_nosmote")
datatable(commonfeatsdf)
#dev.off()







