#Random Forest

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

#Chromosome 1
chr1data_f <- readRDS("chr1data_f.rds")

#Full data
#logitdata_f <- readRDS("logitdata_f.rds")


#set number of bootstrap samples
bootsamps = 3

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

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1data_f$y==0),
                        length(which(chr1data_f$y==1)),
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

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1data_f$y==1))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1data_f$y==1))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1data_f)[2]-1,
                                 ncol=bootsamps))
rownames(rflst[[4]]) <- colnames(chr1data_f)[-1]

#rflst <- list(tpr <- matrix(nrow=ceiling((length(which(logitdata_f$y==1))*2)*.3), 
#                              ncol=bootsamps),
#                fpr <- matrix(nrow=ceiling((length(which(logitdata_f$y==1))*2)*.3), 
#                              ncol=bootsamps),
#                auc <- numeric(bootsamps),
#                varimp <- matrix(nrow=dim(logitdata_f)[2]-1,
#                                 ncol=bootsamps))
#rownames(rflst[[4]]) <- colnames(logitdata_f)[-1]



for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1data_f[which(chr1data_f$y==1),],
                           chr1data_f[sampids[,i],])
  
  #data <- rbind.data.frame(logitdata_f[which(logitdata_f$y==1),],
  #                         logitdata_f[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #turning response into factor variable "no"/"yes"
  train$y <- as.factor(train$y)
  test$y <- as.factor(test$y)
  levels(train$y) <- c("No", "Yes")
  levels(test$y) <- c("No", "Yes")
  
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
  rflst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.rfModel))
  rflst[[4]][,i] <- varImp(rfModel)$importance[,1]
  
}


saveRDS(rflst, "rflst")
  
