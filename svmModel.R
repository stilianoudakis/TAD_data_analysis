
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
#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

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


fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

bootsamps=10

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(chr1data_f$y[which(chr1data_f$y==1)]))

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1data_f$y==0),
                        length(which(chr1data_f$y==1)),
                        replace = TRUE)
}


auclst <- list(auc.svm <- numeric(bootsamps))
performlst <- list(varimp.svm <- matrix(nrow=dim(chr1data_f)[2]-1,ncol=bootsamps))
rownames(performlst[[1]]) <- colnames(chr1data_f)[-1]



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

#set.seed(5354)

grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,
                                     0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                           C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                                 1, 1.5, 2,5))
svmModel <- train(y ~., data = train, 
                  method = "svmRadial",
                  trControl=fitControl,
                  metric="ROC",
                  tuneGrid = grid_radial)
pred.svmModel <- as.vector(predict(svmModel, 
                                   newdata=test, 
                                   type="prob")[,"Yes"])
roc.svmModel <- pROC::roc(test$y, pred.svmModel)
auclst[[1]][i] <- pROC::auc(roc.svmModel)
#svm varimp
performlst[[1]][,i] <- varImp(svmModel)$importance[,1]

}

auclst

varimp.svm <- as.vector(rowMeans(performlst[[1]]))
varimp.svm.df <- data.frame(Feature=rownames(performlst[[1]]),
                            Importance=varimp.svm)
varimp.svm.df <- varimp.svm.df[order(varimp.svm.df$Importance),]
numvarsvm <- dim(varimp.svm.df)[1]
varimp.svm.df <- varimp.svm.df[(numvarsvm-19):numvarsvm,]
varimp.svm.df$Feature <- factor(varimp.svm.df$Feature,levels=varimp.svm.df$Feature)
svmp <- ggplot(varimp.svm.df, aes(x=Feature, 
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

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/varimps_nosmote")
svmp
#dev.off()
