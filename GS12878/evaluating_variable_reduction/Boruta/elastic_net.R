#Evaluating Variable Reduction Techniques using Boruta

#Elastic Net

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)
library(Boruta)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/comparing_normalization/")

#Chromosome 1
chr1data_f <- readRDS("chr1data_f.rds")

#Full data
#logitdata_f <- readRDS("logitdata_f.rds")

##############################################################

#using BORUTA package

set.seed(123)
boruta_chr1 <- Boruta(y ~ ., data=chr1data_f,
                      doTrace=2)
print(boruta_chr1)

plot(boruta_chr1, xlab = "", xaxt = "n")

lz<-lapply(1:ncol(boruta_chr1$ImpHistory),function(i)
  boruta_chr1$ImpHistory[is.finite(boruta_chr1$ImpHistory[,i]),i])

names(lz) <- colnames(boruta_chr1$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta_chr1$ImpHistory), cex.axis = 0.7)

final.boruta <- TentativeRoughFix(boruta_chr1)
print(final.boruta)

getSelectedAttributes(final.boruta, withTentative = F)

###############################################################

#set number of bootstrap samples
bootsamps = 10

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
                  nrow=length(chr1data_f$y[which(chr1data_f$y=="Yes")]))

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(logitdata_f$y[which(logitdata_f$y==1)]))

#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1data_f$y=="No"),
                        length(which(chr1data_f$y=="Yes")),
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

