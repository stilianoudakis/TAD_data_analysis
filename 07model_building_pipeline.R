#Model filtering

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
#library(DT)
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

#Full Data

cols <- c(grep("dist",colnames(logitdata)))
logitdata[,cols] <- apply(logitdata[,cols], 2, function(x){log(x + 1, base=2)})

nzv <- nearZeroVar(logitdata[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])

#Removing zero variance predictors
logitdata_f <- logitdata[, -which(colnames(logitdata) %in% nzvar)]

#check for linear dependencies
comboinfo <- findLinearCombos(logitdata_f)
logitdata_f <- logitdata_f[,-comboinfo$remove]

saveRDS(logitdata_f, "logitdata_f")



#Chromosome 1

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

saveRDS(chr1data_f, "chr1data_f")


