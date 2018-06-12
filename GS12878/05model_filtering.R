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

gs12878 <- readRDS("gs12878.rds")

##################################################################

#Full Data

#Remove CHR variable
gs12878_f <- gs12878[,-which(colnames(gs12878)=="CHR")]

#Taking log2 transform of continous data
cols <- c(grep("dist",colnames(gs12878_f)))
gs12878_f[,cols] <- apply(gs12878_f[,cols], 2, function(x){log(x + 1, base=2)})

#Changing binary variables to factors
cols <- c(grep("dist",colnames(gs12878_f), invert = TRUE))
gs12878_f[,cols] <- lapply(gs12878_f[,cols], factor)

#Changing levels of response (y) to yes no
levels(gs12878_f$y) <- c("No", "Yes")

#Removing zero variance predictors
nzv <- nearZeroVar(gs12878_f[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])
gs12878_f <- gs12878_f[, -which(colnames(gs12878_f) %in% nzvar)]


saveRDS(gs12878_f, "gs12878_f.rds")

##

#Chromosome 1

chr1_gs12878 <- gs12878[which(gs12878$CHR=="chr1"),]

cols <- c(grep("dist",colnames(chr1_gs12878)))
chr1_gs12878[,cols] <- apply(chr1_gs12878[,cols], 2, function(x){log(x + 1, base=2)})

nzv <- nearZeroVar(chr1_gs12878[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])

#Removing zero variance predictors
chr1_gs12878_f <- chr1_gs12878[, -which(colnames(chr1_gs12878) %in% nzvar)]


saveRDS(chr1_gs12878_f, "chr1_gs12878_f.rds")


