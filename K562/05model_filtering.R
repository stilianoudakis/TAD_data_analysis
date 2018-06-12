#Model Filtering

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

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

k562 <- readRDS("k562.rds")

##################################################################

#Full Data

#Remove CHR variable
k562_f <- k562[,-which(colnames(k562)=="CHR")]

#Taking log2 transform of continous data
cols <- c(grep("dist",colnames(k562_f)))
k562_f[,cols] <- apply(k562_f[,cols], 2, function(x){log(x + 1, base=2)})

#Changing binary variables to factors
cols <- c(intersect(grep("score",colnames(k562_f), invert = TRUE),
          grep("dist",colnames(k562_f), invert = TRUE)))
k562_f[,cols] <- lapply(k562_f[,cols], factor)

#Changing levels of response (y) to yes no
levels(k562_f$y) <- c("No", "Yes")

#Removing zero variance predictors
nzv <- nearZeroVar(k562_f[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])
k562_f <- k562_f[, -which(colnames(k562_f) %in% nzvar)]


saveRDS(k562_f, "k562_f.rds")

##

#Chromosome 1

chr1_k562 <- k562[which(k562$CHR=="chr1"),]

#Taking log2 transform of continous data
cols <- c(grep("dist",colnames(chr1_k562)))
chr1_k562[,cols] <- apply(chr1_k562[,cols], 2, function(x){log(x + 1, base=2)})

#Changing binary variables to factors
cols <- c(intersect(grep("score",colnames(chr1_k562), invert = TRUE),
                    grep("dist",colnames(chr1_k562), invert = TRUE)))
chr1_k562[,cols] <- lapply(chr1_k562[,cols], factor)

#Changing levels of response (y) to yes no
levels(chr1_k562$y) <- c("No", "Yes")

#Removing zero variance predictors
nzv <- nearZeroVar(chr1_k562[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])
chr1_k562_f <- chr1_k562[, -which(colnames(chr1_k562) %in% nzvar)]


saveRDS(chr1_k562_f, "chr1_k562_f.rds")