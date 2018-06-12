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

gm12878 <- readRDS("gm12878.rds")

##################################################################

#Full Data

#Remove CHR variable
gm12878_f <- gm12878[,-which(colnames(gm12878)=="CHR")]

#Taking log2 transform of continous data
cols <- c(grep("dist",colnames(gm12878_f)))
gm12878_f[,cols] <- apply(gm12878_f[,cols], 2, function(x){log(x + 1, base=2)})

#Changing binary variables to factors
cols <- c(intersect(grep("score",colnames(gm12878_f), invert = TRUE),
          grep("dist",colnames(gm12878_f), invert = TRUE)))
gm12878_f[,cols] <- lapply(gm12878_f[,cols], factor)

#Changing levels of response (y) to yes no
levels(gm12878_f$y) <- c("No", "Yes")

#Removing zero variance predictors
nzv <- nearZeroVar(gm12878_f[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])
gm12878_f <- gm12878_f[, -which(colnames(gm12878_f) %in% nzvar)]


saveRDS(gm12878_f, "gm12878_f.rds")

##

#Chromosome 1

chr1_gm12878 <- gm12878[which(gm12878$CHR=="chr1"),]

#Taking log2 transform of continous data
cols <- c(grep("dist",colnames(chr1_gm12878)))
chr1_gm12878[,cols] <- apply(chr1_gm12878[,cols], 2, function(x){log(x + 1, base=2)})

#Changing binary variables to factors
cols <- c(intersect(grep("score",colnames(chr1_gm12878), invert = TRUE),
                    grep("dist",colnames(chr1_gm12878), invert = TRUE)))
chr1_gm12878[,cols] <- lapply(chr1_gm12878[,cols], factor)

#Changing levels of response (y) to yes no
levels(chr1_gm12878$y) <- c("No", "Yes")

#Removing zero variance predictors
nzv <- nearZeroVar(chr1_gm12878[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])
chr1_gm12878_f <- chr1_gm12878[, -which(colnames(chr1_gm12878) %in% nzvar)]


saveRDS(chr1_gm12878_f, "chr1_gm12878_f.rds")


