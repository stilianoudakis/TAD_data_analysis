#Evaluating Variable Reduction Techniques using Boruta


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
library(rFerns)

#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/Boruta/")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")


#using BORUTA package

set.seed(123)
boruta_chr1 <- Boruta(y ~ ., data=chr1_gm12878_f,
                      doTrace=2,
                      getImp=getImpFerns)
print(boruta_chr1)

saveRDS(boruta_chr1, "boruta_chr1.rds")

feats <- getSelectedAttributes(boruta_chr1, withTentative = T)

boruta_chr1_gm12878 <- chr1_gm12878_f[,c("y",feats)]

saveRDS(boruta_chr1_gm12878, "boruta_chr1_gm12878.rds")

dim(boruta_chr1_gm12878)
#247632     21

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

boruta_chr1 <- readRDS("boruta_chr1.rds")

par(mar=c(7,4,2,2))
plot(boruta_chr1, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta_chr1$ImpHistory),function(i)
  boruta_chr1$ImpHistory[is.finite(boruta_chr1$ImpHistory[,i]),i])
names(lz) <- colnames(boruta_chr1$ImpHistory)
Labels <- sort(sapply(lz,median))
Labels <- names(Labels)
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
Labels[which(Labels=="novel_sequence_insertion_dist")] <- "NSI_dist"
Labels[which(Labels=="sequence_alteration_dist")] <- "SA_dist"
Labels[which(Labels=="tandem_duplication_dist")] <- "TD_dist"
Labels[which(Labels=="mobile_element_insertion_dist")] <- "MEI_dist"
axis(side = 1,las=2,labels = Labels,
       at = 1:ncol(boruta_chr1$ImpHistory), cex.axis = 0.7)
par(mar=c(4,4,2,2))
