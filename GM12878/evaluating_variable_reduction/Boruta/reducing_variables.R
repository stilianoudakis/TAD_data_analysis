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

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/comparing_normalization/")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")

#randomly sample to reduce dataset
#set.seed(123)
#chr1_gm12878_f <- chr1_gm12878_f[sample(nrow(chr1_gm12878_f),10000),c(1,sample(2:66,20))]


#Full data
#gm12878_f <- readRDS("gm12878_f.rds")

##############################################################

#using BORUTA package

set.seed(123)
boruta_chr1 <- Boruta(y ~ ., data=chr1_gm12878_f,
                      doTrace=2,
                      getImp=getImpFerns)
print(boruta_chr1)

saveRDS(boruta_chr1, "boruta_chr1.rds")

par(mar=c(7, 4, 2, 2))
plot(boruta_chr1, xlab = "", xaxt = "n")

lz<-lapply(1:ncol(boruta_chr1$ImpHistory),function(i)
  boruta_chr1$ImpHistory[is.finite(boruta_chr1$ImpHistory[,i]),i])

names(lz) <- colnames(boruta_chr1$ImpHistory)
Labels <- sort(sapply(lz,median))
Labels <- names(Labels)
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
Labels[which(Labels=="mobile_element_insertion_dist")] <- "MEI_dist"
Labels[which(Labels=="tandem_duplication_dist")] <- "TD_dist"
Labels[which(Labels=="novel_sequence_insertion_dist")] <- "NSI_dist"
Labels[which(Labels=="sequence_alteration_dist")] <- "SA_dist"
axis(side = 1,las=2,labels = Labels,
     at = 1:ncol(boruta_chr1$ImpHistory), cex.axis = 0.7)

par(mar=c(4, 4, 2, 2))

#final.boruta <- TentativeRoughFix(boruta_chr1)
#print(final.boruta)

feats <- getSelectedAttributes(boruta_chr1, withTentative = T)

boruta_chr1_gm12878 <- chr1_gm12878_f[,c("y",feats)]

saveRDS(boruta_chr1_gm12878, "boruta_chr1_gm12878.rds")

##################################################################

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f_r <- rbind.data.frame(chr1_gm12878_f[samps,],
                                     chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])

boruta_chr1_r <- Boruta(y ~ ., data=chr1_gm12878_f_r,
                        doTrace=2,
                        getImp=getImpFerns)
print(boruta_chr1_r)

saveRDS(boruta_chr1_r, "boruta_chr1_r.rds")

par(mar=c(7, 4, 2, 2))
plot(boruta_chr1_r, xlab = "", xaxt = "n")

lz<-lapply(1:ncol(boruta_chr1_r$ImpHistory),function(i)
  boruta_chr1_r$ImpHistory[is.finite(boruta_chr1_r$ImpHistory[,i]),i])

names(lz) <- colnames(boruta_chr1_r$ImpHistory)
Labels <- sort(sapply(lz,median))
Labels <- names(Labels)
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
Labels[which(Labels=="mobile_element_insertion_dist")] <- "MEI_dist"
Labels[which(Labels=="tandem_duplication_dist")] <- "TD_dist"
Labels[which(Labels=="novel_sequence_insertion_dist")] <- "NSI_dist"
Labels[which(Labels=="sequence_alteration_dist")] <- "SA_dist"
axis(side = 1,las=2,labels = Labels,
     at = 1:ncol(boruta_chr1_r$ImpHistory), cex.axis = 0.7)

par(mar=c(4, 4, 2, 2))

#final.boruta <- TentativeRoughFix(boruta_chr1)
#print(final.boruta)

feats_r <- getSelectedAttributes(boruta_chr1_r, withTentative = T)

boruta_chr1_gm12878_r <- chr1_gm12878_f_r[,c("y",feats)]

