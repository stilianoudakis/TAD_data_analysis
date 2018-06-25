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
#library(rFerns)

#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/Boruta/")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")


#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f_r <- rbind.data.frame(chr1_gm12878_f[samps,],
                                     chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])

boruta_chr1_r <- Boruta(y ~ ., data=chr1_gm12878_f_r,
                        doTrace=2) #,
                        #getImp=getImpFerns)
print(boruta_chr1_r)

saveRDS(boruta_chr1_r, "boruta_chr1_r.rds")


feats_r <- getSelectedAttributes(boruta_chr1_r, withTentative = T)

boruta_chr1_gm12878_r <- chr1_gm12878_f[,c("y",feats_r)]

saveRDS(boruta_chr1_gm12878_r, "boruta_chr1_gm12878_r.rds")