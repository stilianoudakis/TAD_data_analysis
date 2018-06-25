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


#using BORUTA package

set.seed(123)
boruta_chr1 <- Boruta(y ~ ., data=chr1_gm12878_f,
                      doTrace=2) #,
                      #getImp=getImpFerns)
print(boruta_chr1)

saveRDS(boruta_chr1, "boruta_chr1.rds")


feats <- getSelectedAttributes(boruta_chr1, withTentative = T)

boruta_chr1_gm12878 <- chr1_gm12878_f[,c("y",feats)]

saveRDS(boruta_chr1_gm12878, "boruta_chr1_gm12878.rds")


