#Evaluating Variable Reduction Techniques using recursive feature selection

#loading packages 

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)
library(leaps)

#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/variable_selection/R.S./")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f <- rbind.data.frame(chr1_gm12878_f[samps,],
                                   chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])


#setting rfe parameters
control <- rfeControl(functions=rfFuncs, method="cv", number=10, repeats=5)

rfeModel <- rfe(chr1_gm12878_f[,-1], 
                chr1_gm12878_f[,1], 
                sizes=c(2:25, 30, 35, 40, 45, 50, 55, 60, 65), 
                metric="Accuracy",
                rfeControl=control)

saveRDS(rfeModel, "rfeModel.rds")