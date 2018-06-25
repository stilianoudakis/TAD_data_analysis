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

#Variables that were removed:

[1] "A_dist"                   "B_dist"
[3] "complex"                  "mobile_element_insertion"
[5] "novel_sequence_insertion" "sequence_alteration"
[7] "tandem_duplication"       "deletion_dist"
[9] "gerp_dist"                "gerp_score"
[11] "DNA"                      "low_complexity"
[13] "other"                    "RC"
[15] "RNA"                      "rRNA"
[17] "satellite"                "scRNA"
[19] "simple_repeat"            "SINE"
[21] "snRNA"                    "srpRNA"
[23] "tRNA"                     "unknown"
[25] "line_dist"                "LTR_dist"
[27] "se_k562"                  "UCNE"
[29] "UCNE_score"               "VMR_dist"
[31] "k562_RepetitiveCNV14"     "k562_RepetitiveCNV15"
[33] "k562_ActivePromoter"      "k562_WeakPromoter"
[35] "k562_PoisedPromoter"      "k562_StrongEnhancer4"
[37] "k562_StrongEnhancer5"     "k562_WeakEnhancer6"
[39] "k562_Insulator"           "k562_TxnTransition"
[41] "k562_Heterochromlo_dist"  "k562_CTCF"
[43] "k562_E"                   "k562_PF"
[45] "k562_TSS"                 "k562_WE"
[47] "k562_R_dist"


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

#Variables that were removed:

[1] "A_dist"                   "B_dist"                  
[3] "complex"                  "inversion"               
[5] "mobile_element_insertion" "novel_sequence_insertion"
[7] "sequence_alteration"      "tandem_duplication"      
[9] "gerp_score"               "DNA"                     
[11] "low_complexity"           "other"                   
[13] "RC"                       "RNA"                     
[15] "rRNA"                     "satellite"               
[17] "scRNA"                    "simple_repeat"           
[19] "snRNA"                    "srpRNA"                  
[21] "tRNA"                     "unknown"                 
[23] "se_k562"                  "UCNE"                    
[25] "UCNE_score"               "k562_RepetitiveCNV14"    
[27] "k562_RepetitiveCNV15"     "k562_ActivePromoter"     
[29] "k562_WeakPromoter"        "k562_PoisedPromoter"     
[31] "k562_StrongEnhancer4"     "k562_StrongEnhancer5"    
[33] "k562_WeakEnhancer6"       "k562_Insulator"          
[35] "k562_TxnTransition"       "k562_CTCF"               
[37] "k562_E"                   "k562_PF"                 
[39] "k562_TSS"                 "k562_WE"                 
[41] "CHR"                     

saveRDS(chr1_k562_f, "chr1_k562_f.rds")