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
#gm12878_f$y <- factor(gm12878_f$y,levels(gm12878_f$y)[c(2,1)])

#Removing zero variance predictors
nzv <- nearZeroVar(gm12878_f[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])
gm12878_f <- gm12878_f[, -which(colnames(gm12878_f) %in% nzvar)]

#Variables that were removed:

[1] "A_dist"                     "B_dist"
[3] "complex"                    "mobile_element_insertion"
[5] "novel_sequence_insertion"   "sequence_alteration"
[7] "tandem_duplication"         "deletion_dist"
[9] "gerp_dist"                  "gerp_score"
[11] "DNA"                        "low_complexity"
[13] "other"                      "RC"
[15] "RNA"                        "rRNA"
[17] "satellite"                  "scRNA"
[19] "simple_repeat"              "SINE"
[21] "snRNA"                      "srpRNA"
[23] "tRNA"                       "unknown"
[25] "line_dist"                  "LTR_dist"
[27] "se_GM12878"                 "UCNE"
[29] "UCNE_score"                 "VMR_dist"
[31] "Gm12878_Repressed"          "Gm12878_RepetitiveCNV14"
[33] "Gm12878_RepetitiveCNV15"    "Gm12878_ActivePromoter"
[35] "Gm12878_WeakPromoter"       "Gm12878_PoisedPromoter"
[37] "Gm12878_StrongEnhancer4"    "Gm12878_StrongEnhancer5"
[39] "Gm12878_WeakEnhancer6"      "Gm12878_Insulator"
[41] "Gm12878_TxnTransition"      "Gm12878_Heterochromlo_dist"
[43] "Gm12878_CTCF"               "Gm12878_E"
[45] "Gm12878_PF"                 "Gm12878_TSS"
[47] "Gm12878_WE"                 "Gm12878_R_dist"
[49] "Gm12878_T_dist"


#with centered regions:

[1] "complex"                  "mobile_element_insertion"
[3] "novel_sequence_insertion" "sequence_alteration"     
[5] "tandem_duplication"       "gerp_score"              
[7] "DNA"                      "low_complexity"          
[9] "other"                    "RC"                      
[11] "RNA"                      "rRNA"                    
[13] "satellite"                "scRNA"                   
[15] "simple_repeat"            "SINE"                    
[17] "snRNA"                    "srpRNA"                  
[19] "tRNA"                     "unknown"                 
[21] "se_GM12878"               "UCNE"                    
[23] "UCNE_score"               "Gm12878_Repressed"       
[25] "Gm12878_RepetitiveCNV14"  "Gm12878_RepetitiveCNV15" 
[27] "Gm12878_ActivePromoter"   "Gm12878_WeakPromoter"    
[29] "Gm12878_PoisedPromoter"   "Gm12878_StrongEnhancer4" 
[31] "Gm12878_StrongEnhancer5"  "Gm12878_WeakEnhancer6"   
[33] "Gm12878_Insulator"        "Gm12878_TxnTransition"   
[35] "Gm12878_CTCF"             "Gm12878_E"               
[37] "Gm12878_PF"               "Gm12878_TSS"             
[39] "Gm12878_WE"  

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
#chr1_gm12878$y <- factor(chr1_gm12878$y,levels(chr1_gm12878$y)[c(2,1)])

#Removing zero variance predictors
nzv <- nearZeroVar(chr1_gm12878[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])
chr1_gm12878_f <- chr1_gm12878[, -which(colnames(chr1_gm12878) %in% nzvar)]

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
[23] "se_GM12878"               "UCNE"                    
[25] "UCNE_score"               "Gm12878_RepetitiveCNV14" 
[27] "Gm12878_RepetitiveCNV15"  "Gm12878_ActivePromoter"  
[29] "Gm12878_WeakPromoter"     "Gm12878_PoisedPromoter"  
[31] "Gm12878_StrongEnhancer4"  "Gm12878_StrongEnhancer5" 
[33] "Gm12878_WeakEnhancer6"    "Gm12878_Insulator"       
[35] "Gm12878_TxnTransition"    "Gm12878_CTCF"            
[37] "Gm12878_E"                "Gm12878_PF"              
[39] "Gm12878_TSS"              "Gm12878_WE"              
[41] "CHR"

# with centered regions

[1] "complex"                  "inversion"               
[3] "mobile_element_insertion" "novel_sequence_insertion"
[5] "sequence_alteration"      "tandem_duplication"      
[7] "gerp_score"               "DNA"                     
[9] "low_complexity"           "other"                   
[11] "RC"                       "RNA"                     
[13] "rRNA"                     "satellite"               
[15] "scRNA"                    "simple_repeat"           
[17] "snRNA"                    "srpRNA"                  
[19] "tRNA"                     "unknown"                 
[21] "se_GM12878"               "UCNE"                    
[23] "UCNE_score"               "Gm12878_RepetitiveCNV14" 
[25] "Gm12878_RepetitiveCNV15"  "Gm12878_ActivePromoter"  
[27] "Gm12878_WeakPromoter"     "Gm12878_PoisedPromoter"  
[29] "Gm12878_StrongEnhancer4"  "Gm12878_StrongEnhancer5" 
[31] "Gm12878_WeakEnhancer6"    "Gm12878_Insulator"       
[33] "Gm12878_TxnTransition"    "Gm12878_CTCF"            
[35] "Gm12878_E"                "Gm12878_PF"              
[37] "Gm12878_TSS"              "Gm12878_WE"              
[39] "CHR"  

saveRDS(chr1_gm12878_f, "chr1_gm12878_f.rds")



