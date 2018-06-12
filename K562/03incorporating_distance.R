#Incorporating distance

# Loading Libraries

#library(MultiAssayExperiment)
library(GenomicRanges)
#library(IRanges)
library(caret)
library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)

# Reading in TAD data that was overlapped within subcompartments

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

tad_subcomp_full_k562 <- readRDS("tad_subcomp_full_k562.rds")
subcompartgr <- readRDS("subcompartgr.rds")
k562 <- readRDS("k562.rds")

# Calculating various distances: Distance to A and Distance to B

Asubcompartgr <- subcompartgr[which(subcompartgr$compartment=="A")]
Bsubcompartgr <- subcompartgr[which(subcompartgr$compartment=="B")]


#A compartment
A_dist <- distanceToNearest(tad_subcomp_full_k562, Asubcompartgr)
A_dist <- mcols(A_dist)$distance

#B compartment
B_dist <- distanceToNearest(tad_subcomp_full_k562, Bsubcompartgr)
B_dist <- mcols(B_dist)$distance

mcols(tad_subcomp_full_k562)$A_dist <- A_dist
mcols(tad_subcomp_full_k562)$B_dist <- B_dist

# Concatinating overlap data and creating update data for modelling

saveRDS(tad_subcomp_full_k562, "tad_subcomp_full_k562.rds")

k562$A_dist <- mcols(tad_subcomp_full_k562)$A_dist
k562$B_dist <- mcols(tad_subcomp_full_k562)$B_dist


saveRDS(k562, "k562.rds")

