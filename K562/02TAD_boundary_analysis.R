#TAD Boundary Analysis

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

#Reading in TAD data

#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/data")
setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data_analysis/data")

domains <- read.table("arrowhead_k562_data.txt", header=T)
dim(domains)
#5975    3

# Creating TAD boundaries

colnames(domains)[2:3] <- c("coordinate", "coordinate")
coords <- rbind.data.frame(domains[,c(1,2)],domains[,c(1,3)])

#Removing the X chromosome from the analysis
coords <- coords[-which(coords$Chromosome=="chrX"),]

#Sorting the numeric chromosome coordinates
coords <- coords[order(as.numeric(substr(coords$Chromosome,4,5)), coords$coordinate),]

#remove duplicates for coordinates that are conjoined
coords <- coords[!duplicated(coords),]
dim(coords)
#10293     2

# flanking either side of the TAD boundary by 500 bases for a 1kb centered boundary region
coords$Chromosome <- as.character(coords$Chromosome)
bounds <- GRanges(seqnames=coords$Chromosome, ranges=IRanges(start=coords$coordinate, width=1))
bounds <- resize(bounds, 1000, fix = "center")
bounds
prop.table(table(seqnames(bounds)))

# Reading in and cleaning genomic feature data (subcompartments)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
subcompartgr <- readRDS("subcompartgr.rds")

# Finding where the 1kb flanked tad boundaries overlap the subcompartments and labeling the overlap as A or B

#For Within:
within <- findOverlaps(bounds, subcompartgr, type="within")
#there are 16593 intervals from the tad boundaries that overlab strictly within flanked subcompartments
subcompoverlapwithin <- bounds[queryHits(within)]
#attaching the subcompartment info to the within overlaps
mcols(subcompoverlapwithin)$compartment <- mcols(subcompartgr)$compartment[subjectHits(within)]
#adding percent overlap for each A abd B compartments
mcols(subcompoverlapwithin)$percentA <- ifelse(mcols(subcompoverlapwithin)$compartment=="A",1,0)
mcols(subcompoverlapwithin)$percentB <- ifelse(mcols(subcompoverlapwithin)$compartment=="B",1,0)

#For none:
subcompoverlapnone <- bounds[which(bounds %outside% subcompartgr)]
#there are 19 intervals that have no overlaps
#attaching the subcompartment info in the form of "N" for no overlap
mcols(subcompoverlapnone)$compartment <- "N"
#adding percent overlap
mcols(subcompoverlapnone)$percentA <- 0
mcols(subcompoverlapnone)$percentB <- 0

#For partial: the difference among the sets of any type of overlap and within overlaps
subcompoverlapany <- findOverlaps(bounds, subcompartgr, type="any")
partial <- setdiff(subcompoverlapany,within)

# Determing the percentage of partial overlaps

#determining the frequency of partial overlaps
table(queryHits(partial))
table(table(queryHits(partial)))
#1   2 
#2 205 

#First consider where the overlap is with only 1 subcompartment interval
as.numeric(rownames(as.matrix(which(table(queryHits(partial))==1))))
oneinterval <- bounds[as.numeric(rownames(as.matrix(which(table(queryHits(partial))==1))))]
poneinterval <- findOverlapPairs(oneinterval, subcompartgr)
poneinterval <- pintersect(poneinterval)

overlaps1 <- findOverlaps(oneinterval, subcompartgr)
mcols(oneinterval)$compartment <- mcols(subcompartgr)$compartment[subjectHits(overlaps1)]
#attching compartment percent
mcols(oneinterval)$percentA <- ifelse(mcols(oneinterval)$compartment=="A",
                                      round(width(poneinterval)/width(oneinterval),2),
                                      0)
mcols(oneinterval)$percentB <- ifelse(mcols(oneinterval)$compartment=="B",
                                      round(width(poneinterval)/width(oneinterval),2),
                                      0)


#Now we consider the overlaps between two subcompartment intervals
twointerval <- bounds[as.numeric(rownames(as.matrix(which(table(queryHits(partial))==2))))]
ptwointerval <- findOverlapPairs(twointerval, subcompartgr)
ptwointerval <- pintersect(ptwointerval)

overlaps2 <- findOverlaps(twointerval, subcompartgr)
mcols(twointerval)$compartment <- "C"
mcols(twointerval)$percentA <- round(width(ptwointerval[!duplicated(queryHits(overlaps2))])/width(twointerval),2)
mcols(twointerval)$percentB <- round(width(ptwointerval[!duplicated(queryHits(overlaps2))])/width(twointerval),2)

# Concatinating all types of overlaps into 1 granges object

tad_subcomp <- GRangesList(subcompoverlapwithin, subcompoverlapnone, oneinterval, twointerval)
tad_subcomp <- unlist(tad_subcomp)
tad_subcomp <- tad_subcomp[order(substr(seqnames(tad_subcomp),4,5), start(tad_subcomp))]

#attaching a class variable y to denote that a tad boundary exists in the flanked interval
mcols(tad_subcomp)$y <- 1

# Binning the genome

bins <- rep( list(GRangesList()), length(unique(coords$Chromosome)) )

for(i in 1:length(unique(coords$Chromosome))){
  seqn <- unique(coords$Chromosome)[i]
  
  midpt <- ((min(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])-500) + 
              (max(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])+500))/2
  
  bins[[i]] <- GRanges(seqnames=seqn, 
                       ranges=IRanges(start=seq(min(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])-500, 
                                                max(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])+500,
                                                1000),
                                      width=1000))
  
}

binslist <- GRangesList(bins)
binslist <- unlist(binslist)

#Creating an indicator variable Y that denotes whether the tad boundary overlaps the genomic bin 

y <- countOverlaps(binslist, bounds)
length(y)
table(y)
#0       1 
#2616179   10293
prop.table(table(y))
#0           1 
#0.996081055 0.003918945 

# adding the y vector to the bin granges
mcols(binslist)$y <- y

# Finding overlaps between binned genome and subcompartment intervals

#For Within:
within_g <- findOverlaps(binslist, subcompartgr, type="within")
#there are 2612294 intervals from the binned genome that overlab strictly within subcompartment intervals
subcompoverlapwithin_g <- binslist[queryHits(within_g)]
#attaching the subcompartment info to the within overlaps
mcols(subcompoverlapwithin_g)$compartment <- mcols(subcompartgr)$compartment[subjectHits(within_g)]
#adding percent overlap for each A abd B compartments
mcols(subcompoverlapwithin_g)$percentA <- ifelse(mcols(subcompoverlapwithin_g)$compartment=="A",1,0)
mcols(subcompoverlapwithin_g)$percentB <- ifelse(mcols(subcompoverlapwithin_g)$compartment=="B",1,0)

#For none:
subcompoverlapnone_g <- binslist[which(binslist %outside% subcompartgr)]
#there are 150986 intervals that have no overlaps
#attaching the subcompartment info in the form of "N" for no overlap
mcols(subcompoverlapnone_g)$compartment <- "N"
#adding percent overlap
mcols(subcompoverlapnone_g)$percentA <- 0
mcols(subcompoverlapnone_g)$percentB <- 0

#For partial: the difference among the sets of any type of overlap and within overlaps
subcompoverlapany_g <- findOverlaps(binslist, subcompartgr, type="any")
partial_g <- setdiff(subcompoverlapany_g,within_g)


# Determing the percentage of partial overlaps

#determining the frequency of partial overlaps
table(queryHits(partial_g))
table(table(queryHits(partial_g)))

#First consider where the overlap is with only 1 subcompartment interval
oneinterval_g <- binslist[as.numeric(rownames(as.matrix(which(table(queryHits(partial_g))==1))))]
poneinterval_g <- findOverlapPairs(oneinterval_g, subcompartgr)
poneinterval_g <- pintersect(poneinterval_g)

overlaps1_g <- findOverlaps(oneinterval_g, subcompartgr)
mcols(oneinterval_g)$compartment <- mcols(subcompartgr)$compartment[subjectHits(overlaps1_g)]
#attching compartment percent
mcols(oneinterval_g)$percentA <- ifelse(mcols(oneinterval_g)$compartment=="A",
                                        round(width(poneinterval_g)/width(oneinterval_g),2),
                                        0)
mcols(oneinterval_g)$percentB <- ifelse(mcols(oneinterval_g)$compartment=="B",
                                        round(width(poneinterval_g)/width(oneinterval_g),2),
                                        0)


#Now we consider the overlaps between two subcompartment intervals
twointerval_g <- binslist[as.numeric(rownames(as.matrix(which(table(queryHits(partial_g))==2))))]
ptwointerval_g <- findOverlapPairs(twointerval_g, subcompartgr)
ptwointerval_g <- pintersect(ptwointerval_g)

overlaps2_g <- findOverlaps(twointerval_g, subcompartgr)
mcols(twointerval_g)$compartment <- "C"
mcols(twointerval_g)$percentA <- round(width(ptwointerval_g[!duplicated(queryHits(overlaps2_g))])/width(twointerval_g),2)
mcols(twointerval_g)$percentB <- round(width(ptwointerval_g[!duplicated(queryHits(overlaps2_g))])/width(twointerval_g),2)

# Concatinating all types of overlaps into 1 granges object

tad_subcomp_full <- GRangesList(subcompoverlapwithin_g, subcompoverlapnone_g, oneinterval_g, twointerval_g)
tad_subcomp_full <- unlist(tad_subcomp_full)
tad_subcomp_full <- tad_subcomp_full[order(substr(seqnames(tad_subcomp_full),4,5), start(tad_subcomp_full))]

table(mcols(tad_subcomp_full)$y)
#0       1 
#2616179   10293
table(mcols(tad_subcomp_full)$percentA)
#      0     0.5       1 
#1695146    2635  928691
table(mcols(tad_subcomp_full)$percentB)
#      0     0.5       1 
#1052668    2784 1571020
table(mcols(tad_subcomp_full)$compartment)
#      A       B       C       N 
# 928783 1571261    2543  123885 

#saving the types of overlaps

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

saveRDS(tad_subcomp_full, "tad_subcomp_full_k562.rds")
saveRDS(subcompoverlapwithin_g, "within_k562.rds")
saveRDS(subcompoverlapnone_g, "none_k562.rds")
saveRDS(oneinterval_g, "partial1_k562.rds")
saveRDS(twointerval_g, "partial2_k562.rds")

# Creating the dataset for model development

#k562 <- data.frame(y = mcols(tad_subcomp_full)$y,
#                      A = mcols(tad_subcomp_full)$percentA,
#                      B = mcols(tad_subcomp_full)$percentB)

#omit the A and B variables created manually and instead use the subcompartment data from merlot folder
k562 <- data.frame(y = mcols(tad_subcomp_full)$y)


saveRDS(k562, "k562.rds")