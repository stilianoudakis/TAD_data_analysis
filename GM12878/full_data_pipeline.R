# Full Data Pipeline applied to Random Forest


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
library(DMwR)
library(gridExtra)
library(ROCR)

# Reading in TAD boundary data

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data_analysis/data")

domains <- read.table("arrowhead_data.txt", header=T)

# Creating TAD boundaries

#Removing the X chromosome from the domain data
domains <- domains[-which(domains$Chromosome=="chrX"),]

#Sorting the numeric chromosome coordinates
domains <- domains[order(as.numeric(substr(domains$Chromosome,4,5)), domains$Start),]

#adding distance variable
domains$distance <- domains$End - domains$Start
#log2 transform of distance
domains$logdist <- log(domains$distance, base = 2)

# Creating a GRanges object out of the tads

coords <- domains
colnames(coords)[2:3] <- c("coordinate", "coordinate")
coords <- rbind.data.frame(coords[,c(1,2)],coords[,c(1,3)])


#remove duplicates for coordinates that are conjoined
coords <- coords[!duplicated(coords),]
coords <- coords[order(as.numeric(substr(coords$Chromosome,4,5))),]


# flanking either side of the TAD boundary by 500 bases for a 1kb centered boundary region
coords$Chromosome <- as.character(coords$Chromosome)
bounds <- GRanges(seqnames=coords$Chromosome, ranges=IRanges(start=coords$coordinate, width=1))
#bounds <- resize(bounds, 1000, fix = "center")


# Binning the genome

#function to bin genome
binFunc <- function(chroms, boundpts, kb){
  #create empty list
  grlist <- rep( list(GRangesList()), length(unique(chroms)) )
  
  for(i in 1:length(unique(chroms))){
    
    #creating seqnames object for granges
    seqn <- unique(chroms)[i]
    
    #identifying minimum and maximum values to establish sequence for iranges object
    min.bound <- min(boundpts[which(chroms==unique(chroms)[i])])-(kb/2)
    max.bound <- max(boundpts[which(chroms==unique(chroms)[i])])+(kb/2)
    
    grlist[[i]] <- GRanges(seqnames=seqn,
                           ranges=IRanges(start=seq(min.bound,
                                                    max.bound,
                                                    kb),
                                          width = kb))
  }
  
  grlist <- GRangesList(grlist)
  grlist <- unlist(grlist)
  
}


#5kb bins
binslist5 <- binFunc(chroms = coords$Chromosome, boundpts = coords$coordinate, kb = 5000)
y <- countOverlaps(binslist5, bounds)


# Creating the dataset more model development

gm12878 <- data.frame(y = mcols(binslist5)$y)

#Adding genomic annotations

#creating a granges object from the center of every range in the binned genome granges object
#this will be used to calculate distances from genomic feature to center of genomic bin
binslist1_center <- resize(binslist1, width = 1, fix = "center")


#3D subcompartments

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/3D_subcompartments/")

subcomp_A <- read.table("GSE63525_GM12878_subcompartments_A.BED",header = FALSE, sep="\t")
subcomp_B <- read.table("GSE63525_GM12878_subcompartments_B.BED",header = FALSE, sep="\t")

A_gr <- GRanges(seqnames=subcomp_A$V1,IRanges(start=subcomp_A$V2, end=subcomp_A$V3))
A_gr_center <- resize(A_gr, width = 1, fix = "center")
B_gr <- GRanges(seqnames=subcomp_B$V1,IRanges(start=subcomp_B$V2, end=subcomp_B$V3))
B_gr_center <- resize(B_gr, width = 1, fix = "center")

gm12878$A <- ifelse(countOverlaps(binslist5,A_gr)>=1,1,0)
gm12878$B <- ifelse(countOverlaps(binslist5,B_gr)>=1,1,0)

gm12878$A_dist <- mcols(distanceToNearest(binslist5_center, A_gr_center))$distance
gm12878$B_dist <- mcols(distanceToNearest(binslist5_center, B_gr_center))$distance


#DGV

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DGV/")

temp = list.files()

complex <- read.table(temp[1],header=FALSE,sep="\t")
deletion <- read.table(temp[2],header=FALSE,sep="\t")
duplication <- read.table(temp[3],header=FALSE,sep="\t")
gain_loss <- read.table(temp[4],header=FALSE,sep="\t")
insertion <- read.table(temp[5],header=FALSE,sep="\t")
inversion <- read.table(temp[6],header=FALSE,sep="\t")
mobile_element_insertion <- read.table(temp[7],header=FALSE,sep="\t")
novel_sequence_insertion <- read.table(temp[8],header=FALSE,sep="\t")
sequence_alteration <- read.table(temp[9],header=FALSE,sep="\t")
tandem_duplication <- read.table(temp[10],header=FALSE,sep="\t")

complex_gr <- GRanges(seqnames=complex$V1,IRanges(start=complex$V2,end=complex$V3))
complex_gr_center <- resize(complex_gr, width = 1, fix = "center")
deletion_gr <- GRanges(seqnames=deletion$V1,IRanges(start=deletion$V2,end=deletion$V3))
deletion_gr_center <- resize(deletion_gr, width = 1, fix = "center")
duplication_gr <- GRanges(seqnames=duplication$V1,IRanges(start=duplication$V2,end=duplication$V3))
duplication_gr_center <- resize(duplication_gr, width = 1, fix = "center")
gain_loss_gr <- GRanges(seqnames=gain_loss$V1,IRanges(start=gain_loss$V2,end=gain_loss$V3))
gain_loss_gr_center <- resize(gain_loss_gr, width = 1, fix = "center")
insertion_gr <- GRanges(seqnames=insertion$V1,IRanges(start=insertion$V2,end=insertion$V3))
insertion_gr_center <- resize(insertion_gr, width = 1, fix = "center")
inversion_gr <- GRanges(seqnames=inversion$V1,IRanges(start=inversion$V2,end=inversion$V3))
inversion_gr_center <- resize(inversion_gr, width = 1, fix = "center")
mobile_element_insertion_gr <- GRanges(seqnames=mobile_element_insertion$V1,IRanges(start=mobile_element_insertion$V2,end=mobile_element_insertion$V3))
mobile_element_insertion_gr_center <- resize(mobile_element_insertion_gr, width = 1, fix = "center")
novel_sequence_insertion_gr <- GRanges(seqnames=novel_sequence_insertion$V1,IRanges(start=novel_sequence_insertion$V2,end=novel_sequence_insertion$V3))
novel_sequence_insertion_gr_center <- resize(novel_sequence_insertion_gr, width = 1, fix = "center")
sequence_alteration_gr <- GRanges(seqnames=sequence_alteration$V1,IRanges(start=sequence_alteration$V2,end=sequence_alteration$V3))
sequence_alteration_gr_center <- resize(sequence_alteration_gr, width = 1, fix = "center")
tandem_duplication_gr <- GRanges(seqnames=tandem_duplication$V1,IRanges(start=tandem_duplication$V2,end=tandem_duplication$V3))
tandem_duplication_gr_center <- resize(tandem_duplication_gr, width = 1, fix = "center")

gm12878$complex <- ifelse(countOverlaps(binslist5,complex_gr)>=1,1,0)
gm12878$deletion <- ifelse(countOverlaps(binslist5,deletion_gr)>=1,1,0)
gm12878$duplication <- ifelse(countOverlaps(binslist5,duplication_gr)>=1,1,0)
gm12878$gain_loss <- ifelse(countOverlaps(binslist5,gain_loss_gr)>=1,1,0)
gm12878$insertion <- ifelse(countOverlaps(binslist5,insertion_gr)>=1,1,0)
gm12878$inversion <- ifelse(countOverlaps(binslist5,inversion_gr)>=1,1,0)
gm12878$mobile_element_insertion <- ifelse(countOverlaps(binslist5,mobile_element_insertion_gr)>=1,1,0)
gm12878$novel_sequence_insertion <- ifelse(countOverlaps(binslist5,novel_sequence_insertion_gr)>=1,1,0)
gm12878$sequence_alteration <- ifelse(countOverlaps(binslist5,sequence_alteration_gr)>=1,1,0)
gm12878$tandem_duplication <- ifelse(countOverlaps(binslist5,tandem_duplication_gr)>=1,1,0)

gm12878$complex_dist <- mcols(distanceToNearest(binslist5_center, complex_gr_center))$distance
gm12878$deletion_dist <- mcols(distanceToNearest(binslist5_center, deletion_gr_center))$distance
gm12878$duplication_dist <- mcols(distanceToNearest(binslist5_center, duplication_gr_center))$distance
gm12878$gain_loss_dist <- mcols(distanceToNearest(binslist5_center, gain_loss_gr_center))$distance
gm12878$insertion_dist <- mcols(distanceToNearest(binslist5_center, insertion_gr_center))$distance
gm12878$inversion_dist <- mcols(distanceToNearest(binslist5_center, inversion_gr_center))$distance
gm12878$mobile_element_insertion_dist <- mcols(distanceToNearest(binslist5_center, mobile_element_insertion_gr_center))$distance
gm12878$novel_sequence_insertion_dist <- mcols(distanceToNearest(binslist5_center, novel_sequence_insertion_gr_center))$distance
gm12878$sequence_alteration_dist <- mcols(distanceToNearest(binslist5_center, sequence_alteration_gr_center))$distance
gm12878$tandem_duplication_dist <- mcols(distanceToNearest(binslist5_center, tandem_duplication_gr_center))$distance


# GERP

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/GERP/")

gerp <- read.table("GERP_hg19.BED",header=FALSE,sep="\t")

gerp_gr <- GRanges(seqnames=gerp$V1,IRanges(start=gerp$V2,end=gerp$V3))
mcols(gerp_gr)$score <- gerp$V5
gerp_gr_center <- resize(gerp_gr, width = 1, fix = "center")

gm12878$gerp <- ifelse(countOverlaps(binslist5,gerp_gr)>=1,1,0)

gm12878$gerp_dist <- mcols(distanceToNearest(binslist5_center, gerp_gr_center))$distance

#finding which flanks overlap the gerp file so that we can add a score variable
#all other flanks will have a score of 0
which(gm12878$gerp==1)
gm12878$gerp_score <- 0
gerpoverlap <- findOverlaps(binslist5,gerp_gr)
gerpoverlapdf <- data.frame(queryHits=queryHits(gerpoverlap), score=gerp_gr[subjectHits(gerpoverlap)]$score)
gerpoverlapmean <- aggregate(gerpoverlapdf$score, list(gerpoverlapdf$queryHits), mean)
gm12878$gerp_score[gerpoverlapmean$Group.1] <- gerpoverlapmean$x


# nestedRepeats

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/nestedRepeats/")

temp = list.files()

DNA <- read.table(temp[1],header=FALSE,sep="\t")
line <- read.table(temp[2],header=FALSE,sep="\t")
low_complexity <- read.table(temp[3],header=FALSE,sep="\t")
LTR  <- read.table(temp[4],header=FALSE,sep="\t")
other <- read.table(temp[5],header=FALSE,sep="\t")
RC <- read.table(temp[6],header=FALSE,sep="\t")
RNA <- read.table(temp[7],header=FALSE,sep="\t")
rRNA <- read.table(temp[8],header=FALSE,sep="\t")
satellite <- read.table(temp[9],header=FALSE,sep="\t")
scRNA <- read.table(temp[10],header=FALSE,sep="\t")
simple_repeat <- read.table(temp[11],header=FALSE,sep="\t")
SINE <- read.table(temp[12],header=FALSE,sep="\t")
snRNA <- read.table(temp[13],header=FALSE,sep="\t")
srpRNA <- read.table(temp[14],header=FALSE,sep="\t")
tRNA <- read.table(temp[15],header=FALSE,sep="\t")
unknown <- read.table(temp[16],header=FALSE,sep="\t")

DNA_gr <- GRanges(seqnames=DNA$V1,IRanges(start=DNA$V2,end=DNA$V3))
DNA_gr_center <- resize(DNA_gr, width = 1, fix = "center")
line_gr <- GRanges(seqnames=line$V1,IRanges(start=line$V2,end=line$V3))
line_gr_center <- resize(line_gr, width = 1, fix = "center")
low_complexity_gr <- GRanges(seqnames=low_complexity$V1,IRanges(start=low_complexity$V2,end=low_complexity$V3))
low_complexity_gr_center <- resize(low_complexity_gr, width = 1, fix = "center")
LTR_gr <- GRanges(seqnames=LTR$V1,IRanges(start=LTR$V2,end=LTR$V3))
LTR_gr_center <- resize(LTR_gr, width = 1, fix = "center")
other_gr <- GRanges(seqnames=other$V1,IRanges(start=other$V2,end=other$V3))
other_gr_center <- resize(other_gr, width = 1, fix = "center")
RC_gr <- GRanges(seqnames=RC$V1,IRanges(start=RC$V2,end=RC$V3))
RC_gr_center <- resize(RC_gr, width = 1, fix = "center")
RNA_gr <- GRanges(seqnames=RNA$V1,IRanges(start=RNA$V2,end=RNA$V3))
RNA_gr_center <- resize(RNA_gr, width = 1, fix = "center")
rRNA_gr <- GRanges(seqnames=rRNA$V1,IRanges(start=rRNA$V2,end=rRNA$V3))
rRNA_gr_center <- resize(rRNA_gr, width = 1, fix = "center")
satellite_gr <- GRanges(seqnames=satellite$V1,IRanges(start=satellite$V2,end=satellite$V3))
satellite_gr_center <- resize(satellite_gr, width = 1, fix = "center")
scRNA_gr <- GRanges(seqnames=scRNA$V1,IRanges(start=scRNA$V2,end=scRNA$V3))
scRNA_gr_center <- resize(scRNA_gr, width = 1, fix = "center")
simple_repeat_gr <- GRanges(seqnames=simple_repeat$V1,IRanges(start=simple_repeat$V2,end=simple_repeat$V3))
simple_repeat_gr_center <- resize(simple_repeat_gr, width = 1, fix = "center")
SINE_gr <- GRanges(seqnames=SINE$V1,IRanges(start=SINE$V2,end=SINE$V3))
SINE_gr_center <- resize(SINE_gr, width = 1, fix = "center")
snRNA_gr <- GRanges(seqnames=snRNA$V1,IRanges(start=snRNA$V2,end=snRNA$V3))
snRNA_gr_center <- resize(snRNA_gr, width = 1, fix = "center")
srpRNA_gr <- GRanges(seqnames=srpRNA$V1,IRanges(start=srpRNA$V2,end=srpRNA$V3))
srpRNA_gr_center <- resize(srpRNA_gr, width = 1, fix = "center")
tRNA_gr <- GRanges(seqnames=tRNA$V1,IRanges(start=tRNA$V2,end=tRNA$V3))
tRNA_gr_center <- resize(tRNA_gr, width = 1, fix = "center")
unknown_gr <- GRanges(seqnames=unknown$V1,IRanges(start=unknown$V2,end=unknown$V3))
unknown_gr_center <- resize(unknown_gr, width = 1, fix = "center")

gm12878$DNA <- ifelse(countOverlaps(binslist5,DNA_gr)>=1,1,0)
gm12878$line <- ifelse(countOverlaps(binslist5,line_gr)>=1,1,0)
gm12878$low_complexity <- ifelse(countOverlaps(binslist5,low_complexity_gr)>=1,1,0)
gm12878$LTR <- ifelse(countOverlaps(binslist5,LTR_gr)>=1,1,0)
gm12878$other <- ifelse(countOverlaps(binslist5,other_gr)>=1,1,0)
gm12878$RC <- ifelse(countOverlaps(binslist5,RC_gr)>=1,1,0)
gm12878$RNA <- ifelse(countOverlaps(binslist5,RNA_gr)>=1,1,0)
gm12878$rRNA <- ifelse(countOverlaps(binslist5,rRNA_gr)>=1,1,0)
gm12878$satellite <- ifelse(countOverlaps(binslist5,satellite_gr)>=1,1,0)
gm12878$scRNA <- ifelse(countOverlaps(binslist5,scRNA_gr)>=1,1,0)
gm12878$simple_repeat <- ifelse(countOverlaps(binslist5,simple_repeat_gr)>=1,1,0)
gm12878$SINE <- ifelse(countOverlaps(binslist5,SINE_gr)>=1,1,0)
gm12878$snRNA <- ifelse(countOverlaps(binslist5,snRNA_gr)>=1,1,0)
gm12878$srpRNA <- ifelse(countOverlaps(binslist5,srpRNA_gr)>=1,1,0)
gm12878$tRNA <- ifelse(countOverlaps(binslist5,tRNA_gr)>=1,1,0)
gm12878$unknown <- ifelse(countOverlaps(binslist5,unknown_gr)>=1,1,0)

gm12878$DNA_dist <- mcols(distanceToNearest(binslist5_center, DNA_gr_center))$distance 
gm12878$line_dist <- mcols(distanceToNearest(binslist5_center, line_gr_center))$distance 
gm12878$low_complexity_dist <- mcols(distanceToNearest(binslist5_center, low_complexity_gr_center))$distance 
gm12878$LTR_dist <- mcols(distanceToNearest(binslist5_center, LTR_gr_center))$distance 
gm12878$other_dist <- mcols(distanceToNearest(binslist5_center, other_gr_center))$distance 
gm12878$RC_dist <- mcols(distanceToNearest(binslist5_center, RC_gr_center))$distance 
#gm12878$RNA_dist <- mcols(distanceToNearest(binslist5_center, RNA_gr_center))$distance 
#gm12878$rRNA_dist <- mcols(distanceToNearest(binslist5_center, rRNA_gr_center))$distance 
gm12878$satellite_dist <- mcols(distanceToNearest(binslist5_center, satellite_gr_center))$distance 
#gm12878$scRNA_dist <- mcols(distanceToNearest(binslist5_center, scRNA_gr_center))$distance 
gm12878$simple_repeat_dist <- mcols(distanceToNearest(binslist5_center, simple_repeat_gr_center))$distance 
gm12878$SINE_dist <- mcols(distanceToNearest(binslist5_center, SINE_gr_center))$distance 
#gm12878$snRNA_dist <- mcols(distanceToNearest(binslist5_center, snRNA_gr_center))$distance 
gm12878$srpRNA_dist <- mcols(distanceToNearest(binslist5_center, srpRNA_gr_center))$distance 
#gm12878$tRNA_dist <- mcols(distanceToNearest(binslist5_center, tRNA_gr_center))$distance 
gm12878$unknown_dist <- mcols(distanceToNearest(binslist5_center, unknown_gr_center))$distance 


# super_enhancers

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/super_enhancers/")

temp = list.files()

se_GM12878 <- read.table(temp[1],header=FALSE,sep="\t")

se_GM12878_gr <- GRanges(seqnames=se_GM12878$V1,IRanges(start=se_GM12878$V2,end=se_GM12878$V3))
se_GM12878_gr_center <- resize(se_GM12878_gr, width = 1, fix = "center")

gm12878$se_GM12878 <- ifelse(countOverlaps(binslist5,se_GM12878_gr)>=1,1,0)

gm12878$se_GM12878_dist <- mcols(distanceToNearest(binslist5_center, se_GM12878_gr_center))$distance


# UCNE

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/UCNEs/")

UCNE <- read.table("UCNE.BED",header=FALSE,sep="\t")

UCNE_gr <- GRanges(seqnames=UCNE$V1,IRanges(start=UCNE$V2,end=UCNE$V3))
UCNE_gr_center <- resize(UCNE_gr, width = 1, fix = "center")
mcols(UCNE_gr)$score <- UCNE$V5

gm12878$UCNE <- ifelse(countOverlaps(binslist5,UCNE_gr)>=1,1,0)

gm12878$UCNE_dist <- mcols(distanceToNearest(binslist5_center, UCNE_gr_center))$distance

#finding which flanks overlap the unce file so that we can add a score variable
#all other flanks will have a score of 0
which(gm12878$UCNE==1)
gm12878$UCNE_score <- 0
UCNEoverlap <- findOverlaps(binslist5,UCNE_gr)
UCNEoverlapdf <- data.frame(queryHits=queryHits(UCNEoverlap), score=UCNE_gr[subjectHits(UCNEoverlap)]$score)
UCNEoverlapmean <- aggregate(UCNEoverlapdf$score, list(UCNEoverlapdf$queryHits), mean)
gm12878$UCNE_score[UCNEoverlapmean$Group.1] <- UCNEoverlapmean$x


# VMR

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/VMRs/")

VMR <- read.table("VMR_hg19.BED",header=FALSE,sep="\t")

VMR_gr <- GRanges(seqnames=VMR$V1,IRanges(start=VMR$V2,end=VMR$V3))
VMR_gr_center <- resize(VMR_gr, width = 1, fix = "center")

gm12878$VMR <- ifelse(countOverlaps(binslist5,VMR_gr)>=1,1,0)

gm12878$VMR_dist <- mcols(distanceToNearest(binslist5_center, VMR_gr_center))$distance


# BroadHMM

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/BroadHMM/")

temp <- list.files()

#Gm12878
Gm12878_TxnElongation <- read.table(temp[1],header=FALSE,sep="\t")
Gm12878_WeakTxn <- read.table(temp[2],header=FALSE,sep="\t")
Gm12878_Repressed <- read.table(temp[3],header=FALSE,sep="\t")
Gm12878_Heterochromlo <- read.table(temp[4],header=FALSE,sep="\t")
Gm12878_RepetitiveCNV14 <- read.table(temp[5],header=FALSE,sep="\t")
Gm12878_RepetitiveCNV15 <- read.table(temp[6],header=FALSE,sep="\t")
Gm12878_ActivePromoter <- read.table(temp[7],header=FALSE,sep="\t")
Gm12878_WeakPromoter <- read.table(temp[8],header=FALSE,sep="\t")
Gm12878_PoisedPromoter <- read.table(temp[9],header=FALSE,sep="\t")
Gm12878_StrongEnhancer4 <- read.table(temp[10],header=FALSE,sep="\t")
Gm12878_StrongEnhancer5 <- read.table(temp[11],header=FALSE,sep="\t") 
Gm12878_WeakEnhancer6 <- read.table(temp[12],header=FALSE,sep="\t")
Gm12878_WeakEnhancer7 <- read.table(temp[13],header=FALSE,sep="\t")
Gm12878_Insulator <- read.table(temp[14],header=FALSE,sep="\t")
Gm12878_TxnTransition <- read.table(temp[15],header=FALSE,sep="\t")

Gm12878_TxnElongation_gr <- GRanges(seqnames=Gm12878_TxnElongation$V1,IRanges(start=Gm12878_TxnElongation$V2,end=Gm12878_TxnElongation$V3))
Gm12878_TxnElongation_gr_center <- resize(Gm12878_TxnElongation_gr, width = 1, fix = "center")
Gm12878_WeakTxn_gr <- GRanges(seqnames=Gm12878_WeakTxn$V1,IRanges(start=Gm12878_WeakTxn$V2,end=Gm12878_WeakTxn$V3))
Gm12878_WeakTxn_gr_center <- resize(Gm12878_WeakTxn_gr, width = 1, fix = "center")
Gm12878_Repressed_gr <- GRanges(seqnames=Gm12878_Repressed$V1,IRanges(start=Gm12878_Repressed$V2,end=Gm12878_Repressed$V3))
Gm12878_Repressed_gr_center <- resize(Gm12878_Repressed_gr, width = 1, fix = "center")
Gm12878_Heterochromlo_gr <- GRanges(seqnames=Gm12878_Heterochromlo$V1,IRanges(start=Gm12878_Heterochromlo$V2,end=Gm12878_Heterochromlo$V3)) 
Gm12878_Heterochromlo_gr_center <- resize(Gm12878_Heterochromlo_gr, width = 1, fix = "center")
Gm12878_RepetitiveCNV14_gr <- GRanges(seqnames=Gm12878_RepetitiveCNV14$V1,IRanges(start=Gm12878_RepetitiveCNV14$V2,end=Gm12878_RepetitiveCNV14$V3)) 
Gm12878_RepetitiveCNV14_gr_center <- resize(Gm12878_RepetitiveCNV14_gr, width = 1, fix = "center")
Gm12878_RepetitiveCNV15_gr <- GRanges(seqnames=Gm12878_RepetitiveCNV15$V1,IRanges(start=Gm12878_RepetitiveCNV15$V2,end=Gm12878_RepetitiveCNV15$V3))
Gm12878_RepetitiveCNV15_gr_center <- resize(Gm12878_RepetitiveCNV15_gr, width = 1, fix = "center")
Gm12878_ActivePromoter_gr <- GRanges(seqnames=Gm12878_ActivePromoter$V1,IRanges(start=Gm12878_ActivePromoter$V2,end=Gm12878_ActivePromoter$V3))
Gm12878_ActivePromoter_gr_center <- resize(Gm12878_ActivePromoter_gr, width = 1, fix = "center")
Gm12878_WeakPromoter_gr <- GRanges(seqnames=Gm12878_WeakPromoter$V1,IRanges(start=Gm12878_WeakPromoter$V2,end=Gm12878_WeakPromoter$V3))
Gm12878_WeakPromoter_gr_center <- resize(Gm12878_WeakPromoter_gr, width = 1, fix = "center")
Gm12878_PoisedPromoter_gr <- GRanges(seqnames=Gm12878_PoisedPromoter$V1,IRanges(start=Gm12878_PoisedPromoter$V2,end=Gm12878_PoisedPromoter$V3)) 
Gm12878_PoisedPromoter_gr_center <- resize(Gm12878_PoisedPromoter_gr, width = 1, fix = "center")
Gm12878_StrongEnhancer4_gr <- GRanges(seqnames=Gm12878_StrongEnhancer4$V1,IRanges(start=Gm12878_StrongEnhancer4$V2,end=Gm12878_StrongEnhancer4$V3))
Gm12878_StrongEnhancer4_gr_center <- resize(Gm12878_StrongEnhancer4_gr, width = 1, fix = "center") 
Gm12878_StrongEnhancer5_gr <- GRanges(seqnames=Gm12878_StrongEnhancer5$V1,IRanges(start=Gm12878_StrongEnhancer5$V2,end=Gm12878_StrongEnhancer5$V3))
Gm12878_StrongEnhancer5_gr_center <- resize(Gm12878_StrongEnhancer5_gr, width = 1, fix = "center")
Gm12878_WeakEnhancer6_gr <- GRanges(seqnames=Gm12878_WeakEnhancer6$V1,IRanges(start=Gm12878_WeakEnhancer6$V2,end=Gm12878_WeakEnhancer6$V3)) 
Gm12878_WeakEnhancer6_gr_center <- resize(Gm12878_WeakEnhancer6_gr, width = 1, fix = "center")
Gm12878_WeakEnhancer7_gr <- GRanges(seqnames=Gm12878_WeakEnhancer7$V1,IRanges(start=Gm12878_WeakEnhancer7$V2,end=Gm12878_WeakEnhancer7$V3))
Gm12878_WeakEnhancer7_gr_center <- resize(Gm12878_WeakEnhancer7_gr, width = 1, fix = "center") 
Gm12878_Insulator_gr <- GRanges(seqnames=Gm12878_Insulator$V1,IRanges(start=Gm12878_Insulator$V2,end=Gm12878_Insulator$V3))
Gm12878_Insulator_gr_center <- resize(Gm12878_Insulator_gr, width = 1, fix = "center")
Gm12878_TxnTransition_gr <- GRanges(seqnames=Gm12878_TxnTransition$V1,IRanges(start=Gm12878_TxnTransition$V2,end=Gm12878_TxnTransition$V3)) 
Gm12878_TxnTransition_gr_center <- resize(Gm12878_TxnTransition_gr, width = 1, fix = "center")

gm12878$Gm12878_TxnElongation <- ifelse(countOverlaps(binslist5,Gm12878_TxnElongation_gr)>=1,1,0)
gm12878$Gm12878_WeakTxn <- ifelse(countOverlaps(binslist5,Gm12878_WeakTxn_gr)>=1,1,0)
gm12878$Gm12878_Repressed <- ifelse(countOverlaps(binslist5,Gm12878_Repressed_gr)>=1,1,0)
gm12878$Gm12878_Heterochromlo <- ifelse(countOverlaps(binslist5,Gm12878_Heterochromlo_gr)>=1,1,0)
gm12878$Gm12878_RepetitiveCNV14 <- ifelse(countOverlaps(binslist5,Gm12878_RepetitiveCNV14_gr)>=1,1,0)
gm12878$Gm12878_RepetitiveCNV15 <- ifelse(countOverlaps(binslist5,Gm12878_RepetitiveCNV15_gr)>=1,1,0)
gm12878$Gm12878_ActivePromoter <- ifelse(countOverlaps(binslist5,Gm12878_ActivePromoter_gr)>=1,1,0)
gm12878$Gm12878_WeakPromoter <- ifelse(countOverlaps(binslist5,Gm12878_WeakPromoter_gr)>=1,1,0)
gm12878$Gm12878_PoisedPromoter <- ifelse(countOverlaps(binslist5,Gm12878_PoisedPromoter_gr)>=1,1,0)
gm12878$Gm12878_StrongEnhancer4 <- ifelse(countOverlaps(binslist5,Gm12878_StrongEnhancer4_gr)>=1,1,0)
gm12878$Gm12878_StrongEnhancer5 <- ifelse(countOverlaps(binslist5,Gm12878_StrongEnhancer5_gr)>=1,1,0)
gm12878$Gm12878_WeakEnhancer6 <- ifelse(countOverlaps(binslist5,Gm12878_WeakEnhancer6_gr)>=1,1,0)
gm12878$Gm12878_WeakEnhancer7 <- ifelse(countOverlaps(binslist5,Gm12878_WeakEnhancer7_gr)>=1,1,0)
gm12878$Gm12878_Insulator <- ifelse(countOverlaps(binslist5,Gm12878_Insulator_gr)>=1,1,0)
gm12878$Gm12878_TxnTransition <- ifelse(countOverlaps(binslist5,Gm12878_TxnTransition_gr)>=1,1,0)

gm12878$Gm12878_TxnElongation_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_TxnElongation_gr_center))$distance
gm12878$Gm12878_WeakTxn_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakTxn_gr_center))$distance
gm12878$Gm12878_Repressed_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_Repressed_gr_center))$distance
gm12878$Gm12878_Heterochromlo_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_Heterochromlo_gr_center))$distance
gm12878$Gm12878_RepetitiveCNV14_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_RepetitiveCNV14_gr_center))$distance
gm12878$Gm12878_RepetitiveCNV15_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_RepetitiveCNV15_gr_center))$distance
gm12878$Gm12878_ActivePromoter_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_ActivePromoter_gr_center))$distance
gm12878$Gm12878_WeakPromoter_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakPromoter_gr_center))$distance
gm12878$Gm12878_PoisedPromoter_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_PoisedPromoter_gr_center))$distance
gm12878$Gm12878_StrongEnhancer4_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_StrongEnhancer4_gr_center))$distance
gm12878$Gm12878_StrongEnhancer5_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_StrongEnhancer5_gr_center))$distance
gm12878$Gm12878_WeakEnhancer6_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakEnhancer6_gr_center))$distance
gm12878$Gm12878_WeakEnhancer7_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakEnhancer7_gr_center))$distance
gm12878$Gm12878_Insulator_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_Insulator_gr_center))$distance
gm12878$Gm12878_TxnTransition_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_TxnTransition_gr_center))$distance


# Combined

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/Combined/")

temp <- list.files()

#Gm12878
Gm12878_CTCF <- read.table(temp[1],header=FALSE,sep="\t") 
Gm12878_E <- read.table(temp[2],header=FALSE,sep="\t")
Gm12878_PF <- read.table(temp[3],header=FALSE,sep="\t")
Gm12878_R <- read.table(temp[4],header=FALSE,sep="\t")
Gm12878_T <- read.table(temp[5],header=FALSE,sep="\t")
Gm12878_TSS <- read.table(temp[6],header=FALSE,sep="\t")
Gm12878_WE <- read.table(temp[7],header=FALSE,sep="\t")

Gm12878_CTCF_gr <- GRanges(seqnames=Gm12878_CTCF$V1,IRanges(start=Gm12878_CTCF$V2,end=Gm12878_CTCF$V3))
Gm12878_CTCF_gr_center <- resize(Gm12878_CTCF_gr, width = 1, fix = "center")
Gm12878_E_gr <- GRanges(seqnames=Gm12878_E$V1,IRanges(start=Gm12878_E$V2,end=Gm12878_E$V3))
Gm12878_E_gr_center <- resize(Gm12878_E_gr, width = 1, fix = "center")
Gm12878_PF_gr <- GRanges(seqnames=Gm12878_PF$V1,IRanges(start=Gm12878_PF$V2,end=Gm12878_PF$V3))
Gm12878_PF_gr_center <- resize(Gm12878_PF_gr, width = 1, fix = "center")
Gm12878_R_gr <- GRanges(seqnames=Gm12878_R$V1,IRanges(start=Gm12878_R$V2,end=Gm12878_R$V3))
Gm12878_R_gr_center <- resize(Gm12878_R_gr, width = 1, fix = "center")
Gm12878_T_gr <- GRanges(seqnames=Gm12878_T$V1,IRanges(start=Gm12878_T$V2,end=Gm12878_T$V3))
Gm12878_T_gr_center <- resize(Gm12878_T_gr, width = 1, fix = "center")
Gm12878_TSS_gr <- GRanges(seqnames=Gm12878_TSS$V1,IRanges(start=Gm12878_TSS$V2,end=Gm12878_TSS$V3))
Gm12878_TSS_gr_center <- resize(Gm12878_TSS_gr, width = 1, fix = "center")
Gm12878_WE_gr <- GRanges(seqnames=Gm12878_WE$V1,IRanges(start=Gm12878_WE$V2,end=Gm12878_WE$V3))
Gm12878_WE_gr_center <- resize(Gm12878_WE_gr, width = 1, fix = "center")

gm12878$Gm12878_CTCF <- ifelse(countOverlaps(binslist5,Gm12878_CTCF_gr)>=1,1,0) 
gm12878$Gm12878_E <- ifelse(countOverlaps(binslist5,Gm12878_E_gr)>=1,1,0)
gm12878$Gm12878_PF <- ifelse(countOverlaps(binslist5,Gm12878_PF_gr)>=1,1,0)
gm12878$Gm12878_R <- ifelse(countOverlaps(binslist5,Gm12878_R_gr)>=1,1,0)
gm12878$Gm12878_T <- ifelse(countOverlaps(binslist5,Gm12878_T_gr)>=1,1,0)
gm12878$Gm12878_TSS <- ifelse(countOverlaps(binslist5,Gm12878_TSS_gr)>=1,1,0)
gm12878$Gm12878_WE <- ifelse(countOverlaps(binslist5,Gm12878_WE_gr)>=1,1,0)

gm12878$Gm12878_CTCF_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_CTCF_gr_center))$distance
gm12878$Gm12878_E_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_E_gr_center))$distance
gm12878$Gm12878_PF_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_PF_gr_center))$distance
gm12878$Gm12878_R_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_R_gr_center))$distance
gm12878$Gm12878_T_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_T_gr_center))$distance
gm12878$Gm12878_TSS_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_TSS_gr_center))$distance
gm12878$Gm12878_WE_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WE_gr_center))$distance


# DNase I

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DNaseI/")

temp <- list.files()

Gm12878_DNaseI <- read.table(temp[1],header=FALSE,sep="\t") 

Gm12878_DNaseI_gr <- GRanges(seqnames=Gm12878_DNaseI$V1,IRanges(start=Gm12878_DNaseI$V2,end=Gm12878_DNaseI$V3))
Gm12878_DNaseI_gr_center <- resize(Gm12878_DNaseI_gr, width = 1, fix = "center")

gm12878$Gm12878_DNaseI <- ifelse(countOverlaps(binslist5,Gm12878_DNaseI_gr)>=1,1,0)

gm12878$Gm12878_DNaseI_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_DNaseI_gr_center))$distance


# Histone Modifications

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/HistoneModifications/")

temp <- list.files()

Gm12878_H2az <- read.table(temp[1],header=FALSE,sep="\t") 
Gm12878_H3k27ac <- read.table(temp[2],header=FALSE,sep="\t")
Gm12878_H3k27me3 <- read.table(temp[3],header=FALSE,sep="\t")
Gm12878_H3k36me3 <- read.table(temp[4],header=FALSE,sep="\t")
Gm12878_H3k4me1 <- read.table(temp[5],header=FALSE,sep="\t")
Gm12878_H3k4me2 <- read.table(temp[6],header=FALSE,sep="\t")
Gm12878_H3k4me3 <- read.table(temp[7],header=FALSE,sep="\t")
Gm12878_H3k79me2 <- read.table(temp[8],header=FALSE,sep="\t")
Gm12878_H3k9ac <- read.table(temp[9],header=FALSE,sep="\t")
Gm12878_H3k9me3 <- read.table(temp[10],header=FALSE,sep="\t")
Gm12878_H4k20me1 <- read.table(temp[11],header=FALSE,sep="\t")

Gm12878_H2az_gr <- GRanges(seqnames=Gm12878_H2az$V1,IRanges(start=Gm12878_H2az$V2,end=Gm12878_H2az$V3))
Gm12878_H2az_gr_center <- resize(Gm12878_H2az_gr, width = 1, fix = "center")
Gm12878_H3k27ac_gr <- GRanges(seqnames=Gm12878_H3k27ac$V1,IRanges(start=Gm12878_H3k27ac$V2,end=Gm12878_H3k27ac$V3))
Gm12878_H3k27ac_gr_center <- resize(Gm12878_H3k27ac_gr, width = 1, fix = "center")
Gm12878_H3k27me3_gr <- GRanges(seqnames=Gm12878_H3k27me3$V1,IRanges(start=Gm12878_H3k27me3$V2,end=Gm12878_H3k27me3$V3))
Gm12878_H3k27me3_gr_center <- resize(Gm12878_H3k27me3_gr, width = 1, fix = "center")
Gm12878_H3k36me3_gr <- GRanges(seqnames=Gm12878_H3k36me3$V1,IRanges(start=Gm12878_H3k36me3$V2,end=Gm12878_H3k36me3$V3))
Gm12878_H3k36me3_gr_center <- resize(Gm12878_H3k36me3_gr, width = 1, fix = "center")
Gm12878_H3k4me1_gr <- GRanges(seqnames=Gm12878_H3k4me1$V1,IRanges(start=Gm12878_H3k4me1$V2,end=Gm12878_H3k4me1$V3))
Gm12878_H3k4me1_gr_center <- resize(Gm12878_H3k4me1_gr, width = 1, fix = "center")
Gm12878_H3k4me2_gr <- GRanges(seqnames=Gm12878_H3k4me2$V1,IRanges(start=Gm12878_H3k4me2$V2,end=Gm12878_H3k4me2$V3))
Gm12878_H3k4me2_gr_center <- resize(Gm12878_H3k4me2_gr, width = 1, fix = "center")
Gm12878_H3k4me3_gr <- GRanges(seqnames=Gm12878_H3k4me3$V1,IRanges(start=Gm12878_H3k4me3$V2,end=Gm12878_H3k4me3$V3))
Gm12878_H3k4me3_gr_center <- resize(Gm12878_H3k4me3_gr, width = 1, fix = "center")
Gm12878_H3k79me2_gr <- GRanges(seqnames=Gm12878_H3k79me2$V1,IRanges(start=Gm12878_H3k79me2$V2,end=Gm12878_H3k79me2$V3))
Gm12878_H3k79me2_gr_center <- resize(Gm12878_H3k79me2_gr, width = 1, fix = "center")
Gm12878_H3k9ac_gr <- GRanges(seqnames=Gm12878_H3k9ac$V1,IRanges(start=Gm12878_H3k9ac$V2,end=Gm12878_H3k9ac$V3))
Gm12878_H3k9ac_gr_center <- resize(Gm12878_H3k9ac_gr, width = 1, fix = "center")
Gm12878_H3k9me3_gr <- GRanges(seqnames=Gm12878_H3k9me3$V1,IRanges(start=Gm12878_H3k9me3$V2,end=Gm12878_H3k9me3$V3))
Gm12878_H3k9me3_gr_center <- resize(Gm12878_H3k9me3_gr, width = 1, fix = "center")
Gm12878_H4k20me1_gr <- GRanges(seqnames=Gm12878_H4k20me1$V1,IRanges(start=Gm12878_H4k20me1$V2,end=Gm12878_H4k20me1$V3))
Gm12878_H4k20me1_gr_center <- resize(Gm12878_H4k20me1_gr, width = 1, fix = "center")

gm12878$Gm12878_H2az <- ifelse(countOverlaps(binslist5,Gm12878_H2az_gr)>=1,1,0) 
gm12878$Gm12878_H3k27ac <- ifelse(countOverlaps(binslist5,Gm12878_H3k27ac_gr)>=1,1,0)
gm12878$Gm12878_H3k27me3 <- ifelse(countOverlaps(binslist5,Gm12878_H3k27me3_gr)>=1,1,0)
gm12878$Gm12878_H3k36me3 <- ifelse(countOverlaps(binslist5,Gm12878_H3k36me3_gr)>=1,1,0)
gm12878$Gm12878_H3k4me1 <- ifelse(countOverlaps(binslist5,Gm12878_H3k4me1_gr)>=1,1,0)
gm12878$Gm12878_H3k4me2 <- ifelse(countOverlaps(binslist5,Gm12878_H3k4me2_gr)>=1,1,0)
gm12878$Gm12878_H3k4me3 <- ifelse(countOverlaps(binslist5,Gm12878_H3k4me3_gr)>=1,1,0)
gm12878$Gm12878_H3k79me2 <- ifelse(countOverlaps(binslist5,Gm12878_H3k79me2_gr)>=1,1,0)
gm12878$Gm12878_H3k9ac <- ifelse(countOverlaps(binslist5,Gm12878_H3k9ac_gr)>=1,1,0)
gm12878$Gm12878_H3k9me3 <- ifelse(countOverlaps(binslist5,Gm12878_H3k9me3_gr)>=1,1,0)
gm12878$Gm12878_H4k20me1 <- ifelse(countOverlaps(binslist5,Gm12878_H4k20me1_gr)>=1,1,0)

gm12878$Gm12878_H2az_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H2az_gr_center))$distance
gm12878$Gm12878_H3k27ac_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k27ac_gr_center))$distance
gm12878$Gm12878_H3k27me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k27me3_gr_center))$distance
gm12878$Gm12878_H3k36me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k36me3_gr_center))$distance
gm12878$Gm12878_H3k4me1_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k4me1_gr_center))$distance
gm12878$Gm12878_H3k4me2_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k4me2_gr_center))$distance
gm12878$Gm12878_H3k4me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k4me3_gr_center))$distance
gm12878$Gm12878_H3k79me2_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k79me2_gr_center))$distance
gm12878$Gm12878_H3k9ac_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k9ac_gr_center))$distance
gm12878$Gm12878_H3k9me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k9me3_gr_center))$distance
gm12878$Gm12878_H4k20me1_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H4k20me1_gr_center))$distance



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


# Comparing SMOTE to bootstraps

#SMOTE

# Splitting the data
set.seed(5228)
inTrainingSet <- createDataPartition(chr1data_f$y,p=.7,list=FALSE)
train <- chr1data_f[inTrainingSet,]
test <- chr1data_f[-inTrainingSet,]

# Establishing tuning/training parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 3,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#Testing eight different combinations of perc.over/perc.under
#100/200, 200/200, 300/200, 400/200, 
#100/300, 200/300, 300/300, 400/300

#########################################################################################

#ENET

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_sm <- list(tpr <- matrix(nrow=length(test$y), 
                                 ncol=8),
                   fpr <- matrix(nrow=length(test$y), 
                                 ncol=8),
                   auc <- numeric(8),
                   varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                    ncol=8))
rownames(enetlst_sm[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_sm <- matrix(nrow = 16, ncol=8)
rownames(enetperf_sm) <- c("TN",
                           "FN",
                           "FP",
                           "TP",
                           "Total",
                           "Sensitivity",
                           "Specificity",
                           "Kappa",
                           "Accuracy",
                           "Precision",
                           "FPR",
                           "FNR",
                           "FOR",
                           "NPV",
                           "MCC",
                           "F1")



for(i in 1:4){
  set.seed(111)
  #100/200, 200/200, 300/200, 400/200
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 200)
  
  #ENET Model
  enetModel_sm <- train(y ~ ., data=train_smote, 
                        method = "glmnet", 
                        metric="ROC", 
                        trControl = fitControl,
                        family="binomial",
                        tuneLength=5,
                        standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(enetModel_sm, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_sm[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_sm[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_sm[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  enetlst_sm[[4]][,i] <- varImp(enetModel_sm)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_sm,
                             newdata=test,
                             type="raw")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  enetperf_sm[1,i] <- confMat$table[1,1]
  enetperf_sm[2,i] <- confMat$table[1,2]
  enetperf_sm[3,i] <- confMat$table[2,1]
  enetperf_sm[4,i] <- confMat$table[2,2]
  enetperf_sm[5,i] <- sum(confMat$table)
  enetperf_sm[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_sm[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_sm[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_sm[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_sm[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_sm[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_sm[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_sm[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_sm[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_sm[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_sm[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_sm[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
  
  #########################################################################################
  
  set.seed(111)
  #100/300, 200/300, 300/300, 400/300
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 300)
  
  #ENET Model
  enetModel_sm <- train(y ~ ., data=train_smote, 
                        method = "glmnet", 
                        metric="ROC", 
                        trControl = fitControl,
                        family="binomial",
                        tuneLength=5,
                        standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.enetModel <- as.vector(predict(enetModel_sm, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_sm[[1]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_sm[[2]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_sm[[3]][i+4] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  enetlst_sm[[4]][,i+4] <- varImp(enetModel_sm)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_sm[1,i+4] <- confMat$table[1,1]
  enetperf_sm[2,i+4] <- confMat$table[1,2]
  enetperf_sm[3,i+4] <- confMat$table[2,1]
  enetperf_sm[4,i+4] <- confMat$table[2,2]
  enetperf_sm[5,i+4] <- sum(confMat$table)
  enetperf_sm[6,i+4] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_sm[7,i+4] <- as.vector(confMat$byClass["Specificity"])
  enetperf_sm[8,i+4] <- as.vector(confMat$overall["Kappa"])
  enetperf_sm[9,i+4] <- as.vector(confMat$overall["Accuracy"])
  enetperf_sm[10,i+4] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_sm[11,i+4] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_sm[12,i+4] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_sm[13,i+4] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_sm[14,i+4] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_sm[15,i+4] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_sm[15,i+4] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_sm[16,i+4] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

saveRDS(enetlst_sm, "enetlst_sm_lns.RDS")

saveRDS(enetperf_sm, "enetperf_sm.rds")


#Bootstraps

#set number of bootstrap samples
bootsamps = 100

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(chr1_gm12878_f$y[which(chr1_gm12878_f$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_gm12878_f$y=="No"),
                        length(which(chr1_gm12878_f$y=="Yes")),
                        replace = TRUE)
}



#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

########################################################################################

#With log transform and no standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                 ncol=bootsamps))
rownames(enetlst[[4]]) <- colnames(chr1_gm12878_f)[-1]


enetperf_b <- matrix(nrow = 16, ncol=1)
rownames(enetperf_b) <- c("TN",
                          "FN",
                          "FP",
                          "TP",
                          "Total",
                          "Sensitivity",
                          "Specificity",
                          "Kappa",
                          "Accuracy",
                          "Precision",
                          "FPR",
                          "FNR",
                          "FOR",
                          "NPV",
                          "MCC",
                          "F1")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),],
                           chr1_gm12878_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_b[1,i] <- confMat$table[1,1]
  enetperf_b[2,i] <- confMat$table[1,2]
  enetperf_b[3,i] <- confMat$table[2,1]
  enetperf_b[4,i] <- confMat$table[2,2]
  enetperf_b[5,i] <- sum(confMat$table)
  enetperf_b[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_b[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_b[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_b[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_b[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_b[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_b[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_b[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_b[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_sm[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_b[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_b[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst[[3]])

saveRDS(enetlst, "enetlst_bs_lns.rds")

saveRDS(enetperf_b, "enetperf_b.rds")


# Evaluating Normalization


#set number of bootstrap samples
bootsamps = 5

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(chr1_gm12878_f$y[which(chr1_gm12878_f$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_gm12878_f$y=="No"),
                        length(which(chr1_gm12878_f$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

########################################################################################

#With log transform and standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_ls <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   auc <- numeric(bootsamps),
                   varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                    ncol=bootsamps))
rownames(enetlst_ls[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_ls <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_ls) <- c("TN",
                           "FN",
                           "FP",
                           "TP",
                           "Total",
                           "Sensitivity",
                           "Specificity",
                           "Kappa",
                           "Accuracy",
                           "Precision",
                           "FPR",
                           "FNR",
                           "FOR",
                           "NPV",
                           "MCC",
                           "F1")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),],
                           chr1_gm12878_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=TRUE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_ls[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_ls[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_ls[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_ls[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_ls[1,i] <- confMat$table[1,1]
  enetperf_ls[2,i] <- confMat$table[1,2]
  enetperf_ls[3,i] <- confMat$table[2,1]
  enetperf_ls[4,i] <- confMat$table[2,2]
  enetperf_ls[5,i] <- sum(confMat$table)
  enetperf_ls[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_ls[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_ls[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_ls[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_ls[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_ls[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_ls[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_ls[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_ls[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_ls[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_ls[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_ls[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_ls[[3]])

saveRDS(enetlst_ls, "enetlst_ls.rds")
saveRDS(enetperf_ls, "enetperf_ls.rds")

#################################################################################

#With log transform and no standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_lns <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    auc <- numeric(bootsamps),
                    varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_lns[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_lns <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_lns) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),],
                           chr1_gm12878_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_lns[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_lns[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_lns[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_lns[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_lns[1,i] <- confMat$table[1,1]
  enetperf_lns[2,i] <- confMat$table[1,2]
  enetperf_lns[3,i] <- confMat$table[2,1]
  enetperf_lns[4,i] <- confMat$table[2,2]
  enetperf_lns[5,i] <- sum(confMat$table)
  enetperf_lns[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_lns[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_lns[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_lns[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_lns[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_lns[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_lns[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_lns[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_lns[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_lns[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_lns[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_lns[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_lns[[3]])

saveRDS(enetlst_lns, "enetlst_lns.rds")
saveRDS(enetperf_lns, "enetperf_lns.rds")



#################################################################################

#Without log transform and standardization

cols <- c(grep("dist",colnames(chr1_gm12878_f)))
chr1_gm12878_f[,cols] <- apply(chr1_gm12878_f[,cols], 2, function(x){2^x})

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_nls <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    auc <- numeric(bootsamps),
                    varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_nls[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_nls <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_nls) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),],
                           chr1_gm12878_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=TRUE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_nls[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_nls[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_nls[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_nls[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_nls[1,i] <- confMat$table[1,1]
  enetperf_nls[2,i] <- confMat$table[1,2]
  enetperf_nls[3,i] <- confMat$table[2,1]
  enetperf_nls[4,i] <- confMat$table[2,2]
  enetperf_nls[5,i] <- sum(confMat$table)
  enetperf_nls[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_nls[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_nls[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_nls[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_nls[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_nls[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_nls[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_nls[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_nls[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_nls[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_nls[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_nls[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_nls[[3]])

saveRDS(enetlst_nls, "enetlst_nls.rds")
saveRDS(enetperf_nls, "enetperf_nls.rds")


####################################################################

#Without log transform and no standardization

enetlst_nlns <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                   ncol=bootsamps),
                     fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_f$y=="Yes"))*2)*.3), 
                                   ncol=bootsamps),
                     auc <- numeric(bootsamps),
                     varimp <- matrix(nrow=dim(chr1_gm12878_f)[2]-1,
                                      ncol=bootsamps))
rownames(enetlst_nlns[[4]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_nlns <- matrix(nrow = 16, ncol=bootsamps)
rownames(enetperf_nlns) <- c("TN",
                             "FN",
                             "FP",
                             "TP",
                             "Total",
                             "Sensitivity",
                             "Specificity",
                             "Kappa",
                             "Accuracy",
                             "Precision",
                             "FPR",
                             "FNR",
                             "FOR",
                             "NPV",
                             "MCC",
                             "F1")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),],
                           chr1_gm12878_f[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_nlns[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_nlns[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_nlns[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  enetlst_nlns[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_nlns[1,i] <- confMat$table[1,1]
  enetperf_nlns[2,i] <- confMat$table[1,2]
  enetperf_nlns[3,i] <- confMat$table[2,1]
  enetperf_nlns[4,i] <- confMat$table[2,2]
  enetperf_nlns[5,i] <- sum(confMat$table)
  enetperf_nlns[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_nlns[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_nlns[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_nlns[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_nlns[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  enetperf_nlns[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  enetperf_nlns[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  enetperf_nlns[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  enetperf_nlns[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #enetperf_nlns[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  enetperf_nlns[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_nlns[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}

mean(enetlst_nlns[[3]])

saveRDS(enetlst_nlns, "enetlst_nlns.rds")
saveRDS(enetperf_nlns, "enetperf_nlns.rds")


# Performing stepwise selection

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f <- rbind.data.frame(chr1_gm12878_f[samps,],
                                   chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])


#center and scaling data to avoid using intercept term
cols <- names(Filter(is.numeric, chr1_gm12878_f))
chr1_gm12878_f[,cols] <- scale(chr1_gm12878_f[,cols], center = TRUE, scale = TRUE)



#Using cross validation (10 fold)

#forward
k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_gm12878_f), replace=TRUE)
cv.preds.fwd=matrix(NA, nrow=ncol(chr1_gm12878_f)-1,ncol=k)
auc.model.fwd <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_gm12878_f[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_gm12878_f[folds!=j,], family = binomial)
  
  best.fit.fwd = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="forward",
                      trace=0)
  
  numpreds <- length(names(best.fit.fwd$coefficients)[-1])
  cv.preds.fwd[(1:numpreds),j] <- names(best.fit.fwd$coefficients)[-1]
  cols <- names(best.fit.fwd$model)
  model <- glm(y ~ . , data = chr1_gm12878_f[folds==j,cols], family = binomial)
  pred.model <- predict(model, newdata=chr1_gm12878_f[folds==j,cols], type="response")
  roc.model <- roc(chr1_gm12878_f[folds==j,"y"], pred.model)
  auc.model.fwd[j] <- pROC::auc(roc.model)
  
}


vars.fwd <- na.omit(cv.preds.fwd[,which.max(auc.model.fwd)])
vars.fwd[grep("_dist",vars.fwd,invert = TRUE)] <- unlist(lapply(vars.fwd[grep("_dist",vars.fwd,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))

chr1_gm12878_fwd <- chr1_gm12878_f[,which((names(chr1_gm12878_f) %in% vars.fwd) | names(chr1_gm12878_f)=="y")]



# Performing recursive feature elimination

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f <- rbind.data.frame(chr1_gm12878_f[samps,],
                                   chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])


#setting rfe parameters
control <- rfeControl(functions=rfFuncs, method="cv", number=10)#, repeats=5)

trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)


#splitting data
set.seed(7215)
inTrainingSet <- sample(length(chr1_gm12878_f2$y),floor(length(chr1_gm12878_f2$y)*.7))
#inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
train <- chr1_gm12878_f2[inTrainingSet,]
test <- chr1_gm12878_f2[-inTrainingSet,]

rfeModel <- rfe(train[,-1], 
                train[,1], 
                sizes=c(2:50), 
                metric="ROC",
                rfeControl=control,
                trControl = trainctrl)

pred.rfeModel <- as.vector(predict(rfeModel, newdata=test, type="prob")[,"Yes"])

pred.rfeModel2 <- predict(rfeModel,
                          newdata=test,
                          type="raw")

vars <- c("y",predictors(rfeModel))

chr1_gm12878_rfe <- chr1_gm12878_f[,which(names(chr1_gm12878_f) %in% vars)]


# Performing Random Forest

#set number of bootstrap samples

bootsamps = 10

#set tuning parameters

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#create a matrix of row ids that represent the zero class

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(chr1_gm12878_fwd$y[which(chr1_gm12878_fwd$y=="Yes")]))

#filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(chr1_gm12878_fwd$y=="No"),
                        length(which(chr1_gm12878_fwd$y=="Yes")),
                        replace = TRUE)
}


#function for roc curves

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

# Random Forest

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_fwd$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_fwd$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              auc <- numeric(bootsamps),
              varimp <- matrix(nrow=dim(chr1_gm12878_fwd)[2]-1,
                               ncol=bootsamps))
rownames(rflst[[4]]) <- colnames(chr1_gm12878_fwd)[-1]

# Random Forest


#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_fwd$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              fpr <- matrix(nrow=ceiling((length(which(chr1_gm12878_fwd$y=="Yes"))*2)*.3), 
                            ncol=bootsamps),
              auc <- numeric(bootsamps),
              varimp <- matrix(nrow=dim(chr1_gm12878_fwd)[2]-1,
                               ncol=bootsamps))
rownames(rflst[[4]]) <- colnames(chr1_gm12878_fwd)[-1]


rfperf <- matrix(nrow = 16, ncol=bootsamps)
rownames(rfperf) <- c("TN",
                      "FN",
                      "FP",
                      "TP",
                      "Total",
                      "Sensitivity",
                      "Specificity",
                      "Kappa",
                      "Accuracy",
                      "Precision",
                      "FPR",
                      "FNR",
                      "FOR",
                      "NPV",
                      "MCC",
                      "F1")

#Performing Random Forest

for(i in 1:bootsamps){
  
  #combining the two classes to create balanced data
  data <- rbind.data.frame(chr1_gm12878_fwd[which(chr1_gm12878_fwd$y=="Yes"),],
                           chr1_gm12878_fwd[sampids[,i],])
  
  #determining the best number of variables randomly sampled as candidates at each split
  set.seed(5430)
  bestmtry <- tuneRF(data[,-1],data$y,
                     improve=.01,trace=0,plot=F) 
  bestmtry <- data.frame(bestmtry)
  bestmtry <- bestmtry[order(bestmtry$OOBError, decreasing = FALSE),]
  #bestmtry$mtry[1]
  
  #splitting the data
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #determining best number of trees
  tunegrid <- expand.grid(.mtry=bestmtry$mtry[1])
  modellist <- list()
  for (ntree in c(50,200,500,1000)) {
    set.seed(333)
    fit <- train(y~., data=train, 
                 method="rf", 
                 metric="Accuracy",
                 tuneGrid=tunegrid,  
                 ntree=ntree)
    key <- toString(ntree)
    modellist[[key]] <- fit
  }
  # compare results
  results <- resamples(modellist)
  #summary(results)
  #dotplot(results)
  results <- data.frame(summary(results)[3]$statistics$Accuracy)
  results <- results[order(results$Mean, decreasing = TRUE),]
  
  set.seed(1006)
  rfModel <- train(y~., data=train, 
                   method="rf", 
                   metric="ROC", 
                   tuneGrid=tunegrid, 
                   trControl=fitControl, 
                   ntree=as.numeric(rownames(results)[1]))
  
  #Prediction vector for ROC and AUC
  pred.rfModel <- as.vector(predict(rfModel, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  rflst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst[[3]][i] <- pROC::auc(pROC::roc(test$y, pred.rfModel))
  rflst[[4]][,i] <- varImp(rfModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.rfModel2 <- predict(rfModel,
                           newdata=test,
                           type="raw")
  confMat <- confusionMatrix(data=pred.rfModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  rfperf[1,i] <- confMat$table[1,1]
  rfperf[2,i] <- confMat$table[1,2]
  rfperf[3,i] <- confMat$table[2,1]
  rfperf[4,i] <- confMat$table[2,2]
  rfperf[5,i] <- sum(confMat$table)
  rfperf[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  rfperf[7,i] <- as.vector(confMat$byClass["Specificity"])
  rfperf[8,i] <- as.vector(confMat$overall["Kappa"])
  rfperf[9,i] <- as.vector(confMat$overall["Accuracy"])
  rfperf[10,i] <- confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])
  rfperf[11,i] <- confMat$table[2,1]/(confMat$table[2,1]+confMat$table[1,1])
  rfperf[12,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[2,2])
  rfperf[13,i] <- confMat$table[1,2]/(confMat$table[1,2]+confMat$table[1,1])
  rfperf[14,i] <- confMat$table[1,1]/(confMat$table[1,1]+confMat$table[1,2])
  #rfperf[15,i] <- mccr(ifelse(test$y=="Yes",1,0),ifelse(pred.enetModel2=="Yes",1,0))
  rfperf[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  rfperf[16,i] <- (2*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))/(((confMat$table[2,2]/(confMat$table[2,2]+confMat$table[1,2]))*(confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1]))) + (confMat$table[2,2]/(confMat$table[2,2]+confMat$table[2,1])))
  
}


saveRDS(rflst, "rflst.rds")
saveRDS(rfperf, "rfperf.rds")
