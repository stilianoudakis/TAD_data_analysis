#Adding Annotations

# Reading in the TAD and Subcompartment overlap data (granges object)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData/K562")

tad_subcomp_full_k562 <- readRDS("tad_subcomp_full_k562.rds")

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData/K562")

k562 <- readRDS("k562.rds")

tad_subcomp_dist <- GRanges(seqnames = seqnames(tad_subcomp_full_k562),
                            ranges = IRanges(start = start(tad_subcomp_full_k562)+500,
                                             width=1),
                            strand = "*")

#Adding other covariates of interest from Merlot folder

#3D subcompartments

#There is no 3D subcompartment data for the k562 cell line

#DGV

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/DGV/DGV_files")

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
deletion_gr <- GRanges(seqnames=deletion$V1,IRanges(start=deletion$V2,end=deletion$V3))
duplication_gr <- GRanges(seqnames=duplication$V1,IRanges(start=duplication$V2,end=duplication$V3))
gain_loss_gr <- GRanges(seqnames=gain_loss$V1,IRanges(start=gain_loss$V2,end=gain_loss$V3))
insertion_gr <- GRanges(seqnames=insertion$V1,IRanges(start=insertion$V2,end=insertion$V3))
inversion_gr <- GRanges(seqnames=inversion$V1,IRanges(start=inversion$V2,end=inversion$V3))
mobile_element_insertion_gr <- GRanges(seqnames=mobile_element_insertion$V1,IRanges(start=mobile_element_insertion$V2,end=mobile_element_insertion$V3))
novel_sequence_insertion_gr <- GRanges(seqnames=novel_sequence_insertion$V1,IRanges(start=novel_sequence_insertion$V2,end=novel_sequence_insertion$V3))
sequence_alteration_gr <- GRanges(seqnames=sequence_alteration$V1,IRanges(start=sequence_alteration$V2,end=sequence_alteration$V3))
tandem_duplication_gr <- GRanges(seqnames=tandem_duplication$V1,IRanges(start=tandem_duplication$V2,end=tandem_duplication$V3))

k562$complex <- ifelse(countOverlaps(tad_subcomp_full_k562,complex_gr)>=1,1,0)
k562$deletion <- ifelse(countOverlaps(tad_subcomp_full_k562,deletion_gr)>=1,1,0)
k562$duplication <- ifelse(countOverlaps(tad_subcomp_full_k562,duplication_gr)>=1,1,0)
k562$gain_loss <- ifelse(countOverlaps(tad_subcomp_full_k562,gain_loss_gr)>=1,1,0)
k562$insertion <- ifelse(countOverlaps(tad_subcomp_full_k562,insertion_gr)>=1,1,0)
k562$inversion <- ifelse(countOverlaps(tad_subcomp_full_k562,inversion_gr)>=1,1,0)
k562$mobile_element_insertion <- ifelse(countOverlaps(tad_subcomp_full_k562,mobile_element_insertion_gr)>=1,1,0)
k562$novel_sequence_insertion <- ifelse(countOverlaps(tad_subcomp_full_k562,novel_sequence_insertion_gr)>=1,1,0)
k562$sequence_alteration <- ifelse(countOverlaps(tad_subcomp_full_k562,sequence_alteration_gr)>=1,1,0)
k562$tandem_duplication <- ifelse(countOverlaps(tad_subcomp_full_k562,tandem_duplication_gr)>=1,1,0)

k562$complex_dist <- mcols(distanceToNearest(tad_subcomp_dist, complex_gr))$distance
k562$deletion_dist <- mcols(distanceToNearest(tad_subcomp_dist, deletion_gr))$distance
k562$duplication_dist <- mcols(distanceToNearest(tad_subcomp_dist, duplication_gr))$distance
k562$gain_loss_dist <- mcols(distanceToNearest(tad_subcomp_dist, gain_loss_gr))$distance
k562$insertion_dist <- mcols(distanceToNearest(tad_subcomp_dist, insertion_gr))$distance
k562$inversion_dist <- mcols(distanceToNearest(tad_subcomp_dist, inversion_gr))$distance
k562$mobile_element_insertion_dist <- mcols(distanceToNearest(tad_subcomp_dist, mobile_element_insertion_gr))$distance
k562$novel_sequence_insertion_dist <- mcols(distanceToNearest(tad_subcomp_dist, novel_sequence_insertion_gr))$distance
k562$sequence_alteration_dist <- mcols(distanceToNearest(tad_subcomp_dist, sequence_alteration_gr))$distance
k562$tandem_duplication_dist <- mcols(distanceToNearest(tad_subcomp_dist, tandem_duplication_gr))$distance

# GERP

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/GERP/GERP_hg19.bed")
gerp <- read.table("GERP_hg19.BED",header=FALSE,sep="\t")
gerp_gr <- GRanges(seqnames=gerp$V1,IRanges(start=gerp$V2,end=gerp$V3))
mcols(gerp_gr)$score <- gerp$V5

k562$gerp <- ifelse(countOverlaps(tad_subcomp_full_k562,gerp_gr)>=1,1,0)

k562$gerp_dist <- mcols(distanceToNearest(tad_subcomp_dist, gerp_gr))$distance

#finding which flanks overlap the gerp file so that we can add a score variable
#all other flanks will have a score of 0
which(k562$gerp==1)
k562$gerp_score <- 0
gerpoverlap <- findOverlaps(tad_subcomp_full_k562,gerp_gr)
gerpoverlapdf <- data.frame(queryHits=queryHits(gerpoverlap), score=gerp_gr[subjectHits(gerpoverlap)]$score)
gerpoverlapmean <- aggregate(gerpoverlapdf$score, list(gerpoverlapdf$queryHits), mean)
k562$gerp_score[gerpoverlapmean$Group.1] <- gerpoverlapmean$x

# nestedRepeats

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/nestedRepeats/nestedRepeats_files")

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
line_gr <- GRanges(seqnames=line$V1,IRanges(start=line$V2,end=line$V3))
low_complexity_gr <- GRanges(seqnames=low_complexity$V1,IRanges(start=low_complexity$V2,end=low_complexity$V3))
LTR_gr <- GRanges(seqnames=LTR$V1,IRanges(start=LTR$V2,end=LTR$V3))
other_gr <- GRanges(seqnames=other$V1,IRanges(start=other$V2,end=other$V3))
RC_gr <- GRanges(seqnames=RC$V1,IRanges(start=RC$V2,end=RC$V3))
RNA_gr <- GRanges(seqnames=RNA$V1,IRanges(start=RNA$V2,end=RNA$V3))
rRNA_gr <- GRanges(seqnames=rRNA$V1,IRanges(start=rRNA$V2,end=rRNA$V3))
satellite_gr <- GRanges(seqnames=satellite$V1,IRanges(start=satellite$V2,end=satellite$V3))
scRNA_gr <- GRanges(seqnames=scRNA$V1,IRanges(start=scRNA$V2,end=scRNA$V3))
simple_repeat_gr <- GRanges(seqnames=simple_repeat$V1,IRanges(start=simple_repeat$V2,end=simple_repeat$V3))
SINE_gr <- GRanges(seqnames=SINE$V1,IRanges(start=SINE$V2,end=SINE$V3))
snRNA_gr <- GRanges(seqnames=snRNA$V1,IRanges(start=snRNA$V2,end=snRNA$V3))
srpRNA_gr <- GRanges(seqnames=srpRNA$V1,IRanges(start=srpRNA$V2,end=srpRNA$V3))
tRNA_gr <- GRanges(seqnames=tRNA$V1,IRanges(start=tRNA$V2,end=tRNA$V3))
unknown_gr <- GRanges(seqnames=unknown$V1,IRanges(start=unknown$V2,end=unknown$V3))

k562$DNA <- ifelse(countOverlaps(tad_subcomp_full_k562,DNA_gr)>=1,1,0)
k562$line <- ifelse(countOverlaps(tad_subcomp_full_k562,line_gr)>=1,1,0)
k562$low_complexity <- ifelse(countOverlaps(tad_subcomp_full_k562,low_complexity_gr)>=1,1,0)
k562$LTR <- ifelse(countOverlaps(tad_subcomp_full_k562,LTR_gr)>=1,1,0)
k562$other <- ifelse(countOverlaps(tad_subcomp_full_k562,other_gr)>=1,1,0)
k562$RC <- ifelse(countOverlaps(tad_subcomp_full_k562,RC_gr)>=1,1,0)
k562$RNA <- ifelse(countOverlaps(tad_subcomp_full_k562,RNA_gr)>=1,1,0)
k562$rRNA <- ifelse(countOverlaps(tad_subcomp_full_k562,rRNA_gr)>=1,1,0)
k562$satellite <- ifelse(countOverlaps(tad_subcomp_full_k562,satellite_gr)>=1,1,0)
k562$scRNA <- ifelse(countOverlaps(tad_subcomp_full_k562,scRNA_gr)>=1,1,0)
k562$simple_repeat <- ifelse(countOverlaps(tad_subcomp_full_k562,simple_repeat_gr)>=1,1,0)
k562$SINE <- ifelse(countOverlaps(tad_subcomp_full_k562,SINE_gr)>=1,1,0)
k562$snRNA <- ifelse(countOverlaps(tad_subcomp_full_k562,snRNA_gr)>=1,1,0)
k562$srpRNA <- ifelse(countOverlaps(tad_subcomp_full_k562,srpRNA_gr)>=1,1,0)
k562$tRNA <- ifelse(countOverlaps(tad_subcomp_full_k562,tRNA_gr)>=1,1,0)
k562$unknown <- ifelse(countOverlaps(tad_subcomp_full_k562,unknown_gr)>=1,1,0)

k562$DNA_dist <- mcols(distanceToNearest(tad_subcomp_dist, DNA_gr))$distance 
k562$line_dist <- mcols(distanceToNearest(tad_subcomp_dist, line_gr))$distance 
k562$low_complexity_dist <- mcols(distanceToNearest(tad_subcomp_dist, low_complexity_gr))$distance 
k562$LTR_dist <- mcols(distanceToNearest(tad_subcomp_dist, LTR_gr))$distance 
k562$other_dist <- mcols(distanceToNearest(tad_subcomp_dist, other_gr))$distance 
k562$RC_dist <- mcols(distanceToNearest(tad_subcomp_dist, RC_gr))$distance 
#k562$RNA_dist <- mcols(distanceToNearest(tad_subcomp_dist, RNA_gr))$distance 
#k562$rRNA_dist <- mcols(distanceToNearest(tad_subcomp_dist, rRNA_gr))$distance 
k562$satellite_dist <- mcols(distanceToNearest(tad_subcomp_dist, satellite_gr))$distance 
#k562$scRNA_dist <- mcols(distanceToNearest(tad_subcomp_dist, scRNA_gr))$distance 
k562$simple_repeat_dist <- mcols(distanceToNearest(tad_subcomp_dist, simple_repeat_gr))$distance 
k562$SINE_dist <- mcols(distanceToNearest(tad_subcomp_dist, SINE_gr))$distance 
#k562$snRNA_dist <- mcols(distanceToNearest(tad_subcomp_dist, snRNA_gr))$distance 
k562$srpRNA_dist <- mcols(distanceToNearest(tad_subcomp_dist, srpRNA_gr))$distance 
k562$tRNA_dist <- mcols(distanceToNearest(tad_subcomp_dist, tRNA_gr))$distance 
k562$unknown_dist <- mcols(distanceToNearest(tad_subcomp_dist, unknown_gr))$distance 

#super_enhancers

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/super_enhancers")

temp = list.files()

se_k562 <- read.table(temp[2],header=FALSE,sep="\t")

se_k562_gr <- GRanges(seqnames=se_k562$V1,IRanges(start=se_k562$V2,end=se_k562$V3))

k562$se_k562 <- ifelse(countOverlaps(tad_subcomp_full_k562,se_k562_gr)>=1,1,0)

k562$se_k562_dist <- mcols(distanceToNearest(tad_subcomp_dist, se_k562_gr))$distance

# UCNE

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/UCNEs/UCNE.bed")

UCNE <- read.table("UCNE.BED",header=FALSE,sep="\t")
UCNE_gr <- GRanges(seqnames=UCNE$V1,IRanges(start=UCNE$V2,end=UCNE$V3))
mcols(UCNE_gr)$score <- UCNE$V5

k562$UCNE <- ifelse(countOverlaps(tad_subcomp_full_k562,UCNE_gr)>=1,1,0)

k562$UCNE_dist <- mcols(distanceToNearest(tad_subcomp_dist, UCNE_gr))$distance

#finding which flanks overlap the unce file so that we can add a score variable
#all other flanks will have a score of 0
which(k562$UCNE==1)
k562$UCNE_score <- 0
UCNEoverlap <- findOverlaps(tad_subcomp_full_k562,UCNE_gr)
UCNEoverlapdf <- data.frame(queryHits=queryHits(UCNEoverlap), score=UCNE_gr[subjectHits(UCNEoverlap)]$score)
UCNEoverlapmean <- aggregate(UCNEoverlapdf$score, list(UCNEoverlapdf$queryHits), mean)
k562$UCNE_score[UCNEoverlapmean$Group.1] <- UCNEoverlapmean$x

# VMR

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/VMRs")

VMR <- read.table("VMR_hg19.BED",header=FALSE,sep="\t")
VMR_gr <- GRanges(seqnames=VMR$V1,IRanges(start=VMR$V2,end=VMR$V3))

k562$VMR <- ifelse(countOverlaps(tad_subcomp_full_k562,VMR_gr)>=1,1,0)

k562$VMR_dist <- mcols(distanceToNearest(tad_subcomp_dist, VMR_gr))$distance

# BroadHMM

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/BroadHmm")

temp <- list.files()
temp <- temp[grep("K562", temp)]

#k562
k562_TxnElongation <- read.table(temp[2],header=FALSE,sep="\t")
k562_WeakTxn <- read.table(temp[3],header=FALSE,sep="\t")
k562_Repressed <- read.table(temp[4],header=FALSE,sep="\t")
k562_Heterochromlo <- read.table(temp[5],header=FALSE,sep="\t")
k562_RepetitiveCNV14 <- read.table(temp[6],header=FALSE,sep="\t")
k562_RepetitiveCNV15 <- read.table(temp[7],header=FALSE,sep="\t")
k562_ActivePromoter <- read.table(temp[8],header=FALSE,sep="\t")
k562_WeakPromoter <- read.table(temp[9],header=FALSE,sep="\t")
k562_PoisedPromoter <- read.table(temp[10],header=FALSE,sep="\t")
k562_StrongEnhancer4 <- read.table(temp[11],header=FALSE,sep="\t")
k562_StrongEnhancer5 <- read.table(temp[12],header=FALSE,sep="\t") 
k562_WeakEnhancer6 <- read.table(temp[13],header=FALSE,sep="\t")
k562_WeakEnhancer7 <- read.table(temp[14],header=FALSE,sep="\t")
k562_Insulator <- read.table(temp[15],header=FALSE,sep="\t")
k562_TxnTransition <- read.table(temp[16],header=FALSE,sep="\t")

k562_TxnElongation_gr <- GRanges(seqnames=k562_TxnElongation$V1,IRanges(start=k562_TxnElongation$V2,end=k562_TxnElongation$V3))
k562_WeakTxn_gr <- GRanges(seqnames=k562_WeakTxn$V1,IRanges(start=k562_WeakTxn$V2,end=k562_WeakTxn$V3))
k562_Repressed_gr <- GRanges(seqnames=k562_Repressed$V1,IRanges(start=k562_Repressed$V2,end=k562_Repressed$V3))
k562_Heterochromlo_gr <- GRanges(seqnames=k562_Heterochromlo$V1,IRanges(start=k562_Heterochromlo$V2,end=k562_Heterochromlo$V3)) 
k562_RepetitiveCNV14_gr <- GRanges(seqnames=k562_RepetitiveCNV14$V1,IRanges(start=k562_RepetitiveCNV14$V2,end=k562_RepetitiveCNV14$V3)) 
k562_RepetitiveCNV15_gr <- GRanges(seqnames=k562_RepetitiveCNV15$V1,IRanges(start=k562_RepetitiveCNV15$V2,end=k562_RepetitiveCNV15$V3))
k562_ActivePromoter_gr <- GRanges(seqnames=k562_ActivePromoter$V1,IRanges(start=k562_ActivePromoter$V2,end=k562_ActivePromoter$V3))
k562_WeakPromoter_gr <- GRanges(seqnames=k562_WeakPromoter$V1,IRanges(start=k562_WeakPromoter$V2,end=k562_WeakPromoter$V3))
k562_PoisedPromoter_gr <- GRanges(seqnames=k562_PoisedPromoter$V1,IRanges(start=k562_PoisedPromoter$V2,end=k562_PoisedPromoter$V3)) 
k562_StrongEnhancer4_gr <- GRanges(seqnames=k562_StrongEnhancer4$V1,IRanges(start=k562_StrongEnhancer4$V2,end=k562_StrongEnhancer4$V3)) 
k562_StrongEnhancer5_gr <- GRanges(seqnames=k562_StrongEnhancer5$V1,IRanges(start=k562_StrongEnhancer5$V2,end=k562_StrongEnhancer5$V3))
k562_WeakEnhancer6_gr <- GRanges(seqnames=k562_WeakEnhancer6$V1,IRanges(start=k562_WeakEnhancer6$V2,end=k562_WeakEnhancer6$V3)) 
k562_WeakEnhancer7_gr <- GRanges(seqnames=k562_WeakEnhancer7$V1,IRanges(start=k562_WeakEnhancer7$V2,end=k562_WeakEnhancer7$V3)) 
k562_Insulator_gr <- GRanges(seqnames=k562_Insulator$V1,IRanges(start=k562_Insulator$V2,end=k562_Insulator$V3))
k562_TxnTransition_gr <- GRanges(seqnames=k562_TxnTransition$V1,IRanges(start=k562_TxnTransition$V2,end=k562_TxnTransition$V3)) 

k562$k562_TxnElongation <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_TxnElongation_gr)>=1,1,0)
k562$k562_WeakTxn <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_WeakTxn_gr)>=1,1,0)
k562$k562_Repressed <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_Repressed_gr)>=1,1,0)
k562$k562_Heterochromlo <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_Heterochromlo_gr)>=1,1,0)
k562$k562_RepetitiveCNV14 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_RepetitiveCNV14_gr)>=1,1,0)
k562$k562_RepetitiveCNV15 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_RepetitiveCNV15_gr)>=1,1,0)
k562$k562_ActivePromoter <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_ActivePromoter_gr)>=1,1,0)
k562$k562_WeakPromoter <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_WeakPromoter_gr)>=1,1,0)
k562$k562_PoisedPromoter <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_PoisedPromoter_gr)>=1,1,0)
k562$k562_StrongEnhancer4 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_StrongEnhancer4_gr)>=1,1,0)
k562$k562_StrongEnhancer5 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_StrongEnhancer5_gr)>=1,1,0)
k562$k562_WeakEnhancer6 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_WeakEnhancer6_gr)>=1,1,0)
k562$k562_WeakEnhancer7 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_WeakEnhancer7_gr)>=1,1,0)
k562$k562_Insulator <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_Insulator_gr)>=1,1,0)
k562$k562_TxnTransition <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_TxnTransition_gr)>=1,1,0)

k562$k562_TxnElongation_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_TxnElongation_gr))$distance
k562$k562_WeakTxn_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_WeakTxn_gr))$distance
k562$k562_Repressed_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_Repressed_gr))$distance
k562$k562_Heterochromlo_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_Heterochromlo_gr))$distance
k562$k562_RepetitiveCNV14_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_RepetitiveCNV14_gr))$distance
k562$k562_RepetitiveCNV15_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_RepetitiveCNV15_gr))$distance
k562$k562_ActivePromoter_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_ActivePromoter_gr))$distance
k562$k562_WeakPromoter_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_WeakPromoter_gr))$distance
k562$k562_PoisedPromoter_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_PoisedPromoter_gr))$distance
k562$k562_StrongEnhancer4_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_StrongEnhancer4_gr))$distance
k562$k562_StrongEnhancer5_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_StrongEnhancer5_gr))$distance
k562$k562_WeakEnhancer6_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_WeakEnhancer6_gr))$distance
k562$k562_WeakEnhancer7_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_WeakEnhancer7_gr))$distance
k562$k562_Insulator_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_Insulator_gr))$distance
k562$k562_TxnTransition_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_TxnTransition_gr))$distance

#Combined

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/Combined")

temp <- list.files()
temp <- temp[grep("K562", temp)]

#k562
k562_CTCF <- read.table(temp[2],header=FALSE,sep="\t") 
k562_E <- read.table(temp[3],header=FALSE,sep="\t")
k562_PF <- read.table(temp[4],header=FALSE,sep="\t")
k562_R <- read.table(temp[5],header=FALSE,sep="\t")
k562_T <- read.table(temp[6],header=FALSE,sep="\t")
k562_TSS <- read.table(temp[7],header=FALSE,sep="\t")
k562_WE <- read.table(temp[8],header=FALSE,sep="\t")

k562_CTCF_gr <- GRanges(seqnames=k562_CTCF$V1,IRanges(start=k562_CTCF$V2,end=k562_CTCF$V3))
k562_E_gr <- GRanges(seqnames=k562_E$V1,IRanges(start=k562_E$V2,end=k562_E$V3))
k562_PF_gr <- GRanges(seqnames=k562_PF$V1,IRanges(start=k562_PF$V2,end=k562_PF$V3))
k562_R_gr <- GRanges(seqnames=k562_R$V1,IRanges(start=k562_R$V2,end=k562_R$V3))
k562_T_gr <- GRanges(seqnames=k562_T$V1,IRanges(start=k562_T$V2,end=k562_T$V3))
k562_TSS_gr <- GRanges(seqnames=k562_TSS$V1,IRanges(start=k562_TSS$V2,end=k562_TSS$V3))
k562_WE_gr <- GRanges(seqnames=k562_WE$V1,IRanges(start=k562_WE$V2,end=k562_WE$V3))

k562$k562_CTCF <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_CTCF_gr)>=1,1,0) 
k562$k562_E <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_E_gr)>=1,1,0)
k562$k562_PF <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_PF_gr)>=1,1,0)
k562$k562_R <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_R_gr)>=1,1,0)
k562$k562_T <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_T_gr)>=1,1,0)
k562$k562_TSS <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_TSS_gr)>=1,1,0)
k562$k562_WE <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_WE_gr)>=1,1,0)

k562$k562_CTCF_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_CTCF_gr))$distance
k562$k562_E_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_E_gr))$distance
k562$k562_PF_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_PF_gr))$distance
k562$k562_R_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_R_gr))$distance
k562$k562_T_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_T_gr))$distance
k562$k562_TSS_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_TSS_gr))$distance
k562$k562_WE_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_WE_gr))$distance

# DNase I

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/DNaseI")

temp <- list.files()

k562_DNaseI <- read.table(temp[2],header=FALSE,sep="\t") 

k562_DNaseI_gr <- GRanges(seqnames=k562_DNaseI$V1,IRanges(start=k562_DNaseI$V2,end=k562_DNaseI$V3))

k562$k562_DNaseI <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_DNaseI_gr)>=1,1,0)

k562$k562_DNaseI_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_DNaseI_gr))$distance


# Histone Modifications

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/HistoneModifications")

temp <- list.files()

k562_H2az <- read.table(temp[12],header=FALSE,sep="\t") 
k562_H3k27ac <- read.table(temp[13],header=FALSE,sep="\t")
k562_H3k27me3 <- read.table(temp[14],header=FALSE,sep="\t")
k562_H3k36me3 <- read.table(temp[15],header=FALSE,sep="\t")
k562_H3k4me1 <- read.table(temp[16],header=FALSE,sep="\t")
k562_H3k4me2 <- read.table(temp[17],header=FALSE,sep="\t")
k562_H3k4me3 <- read.table(temp[18],header=FALSE,sep="\t")
k562_H3k79me2 <- read.table(temp[19],header=FALSE,sep="\t")
k562_H3k9ac <- read.table(temp[20],header=FALSE,sep="\t")
k562_H3k9me1 <- read.table(temp[21],header=FALSE,sep="\t")
k562_H3k9me3 <- read.table(temp[23],header=FALSE,sep="\t")
k562_H4k20me1 <- read.table(temp[23],header=FALSE,sep="\t")

k562_H2az_gr <- GRanges(seqnames=k562_H2az$V1,IRanges(start=k562_H2az$V2,end=k562_H2az$V3))
k562_H3k27ac_gr <- GRanges(seqnames=k562_H3k27ac$V1,IRanges(start=k562_H3k27ac$V2,end=k562_H3k27ac$V3))
k562_H3k27me3_gr <- GRanges(seqnames=k562_H3k27me3$V1,IRanges(start=k562_H3k27me3$V2,end=k562_H3k27me3$V3))
k562_H3k36me3_gr <- GRanges(seqnames=k562_H3k36me3$V1,IRanges(start=k562_H3k36me3$V2,end=k562_H3k36me3$V3))
k562_H3k4me1_gr <- GRanges(seqnames=k562_H3k4me1$V1,IRanges(start=k562_H3k4me1$V2,end=k562_H3k4me1$V3))
k562_H3k4me2_gr <- GRanges(seqnames=k562_H3k4me2$V1,IRanges(start=k562_H3k4me2$V2,end=k562_H3k4me2$V3))
k562_H3k4me3_gr <- GRanges(seqnames=k562_H3k4me3$V1,IRanges(start=k562_H3k4me3$V2,end=k562_H3k4me3$V3))
k562_H3k79me2_gr <- GRanges(seqnames=k562_H3k79me2$V1,IRanges(start=k562_H3k79me2$V2,end=k562_H3k79me2$V3))
k562_H3k9ac_gr <- GRanges(seqnames=k562_H3k9ac$V1,IRanges(start=k562_H3k9ac$V2,end=k562_H3k9ac$V3))
k562_H3k9me1_gr <- GRanges(seqnames=k562_H3k9me1$V1,IRanges(start=k562_H3k9me1$V2,end=k562_H3k9me1$V3))
k562_H3k9me3_gr <- GRanges(seqnames=k562_H3k9me3$V1,IRanges(start=k562_H3k9me3$V2,end=k562_H3k9me3$V3))
k562_H4k20me1_gr <- GRanges(seqnames=k562_H4k20me1$V1,IRanges(start=k562_H4k20me1$V2,end=k562_H4k20me1$V3))

k562$k562_H2az <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H2az_gr)>=1,1,0) 
k562$k562_H3k27ac <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k27ac_gr)>=1,1,0)
k562$k562_H3k27me3 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k27me3_gr)>=1,1,0)
k562$k562_H3k36me3 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k36me3_gr)>=1,1,0)
k562$k562_H3k4me1 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k4me1_gr)>=1,1,0)
k562$k562_H3k4me2 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k4me2_gr)>=1,1,0)
k562$k562_H3k4me3 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k4me3_gr)>=1,1,0)
k562$k562_H3k79me2 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k79me2_gr)>=1,1,0)
k562$k562_H3k9ac <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k9ac_gr)>=1,1,0)
k562$k562_H3k9me1 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k9me1_gr)>=1,1,0)
k562$k562_H3k9me3 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H3k9me3_gr)>=1,1,0)
k562$k562_H4k20me1 <- ifelse(countOverlaps(tad_subcomp_full_k562,k562_H4k20me1_gr)>=1,1,0)

k562$k562_H2az_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H2az_gr))$distance
k562$k562_H3k27ac_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k27ac_gr))$distance
k562$k562_H3k27me3_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k27me3_gr))$distance
k562$k562_H3k36me3_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k36me3_gr))$distance
k562$k562_H3k4me1_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k4me1_gr))$distance
k562$k562_H3k4me2_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k4me2_gr))$distance
k562$k562_H3k4me3_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k4me3_gr))$distance
k562$k562_H3k79me2_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k79me2_gr))$distance
k562$k562_H3k9ac_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k9ac_gr))$distance
k562$k562_H3k9me1_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k9me1_gr))$distance
k562$k562_H3k9me3_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H3k9me3_gr))$distance
k562$k562_H4k20me1_dist <- mcols(distanceToNearest(tad_subcomp_dist, k562_H4k20me1_gr))$distance


# Adding Chromosome information to the data

k562$CHR <- seqnames(tad_subcomp_full_k562)

k562$CHR <- as.character(k562$CHR)

# Saving the data

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData/K562")

saveRDS(k562, "k562.rds")