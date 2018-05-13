#################################################################################

#Creating 1kb bin from min-500 to max+500 for each chromosome
#A flanked TAD boundary will be represented as these bins

#bins <- rep( list(GRangesList()), length(unique(coords$Chromosome)) )

#for(i in 1:length(unique(coords$Chromosome))){
#seqn <- unique(coords$Chromosome)[i]

#midpt <- ((min(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])-500) + 
#            (max(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])+500))/2

#bins[[i]] <- GRanges(seqnames=seqn, 
#                ranges=IRanges(start=seq(min(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])-500, 
#                                         max(coords$coordinate[which(coords$Chromosome==unique(coords$Chromosome)[i])])+500,
#                                         1000),
#                               width=1000))

#}

#binslist <- GRangesList(bins)
#binslist <- unlist(binslist)

#Creating an indicator variable Y that denotes whether the tad boundary overlaps the genomic bin 

#y <- countOverlaps(binslist, bounds)
#length(y)
#table(y)
#prop.table(table(y))

# adding the y vector to the bin granges
#mcols(binslist)$y <- y

#################################################################################
X_A <- countOverlaps(chr22_bins, Asubcompint, type = "any")
table(X_A)
X_A2 <- countOverlaps(chr22_bins, Asubcompint, type = "within")
table(X_A2)
#there are 410015-410000=15 intervals for subcompartment A that
#have partial overlapping with genomic bins

#finding percentage of partial overlaps for compartment A
partialoverlaps <- setdiff(findOverlaps(chr22_bins, Asubcompint, type = "any"),findOverlaps(chr22_bins, Asubcompint, type = "within"))
hits <- pintersect(chr22_bins[queryHits(partialoverlaps)], Asubcompint[subjectHits(partialoverlaps)])
hits
#it appears all partial overlaps land on the very end coordinate of the
#genomic bins
#thus the percentage overlap would be 1/50=.02
percentOverlap <- width(hits) / width(chr22_bins[queryHits(partialoverlaps)])
percentOverlap

#changing the partial overlaps to percentages
X_A[queryHits(partialoverlaps)] <- .02
table(X_A)
```

#Creating predictor vectors corresponding to if there is an overlap for compartment B
```{r}

X_B <- countOverlaps(chr22_bins, Bsubcompint, type = "any")
table(X_B)
#some intervals for subcompartment B have multiple overlaps in one genomic bin

X_B2 <- countOverlaps(chr22_bins, Bsubcompint, type = "within")
table(X_B2)
#there are 258007+7*2-258000=21 intervals for subcompartment B 
#that have partial overlapping with genomic bins

#Finding the percentage of partial overlap for compartment B
partialoverlaps <- setdiff(findOverlaps(chr22_bins, Bsubcompint, type = "any"),
                           findOverlaps(chr22_bins, Bsubcompint, type = "within"))
hits <- pintersect(chr22_bins[queryHits(partialoverlaps)], 
                   Bsubcompint[subjectHits(partialoverlaps)])
hits
#it appears all partial overlaps land on the very end coordinate of the
#genomic bins
#thus the percentage overlap would be 1/50=.02
percentOverlap <- width(hits) / width(chr22_bins[queryHits(partialoverlaps)])
percentOverlap

#changing the partial overlaps to percentages
X_B[queryHits(partialoverlaps)] <- .02
table(X_B)