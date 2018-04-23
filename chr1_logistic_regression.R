#The data comes from the paper by Rao, Huntley, et. al titled,
#"A three-dimensional map of the human genome at kilobase resolution reveals prinicples of chromatin looping"
#GEO accession: GSE63525

#The contact domain annoation file is named: GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt
#The arrowhead algorithm was used to determine domains

#The columns are defined as follows: 
#chromosome1    x1    x2    chromosome2    y1    y2    color    corner_score    Uvar    Lvar    Usign    Lsign
#Explanations of each field are as follows:
#chromosome = the chromosome that the domain is located on
#x1,x2/y1,y2 = the interval spanned by the domain (contact domains manifest as squares on the diagonal of a Hi-C matrix and as such: x1=y1, x2=y2)
#color = the color that the feature will be rendered as if loaded in Juicebox 

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data_analysis")

domains <- read.table("GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", header=T)
domains$chr1 <- as.character(domains$chr1)
domains$chr2 <- as.character(domains$chr2)
head(domains)
dim(domains)
domains_order <- domains[order(domains$chr1,domains$x2),] 

#Focusing hust on chromosome 1
chr1_domains <- domains_order[which(domains_order$chr1=="1"),]
dim(chr1_domains)
#889 12

#Creating bins of 50 kb for chromosome 1 as suggested in Mourad paper
#the last coordinate of chromosome 1 from the 1kb contact matrix was 249238000
#???How do we know what coordinate to go to
start <- seq(0, 249237000, by=50)
chr1_bins <- GRanges(seqnames="chr1", ranges=IRanges(start=start,width=50))

#Creating a GRanges object out of the tads for chromosome 1
#Flanking boundaries by 1kb (500 kb on either side of boundary coordinate)
coords <- c(chr1_domains[,2],chr1_domains[,3])
chr1_coords <- data.frame(chr="chr1", coordinates=coords)
chr1_coords <- chr1_coords[order(chr1_coords$coordinates),]
dim(chr1_coords)
#remove duplicates for coordinates that are conjoined
chr1_coords <- chr1_coords[!duplicated(chr1_coords), ]
dim(chr1_coords)
chr1_bounds <- GRanges(seqnames=chr1_coords$chr, ranges=IRanges(start=chr1_coords$coordinates, width=1))
chr1_bounds <- resize(chr1_bounds, 1001, fix = "center")

#Creating an indicator variable Y that denotes whether the tad boundary is inside the 1 kb bin
#for chr 1
y <- countOverlaps(chr1_bins, chr1_bounds)
length(y)
table(y)
prop.table(table(y))

#Creating predictor variables from subcompartments
subcompart <- read.table("GSE63525_GM12878_subcompartments.BED",header = F)
subcompart <- subcompart[,1:4]
names(subcompart) <- c("chr", "x1", "x2", "compartment")
dim(subcompart)
subcompart$compartment <- as.character(subcompart$compartment)
subcompart$compartment[which(is.na(subcompart$compartment))] <- "N"
table(subcompart$compartment)
#???Should NA compartments be removed or coded as N

#combining compartments
subcompart$compartment[grep("A",subcompart$compartment)] <- "A"
subcompart$compartment[grep("B",subcompart$compartment)] <- "B"
table(subcompart$compartment)

subcompart$compartment <- factor(subcompart$compartment)
subcompart$compartment <- relevel(subcompart$compartment, "N")

#Crearing indicator variables from levels of subcompartments
dummy <- dummyVars(~compartment,
                   data=subcompart,
                   fullRank=T)
subcompart2 <- data.frame(predict(dummy,newdata=subcompart))

subcompart3 <- cbind.data.frame(subcompart[,1:3],subcompart2)
head(subcompart3)

#for chromosome 1
chr1_subcomp <- subcompart3[which(subcompart3$chr=="chr1"),]
Acomp <- chr1_subcomp[which(chr1_subcomp$compartment.A==1),]
Bcomp <- chr1_subcomp[which(chr1_subcomp$compartment.B==1),]

Acompint <- GRanges(seqnames = "chr1", ranges = IRanges(start=Acomp$x1,end=Acomp$x2))
Bcompint <- GRanges(seqnames = "chr1", ranges = IRanges(start=Bcomp$x1,end=Bcomp$x2))

#Creating predictor vectors corresponding to if there 
#is an overlap for compartment A
X_A <- countOverlaps(chr1_bins, Acompint, type = "any")
table(X_A)
#0       1       2 
#2771866 2212851      24 
#some intervals for subcompartment A have multiple overlaps in one 
#genomic bin
X_A2 <- countOverlaps(chr1_bins, Acompint, type = "within")
table(X_A2)
#0       1 
#2772000 2212741 
#there are 2212851+24*2-2212741=158 intervals for subcompartment A 
#that have partial overlapping with genomic bins

#Finding the percentage of partial overlap for compartment A
partialoverlaps <- setdiff(findOverlaps(chr1_bins, Acompint, type = "any"),
                           findOverlaps(chr1_bins, Acompint, type = "within"))
hits <- pintersect(chr1_bins[queryHits(partialoverlaps)], 
                   Acompint[subjectHits(partialoverlaps)])
percentOverlap <- width(hits) / width(chr1_bins[queryHits(partialoverlaps)])
percentOverlap
#???it appears all partial overlaps land on the very end coordinate of the
#genomic bins
#thus the percentage overlap would be 1/50=.02

#changing the partial overlaps to percentages
X_A[queryHits(partialoverlaps)] <- .02
table(X_A)


#Creating predictor vectors corresponding to if there 
#is an overlap for compartment B
X_B <- countOverlaps(chr1_bins, Bcompint, type = "any")
table(X_B)
#0       1       2 
#2796597 2188127      17
#some intervals for subcompartment B also have multiple overlaps in one 
#genomic bin
X_B2 <- countOverlaps(chr1_bins, Bcompint, type = "within")
table(X_B2)
#0       1 
#2796741 2188000 
#there are 2188127+17*2-2188000=161 intervals for subcompartment B 
#that have partial overlapping with genomic bins

#Finding the percentage of partial overlap for compartment B
partialoverlaps <- setdiff(findOverlaps(chr1_bins, Bcompint, type = "any"),
                           findOverlaps(chr1_bins, Bcompint, type = "within"))
hits <- pintersect(chr1_bins[queryHits(partialoverlaps)], 
                   Bcompint[subjectHits(partialoverlaps)])
percentOverlap <- width(hits) / width(chr1_bins[queryHits(partialoverlaps)])
percentOverlap
#???it appears all partial overlaps land on the very end coordinate of the
#genomic bins
#thus the percentage overlap would be 1/50=.02

#changing the partial overlaps to percentages
X_B[queryHits(partialoverlaps)] <- .02
table(X_B)


#Performing logistic regression
tad_data <- data.frame(class=y, A=X_A, B=X_B)
dim(tad_data)
#shuffle the data
set.seed(123)
g <- runif(nrow(tad_data))
tad_data <- tad_data[order(g),]
dim(tad_data)
#Splitting the data
inTrainingSet <- createDataPartition(tad_data$class,
                                     p=.7,
                                     list=FALSE)
tad_train <- tad_data[inTrainingSet,]
tad_test <- tad_data[-inTrainingSet,]
prop.table(table(tad_train$class))   
#0           1 
#0.993145367 0.006854633
dim(tad_train)  
#3489319       3

ctrl <- trainControl(method = "repeatedcv", 
                     number = 10, 
                     savePredictions = TRUE)

modelfit <- glm(class ~ A + B + A:B,
                data = tad_train,
                family = "binomial")
modelfit
summary(modelfit)
pred = predict(modelfit, newdata=tad_test, type="response")





#Using caret
#tad_train$interaction <- tad_train$A*tad_train$B
#tad_train$interaction <- tad_train$A*tad_train$B
tad_train$class <- factor(tad_train$class)
tad_test$class <- factor(tad_test$class)
modelfit2 <- train(class~A+B,  data=tad_train, 
                 method="glm", 
                 family="binomial") #,
                 #trControl = ctrl,
                 #tuneLength = 5)
pred = predict(modelfit2, newdata=tad_test)
accuracy <- table(pred, tad_test[,1])
sum(diag(accuracy))/sum(accuracy)

pred = predict(mod_fit, newdata=tad_test)
confusionMatrix(data=pred, tad_test$class)

library(pROC)
f1 = roc(class ~ A+B, data=tad_train) 
plot(f1, col="red")

library(ROCR)
prob <- predict(modelfit, newdata=tad_test, type="response")
pred <- prediction(prob, tad_test$class)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
