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
#0.99061555 0.00938445
dim(tad_train)  
#492730      3

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
