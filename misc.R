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



#looking at within overlaps to find interval for distance calculation

#A compartment
wA <- distanceToNearest(within, Asubcompartgr)
A_dist <- mcols(wA)$distance

#B compartment
wB <- distanceToNearest(within, Bsubcompartgr)
B_dist <- mcols(wB)$distance

#closest distance
#wClosest <- distanceToNearest(within, subcompartgr)
#Closest_dist <- mcols(wClosest)$distance

#assigning the distances to the within granges object
mcols(within)$A_dist <- A_dist
mcols(within)$B_dist <- B_dist
#mcols(within)$Closest_dist <- Closest_dist


#looking at the bins that had no overlaps and calculating distances
#none_mid <- (start(none)+end(none))/2
#none_mid_gr <- GRanges(seqnames = seqnames(none),
#                       ranges = IRanges(start = none_mid, end = none_mid))

#left distance (using precede function)
pre <- precede(none_mid_gr,subcompartgr, select="all")
pre <- pre[order(queryHits(pre), -subjectHits(pre))]
pre <- pre[!duplicated(queryHits(pre))]
pre_sub <- subcompartgr[subjectHits(pre)]
pre_dist <- distance(none_mid_gr,pre_sub)

#right distance
foll <- follow(none_mid_gr,subcompartgr, select="all")
foll <- foll[order(queryHits(foll), -subjectHits(foll))]
foll <- foll[!duplicated(queryHits(foll))]
foll_sub <- subcompartgr[subjectHits(foll)]
foll_dist <- distance(none_mid_gr, foll_sub)


#closest distance
closest_dist <- distanceToNearest(none_mid_gr, subcompartgr)
closest_dist <- mcols(closest_dist)$distance

#assigning the distances to the none granges object
mcols(none)$left_dist <- foll_dist
mcols(none)$right_dist <- pre_dist
mcols(none)$closest_dist <- closest_dist


#assigning 0 for distances for partial overlaps
mcols(partial1)$left_dist <- 0
mcols(partial1)$right_dist <- 0
mcols(partial1)$closest_dist <- 0

mcols(partial2)$left_dist <- 0
mcols(partial2)$right_dist <- 0
mcols(partial2)$closest_dist <- 0



grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,
                                     0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                           C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                                 1, 1.5, 2,5))
svmModel <- train(y ~., data = train, 
                  method = "svmRadial",
                  trControl=fitControl,
                  metric="ROC",
                  tuneGrid = grid_radial)
pred.svmModel <- as.vector(predict(svmModel, 
                                   newdata=test, 
                                   type="prob")[,"Yes"])
roc.svmModel <- pROC::roc(test$y, pred.svmModel)
auclst[[1]][i] <- pROC::auc(roc.svmModel)
#svm varimp
performlst[[1]][,i] <- varImp(svmModel)$importance[,1]



