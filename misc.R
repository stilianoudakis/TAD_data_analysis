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





# Decision Tree

```{r include=FALSE}
#form <- as.formula(y ~ .)
#tree.1 <- rpart(form,
#                data=train) #, control=rpart.control(minsplit=20,cp=0))

#prp(tree.1)					# Will plot the tree
#prp(tree.1,varlen=3)

#new.tree.1 <- prp(tree.1,snip=TRUE)$obj # interactively trim the tree
#prp(new.tree.1) # display the new tree
#fancyRpartPlot(new.tree.1)

#tree.2 <- rpart(form,data=train,control=rpart.control(minsplit=20,cp=0.007))			# A more reasonable tree
#prp(tree.2)                                     # A fast plot													
#fancyRpartPlot(tree.2)

```

# CART model

```{r}

set.seed(2014)

#train_smote$y <- as.numeric(as.character(train_smote$y))
#train_smote$y <- factor(train_smote$y)
#levels(train$y) <- c("No", "Yes")

cartModel <- train(y ~ ., data=train, method = "rpart", metric="ROC", trControl = fitControl, tuneLength=5)

pred.cartModel <- as.vector(predict(cartModel, newdata=test, type="prob")[,"Yes"])


roc.cartModel <- pROC::roc(test$y, pred.cartModel)

auc.cartModel <- pROC::auc(roc.cartModel)
#0.7476

```




# SVM

```{r}

#set.seed(5354)

#svmModel <- train(y ~., data = train, method = "svmRadial",
#  trControl=fitControl,
#  metric="ROC",
#  tuneLength = 10)

#grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,
# 0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
# C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
# 1, 1.5, 2,5))
#svmModel <- train(y ~., data = train, method = "svmRadial",
#                    trControl=fitControl,
#                    metric="ROC",
#                    tuneGrid = grid_radial,
#                   tuneLength = 10)


#svm.c <- train(y ~., data = train,
#               method='svmRadial',
#               trControl=fitControl,
#               metric="ROC",tuneLength = 10)

#plot(svm.c)

#pred.svmModel <- as.vector(predict(svm.c, newdata=test, type="prob")[,"Yes"])

#roc.svmModel <- pROC::roc(test$y, pred.svmModel)

#auc.svmModel <- pROC::auc(roc.svmModel)
#0.7659

```


# Neural Network

```{r}

#nnetGrid <- expand.grid(.size=1:10,
#                        .decay=c(0,.1,1,2))
#maxSize <- max(nnetGrid$.size)
#numWts <- 1*(maxSize*(length(dim(train)[2])) + maxSize + 1)

#set.seed(5346)
#nnetModel <- train(y ~ ., data=train,
#                   method="nnet",
#                   metric="ROC",
#                   trControl=fitControl,
#                   tuneLength=5,
#                   trace=FALSE,
#                   maxit=2000)

#plot(nnetModel)

#pred.nnetModel <- as.vector(predict(nnetModel, newdata=test, type="prob")[,"Yes"])

#roc.nnetModel <- pROC::roc(test$y, pred.nnetModel)

#auc.nnetModel <- pROC::auc(roc.nnetModel)

```




# Plotting model performance

```{r}
test.auc <- data.frame(model=c("glm","gbm","glmnet","rForest"),auc=c(auc.glmModel, auc.gbmModel, auc.eNetModel, auc.rfModel))

test.auc <- test.auc[order(test.auc$auc, decreasing=TRUE),]

test.auc$model <- factor(test.auc$model, levels=test.auc$model)

test.auc

p<-ggplot(data=test.auc, aes(x=model, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal()
p

plot(roc.rfModel, legacy.axes=TRUE, col="red")
lines(roc.gbmModel, col="blue")
lines(roc.eNetModel, col="green")



```

# Comparing results when using smote

```{r}
# Using SMOTE


# Splitting the data
set.seed(5228)
inTrainingSet <- createDataPartition(chr1data_f$y,p=.7,list=FALSE)
train <- chr1data_f[inTrainingSet,]
test <- chr1data_f[-inTrainingSet,]

prop.table(table(train$y))   # 0.993434982 0.006565018 
dim(train)  #173343     64

train$y <- as.factor(train$y)
test$y <- as.factor(test$y)
levels(train$y) <- c("No", "Yes")
levels(test$y) <- c("No", "Yes")
#all categorical variables must be factors
str(train)
set.seed(111)
train_smote <- SMOTE(y ~ ., data=train, perc.over = 100, perc.under = 200)
table(train_smote$y)
#  No  Yes 
#2224 2224 
prop.table(table(train_smote$y))   #0.5 0.5 
dim(train_smote)
#4448   62

# Classic GLM method
glmModel_sm <- glm(y ~ ., 
                   data = train_smote, 
                   family = binomial)
pred.glmModel_sm <- predict(glmModel_sm, newdata=test, type="response")
roc.glmModel_sm <- roc(test$y, pred.glmModel_sm)
auc.glmModel_sm <- pROC::auc(roc.glmModel_sm)
#0.7884

# Establishing tuning/training parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

# CART model
#set.seed(2014)
#cartModel_sm <- train(y ~ ., data=train_smote, method = "rpart", metric="ROC", trControl = fitControl, tuneLength=5)
#pred.cartModel_sm <- as.vector(predict(cartModel_sm, newdata=test, type="prob")[,"Yes"])
#roc.cartModel_sm <- pROC::roc(test$y, pred.cartModel_sm)
#auc.cartModel_sm <- pROC::auc(roc.cartModel_sm)
#0.771

# Elastic Net
set.seed(2014)
eNetModel_sm <- train(y ~ ., data=train_smote, method = "glmnet", metric="ROC", trControl = fitControl, family="binomial", tuneLength=5)
pred.eNetModel_sm <- as.vector(predict(eNetModel_sm, newdata=test, type="prob")[,"Yes"])
roc.eNetModel_sm <- pROC::roc(test$y, pred.eNetModel_sm)
auc.eNetModel_sm <- pROC::auc(roc.eNetModel_sm)
#0.7916

# GBM
set.seed(2014)
gbmModel_sm <- train(y ~ ., data=train_smote, method = "gbm", metric="ROC", trControl = fitControl, verbose=FALSE, tuneLength=5)
pred.gbmModel_sm <- as.vector(predict(gbmModel_sm, newdata=test, type="prob")[,"Yes"])
roc.gbmModel_sm <- pROC::roc(test$y, pred.gbmModel_sm)
auc.gbmModel_sm <- pROC::auc(roc.gbmModel_sm)
#0.7815

# Random Forest
set.seed(2014)
rfModel_sm <- train(y ~ ., data=train_smote, method = "rf", metric="ROC", trControl = fitControl, verbose=FALSE, tuneLength=5)
pred.rfModel_sm <- as.vector(predict(rfModel_sm, newdata=test, type="prob")[,"Yes"])
roc.rfModel_sm <- pROC::roc(test$y, pred.rfModel_sm)
auc.rfModel_sm <- pROC::auc(roc.rfModel_sm)
#0.784

#SVM
#set.seed(5354)
#svmModel <- train(y ~., data = train_smote, method = "svmRadial",
#  trControl=fitControl,
#  metric="ROC",
#  tuneLength = 10)
#grid_radial <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,
# 0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
# C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
# 1, 1.5, 2,5))
#svm.c <- train(y ~., data = train_smote,
#               method='svmRadial',
#               trControl=fitControl,
#               metric="ROC",tuneLength = 5)

#plot(svmModel)
#pred.svmModel <- as.vector(predict(svm.c, newdata=test, type="prob")[,"Yes"])
#roc.svmModel <- pROC::roc(test$y, pred.svmModel)
#auc.svmModel <- pROC::auc(roc.svmModel)

#Neural Network
#nnetGrid <- expand.grid(.size=1:10,
#                        .decay=c(0,.1,1,2))
#maxSize <- max(nnetGrid$.size)
#numWts <- 1*(maxSize*(length(dim(train)[2])) + maxSize + 1)
#set.seed(5346)
#nnetModel <- train(y ~ ., data=train_smote,
#                   method="nnet",
#                   metric="ROC",
#                   trControl=fitControl,
#                   tuneLength=5,
#                   trace=FALSE,
#                   maxit=2000)
#pred.nnetModel <- as.vector(predict(nnetModel, newdata=test, type="prob")[,"Yes"])
#roc.nnetModel <- pROC::roc(test$y, pred.nnetModel)
#auc.nnetModel <- pROC::auc(roc.nnetModel)

# Plotting model performance
test.auc_sm <- data.frame(model=c("glm","gbm","glmnet","rForest"),auc=c(auc.glmModel_sm, auc.gbmModel_sm, auc.eNetModel_sm, auc.rfModel_sm))
test.auc_sm <- test.auc_sm[order(test.auc_sm$auc, decreasing=TRUE),]
test.auc_sm$model <- factor(test.auc_sm$model, levels=test.auc_sm$model)
test.auc_sm
p_sm<-ggplot(data=test.auc_sm, aes(x=model, y=auc)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()
p_sm

grid.arrange(p, p_sm, ncol=2)

plot(roc.rfModel_sm, legacy.axes=TRUE, col="red")
lines(roc.gbmModel_sm, col="blue")
lines(roc.eNetModel_sm, col="green")

#feature importance plots
#RF
rf.varplot_sm <- cbind.data.frame(Predictors=rownames(data.frame(varImp(rfModel_sm,scale=TRUE)[1])),
                                  Importance=data.frame(varImp(rfModel_sm,scale=TRUE)[1])$Overall)
rf.varplot_sm <- rf.varplot_sm[order(rf.varplot_sm$Importance),]
dim(rf.varplot_sm)
rf.varplot2_sm <- rf.varplot_sm[51:70,]
rf.varplot2_sm$Predictors <- factor(rf.varplot2_sm$Predictors,levels=rf.varplot2_sm$Predictors)
(p <- ggplot(rf.varplot2_sm) + geom_point(aes(Importance,Predictors)) +ggtitle("Variable Importance Plot for Random Forest")) 

ggplot(rf.varplot2_sm, aes(x=Predictors, 
                           y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Random Forest") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()

#GBM
varImp(gbmModel_sm)
gbm.data_sm <- data.frame(varImp(gbmModel_sm)[1])
gbm.data2_sm <- cbind.data.frame(Predictors=rownames(gbm.data_sm),Importance=gbm.data_sm$Overall)
gbm.data2_sm <- gbm.data2_sm[order(gbm.data2_sm$Importance),]
gbm.data2_sm <- gbm.data2_sm[51:70,]
gbm.data2_sm <- gbm.data2_sm[order(gbm.data2_sm$Importance),]
gbm.data2_sm$Predictors <- factor(gbm.data2_sm$Predictors,levels=gbm.data2_sm$Predictors)
(p <- ggplot(gbm.data2_sm) + geom_point(aes(Importance,Predictors)) + ggtitle("Variable Importance Plot for GBM Model"))

ggplot(gbm.data2_sm, aes(x=Predictors, 
                         y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()

#ElasticNet
enetdata_sm <- varImp(eNetModel_sm)$importance
enetdata2_sm <- data.frame(Feature=rownames(enetdata_sm), Importance=enetdata_sm$Overall)
enetdata2_sm$Feature <- as.character(enetdata2_sm$Feature)
enetdata2_sm <- enetdata2_sm[order(enetdata2_sm$Importance),]
enetdata2_sm <- enetdata2_sm[51:70,]
enetdata2_sm$Feature <- factor(enetdata2_sm$Feature, levels=enetdata2_sm$Feature[order(enetdata2_sm$Importance)])

ggplot(enetdata2_sm, aes(x=Feature, 
                         y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="green") +
  coord_flip()

```


# Comparing Models

```{r include=FALSE}
#Without SMOTE

#finding common features between the models
intersect(rf.varplot2$Predictors,gbm.data2$Predictors)
commonfeats <- intersect(intersect(rf.varplot2$Predictors,gbm.data2$Predictors),enetdata2$Feature)

rf.varplot2$ranking <- rank(-rf.varplot2$Importance)
gbm.data2$ranking <- rank(-gbm.data2$Importance)
enetdata2$ranking <- rank(-enetdata2$Importance)

#rankings between rf, gbm, and elastic net
commonfeatsdf <- data.frame(Features = commonfeats,
                            RandomForest = rf.varplot2$ranking[rf.varplot2$Predictors %in% commonfeats],
                            RFImp = rf.varplot2$Importance[rf.varplot2$Predictors %in% commonfeats],
                            GBM = gbm.data2[order(match(gbm.data2$Predictors, commonfeats)),]$ranking[gbm.data2[order(match(gbm.data2$Predictors, commonfeats)),]$Predictors %in% commonfeats],
                            GBMImp = gbm.data2[order(match(gbm.data2$Predictors, commonfeats)),]$Importance[gbm.data2[order(match(gbm.data2$Predictors, commonfeats)),]$Predictors %in% commonfeats],
                            ElasticNet = enetdata2[order(match(enetdata2$Feature, commonfeats)),]$ranking[enetdata2[order(match(enetdata2$Feature, commonfeats)),]$Feature %in% commonfeats],
                            ENetImp = enetdata2[order(match(enetdata2$Feature, commonfeats)),]$Importance[enetdata2[order(match(enetdata2$Feature, commonfeats)),]$Feature %in% commonfeats]
)



#With SMOTE

#finding common features between the models
intersect(rf.varplot2_sm$Predictors,gbm.data2_sm$Predictors)
commonfeats_sm <- intersect(intersect(rf.varplot2_sm$Predictors,gbm.data2_sm$Predictors),enetdata2_sm$Feature)

rf.varplot2_sm$ranking <- rank(-rf.varplot2_sm$Importance)
gbm.data2_sm$ranking <- rank(-gbm.data2_sm$Importance)
enetdata2_sm$ranking <- rank(-enetdata2_sm$Importance)

#rankings between rf, gbm, and elastic net
commonfeatsdf_sm <- data.frame(Features = commonfeats_sm,
                               RandomForest = rf.varplot2_sm$ranking[rf.varplot2_sm$Predictors %in% commonfeats_sm],
                               RFImp = rf.varplot2_sm$Importance[rf.varplot2_sm$Predictors %in% commonfeats_sm],
                               GBM = gbm.data2_sm[order(match(gbm.data2_sm$Predictors, commonfeats_sm)),]$ranking[gbm.data2_sm[order(match(gbm.data2_sm$Predictors, commonfeats_sm)),]$Predictors %in% commonfeats_sm],
                               GBMImp = gbm.data2_sm[order(match(gbm.data2_sm$Predictors, commonfeats_sm)),]$Importance[gbm.data2_sm[order(match(gbm.data2_sm$Predictors, commonfeats_sm)),]$Predictors %in% commonfeats_sm],
                               ElasticNet = enetdata2_sm[order(match(enetdata2_sm$Feature, commonfeats_sm)),]$ranking[enetdata2_sm[order(match(enetdata2_sm$Feature, commonfeats_sm)),]$Feature %in% commonfeats_sm],
                               ENetImp = enetdata2_sm[order(match(enetdata2_sm$Feature, commonfeats_sm)),]$Importance[enetdata2_sm[order(match(enetdata2_sm$Feature, commonfeats_sm)),]$Feature %in% commonfeats_sm]
)


```




