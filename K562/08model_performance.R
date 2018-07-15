# Model Performance

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData/K562")

roc.glmModel <- readRDS("roc.glmModel.rds")
enetlst <- readRDS("k562enetlst.rds")
rflst <- readRDS("k562rflst.rds")
gbmlst <- readRDS("k562gbmlst.rds")

auc.glm <- pROC::auc(roc.glmModel)
auc.enet <- mean(enetlst[[3]])
auc.rf <- mean(rflst[[3]])
auc.gbm <- mean(gbmlst[[3]])

test.auc <- data.frame(model=c("GLM","ElasticNet","RForest","GBM"),
                       auc=c(auc.glm,
                             auc.enet, 
                             auc.rf, 
                             auc.gbm))

test.auc <- test.auc[order(test.auc$auc, decreasing=TRUE),]

test.auc$model <- factor(test.auc$model, levels=test.auc$model)

test.auc

p<-ggplot(data=test.auc, aes(x=model, y=auc)) + 
  xlab("Model") + ylab("AUC") +
  geom_bar(stat="identity", fill="steelblue") + ylim(0,1) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p


plot(roc.glmModel)
lines(1-rowMeans(enetlst[[2]]),
      rowMeans(enetlst[[1]]),
      type="l", col="blue")
lines(1-rowMeans(rflst[[2]]),
      rowMeans(rflst[[1]]),
      type="l", col="green")
lines(1-rowMeans(gbmlst[[2]]),
      rowMeans(gbmlst[[1]]),
      type="l", col="red")
legend("bottomright", legend = c("GLM", 
                                 "Elastic Net",
                                 "Random Forest",
                                 "GBM"),
       fill=c("black","red","blue","green"),
       cex=.75)


#Elastic Net
varimp.enet <- as.vector(rowMeans(enetlst[[4]]))
Labels <- rownames(enetlst[[4]])
Labels[grep("k562_", Labels)] <- gsub("k562_","",Labels[grep("k562_", Labels)])
varimp.enet.df <- data.frame(Feature=Labels,
                                 Importance=varimp.enet)
varimp.enet.df <- varimp.enet.df[order(varimp.enet.df$Importance),]
varimp.enet.df$Feature <- factor(varimp.enet.df$Feature,
                                     levels=varimp.enet.df$Feature)
#numvarenet <- dim(varimp.enet.df)[1]
#varimp.enet.df <- varimp.enet.df[(numvarenet-19):numvarenet,]
enetp <- ggplot(varimp.enet.df, aes(x=Feature, 
                                            y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()


#Random Forest
varimp.rf <- as.vector(rowMeans(rflst[[4]]))
Labels <- rownames(rflst[[4]])
Labels[grep("k562_", Labels)] <- gsub("k562_","",Labels[grep("k562_", Labels)])
varimp.rf.df <- data.frame(Feature=Labels,
                             Importance=varimp.rf)
varimp.rf.df <- varimp.rf.df[order(varimp.rf.df$Importance),]
varimp.rf.df$Feature <- factor(varimp.rf.df$Feature,
                                 levels=varimp.rf.df$Feature)
#numvarrf <- dim(varimp.rf.df)[1]
#varimp.rf.df <- varimp.rf.df[(numvarrf-19):numvarrf,]
rfp <- ggplot(varimp.rf.df, aes(x=Feature, 
                                    y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()


#GBM
varimp.gbm <- as.vector(rowMeans(gbmlst[[4]]))
Labels <- rownames(gbmlst[[4]])
Labels[grep("k562_", Labels)] <- gsub("k562_","",Labels[grep("k562_", Labels)])
varimp.gbm.df <- data.frame(Feature=Labels,
                           Importance=varimp.gbm)
varimp.gbm.df <- varimp.gbm.df[order(varimp.gbm.df$Importance),]
varimp.gbm.df$Feature <- factor(varimp.gbm.df$Feature,
                               levels=varimp.gbm.df$Feature)
#numvargbm <- dim(varimp.gbm.df)[1]
#varimp.gbm.df <- varimp.gbm.df[(numvargbm-19):numvargbm,]
gbmp <- ggplot(varimp.gbm.df, aes(x=Feature, 
                                y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="green") +
  coord_flip()


#Comparing Results
#finding common features between the models

x <- intersect(varimp.rf.df$Feature,varimp.gbm.df$Feature)
z <- intersect(x,varimp.enet.df$Feature)
#z <- intersect(y,varimp.svm.df$Feature)
#z


varimp.rf.df$ranking <- rank(-varimp.rf.df$Importance)
varimp.gbm.df$ranking <- rank(-varimp.gbm.df$Importance)
varimp.enet.df$ranking <- rank(-varimp.enet.df$Importance)
#varimp.svm.df$ranking <- rank(-varimp.svm.df$Importance)


#rankings between rf, gbm, svm, and elastic net
commonfeatsdf <- data.frame(Features = z,
                            RandomForest = varimp.rf.df$ranking[varimp.rf.df$Feature %in% z],
                            RFImp = varimp.rf.df$Importance[varimp.rf.df$Feature %in% z],
                            GBM = varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$ranking[varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$Feature %in% z],
                            GBMImp = varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$Importance[varimp.gbm.df[order(match(varimp.gbm.df$Feature, z)),]$Feature %in% z],
                            ElasticNet = varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$ranking[varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$Feature %in% z],
                            ENetImp = varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$Importance[varimp.enet.df[order(match(varimp.enet.df$Feature, z)),]$Feature %in% z]
)

commonfeatsdf

#saveRDS(commonfeatsdf, "/home/stilianoudakisc/TAD_data_analysis/output/commonfeats_nosmote")

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/common_feats_nosmote")
datatable(commonfeatsdf)
#dev.off()

