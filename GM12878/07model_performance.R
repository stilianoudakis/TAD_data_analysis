#Model Performance

#plotting performance

#mourad model
setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData/GM12878/mourad_model")

binMat.mati <- readRDS("binMat.mati.rds")
Analysisi <- readRDS("Analysisi.rds")

prob=predict(Analysisi,type=c("response"))
binMat.mati$prob=prob

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}
prob=predict(Analysisi,type=c("response"))
head(prob)
roc <- simple_roc(binMat.mati$Border, binMat.mati$prob)
head(roc)

auc.glmModel <- 0.7884

test.auc <- data.frame(model=c("GLM","ElasticNet","RForest","GBM"),
                       auc=c(auc.glmModel,
                             mean(auclst[[1]]), 
                             mean(auclst[[2]]), 
                             mean(auclst[[3]])))

test.auc <- test.auc[order(test.auc$auc, decreasing=TRUE),]

test.auc$model <- factor(test.auc$model, levels=test.auc$model)

test.auc

perf_nosmote<-ggplot(data=test.auc, aes(x=model, y=auc)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal()

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/perf_nosmote")
perf_nosmote
#dev.off()


#Variable Importance Plots

#ENET
varimp.enet <- as.vector(rowMeans(performlst[[1]]))
varimp.enet.df <- data.frame(Feature=rownames(performlst[[1]]),
                             Importance=varimp.enet)
varimp.enet.df <- varimp.enet.df[order(varimp.enet.df$Importance),]
numvarenet <- dim(varimp.enet.df)[1]
varimp.enet.df <- varimp.enet.df[(numvarenet-19):numvarenet,]
varimp.enet.df$Feature <- factor(varimp.enet.df$Feature,levels=varimp.enet.df$Feature)
enetp <- ggplot(varimp.enet.df, aes(x=Feature, 
                                    y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="green") +
  coord_flip()
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/enet_varimp_nosmote")
enetp
#dev.off()


#RF
varimp.rf <- as.vector(rowMeans(performlst[[2]]))
varimp.rf.df <- data.frame(Feature=rownames(performlst[[2]]),
                           Importance=varimp.rf)
varimp.rf.df <- varimp.rf.df[order(varimp.rf.df$Importance),]
numvarrf <- dim(varimp.rf.df)[1]
varimp.rf.df <- varimp.rf.df[(numvarrf-19):numvarrf,]
varimp.rf.df$Feature <- factor(varimp.rf.df$Feature,levels=varimp.rf.df$Feature)
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
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/rf_varimp_nosmote")
rfp
#dev.off()

#GBM
varimp.gbm <- as.vector(rowMeans(performlst[[3]]))
varimp.gbm.df <- data.frame(Feature=rownames(performlst[[3]]),
                            Importance=varimp.gbm)
varimp.gbm.df <- varimp.gbm.df[order(varimp.gbm.df$Importance),]
numvargbm <- dim(varimp.gbm.df)[1]
varimp.gbm.df <- varimp.gbm.df[(numvargbm-19):numvargbm,]
varimp.gbm.df$Feature <- factor(varimp.gbm.df$Feature,levels=varimp.gbm.df$Feature)
gbmp <- ggplot(varimp.gbm.df, aes(x=Feature, 
                                  y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/gbm_varimp_nosmote")
gbmp
#dev.off()

#SVM
#varimp.svm <- as.vector(rowMeans(performlst[[4]]))
#varimp.svm.df <- data.frame(Feature=rownames(performlst[[4]]),
#                            Importance=varimp.svm)
#varimp.svm.df <- varimp.svm.df[order(varimp.svm.df$Importance),]
#numvarsvm <- dim(varimp.svm.df)[1]
#varimp.svm.df <- varimp.svm.df[(numvarsvm-19):numvarsvm,]
#varimp.svm.df$Feature <- factor(varimp.svm.df$Feature,levels=varimp.svm.df$Feature)
#svmp <- ggplot(varimp.svm.df, aes(x=Feature, 
#                                  y=Importance)) +
#  xlab("Predictors") +
#  ylab("Importance") +
#  #ggtitle("Importance Plot for Gradient Boosting Machine") +
#  geom_bar(stat="identity", 
#           width=.5, 
#           position="dodge",
#           fill="red") +
#  coord_flip()
#theme(axis.text.x = element_text(angle = 90, hjust = 1))

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/varimps_nosmote")
#svmp
#dev.off()

#jpeg("/home/stilianoudakisc/TAD_data_analysis/output/svm_varimp_nosmote")
grid.arrange(enetp,rfp,gbmp, ncol=3)
#dev.off()



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

