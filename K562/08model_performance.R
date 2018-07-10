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