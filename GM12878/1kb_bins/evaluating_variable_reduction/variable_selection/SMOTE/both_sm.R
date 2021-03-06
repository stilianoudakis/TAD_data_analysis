#Evaluating Variable Reduction Techniques using variable selection 

#loading packages 

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
library(DMwR)
library(gridExtra)
library(ggplot2)
library(leaps)

#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")

#center and scaling data to avoid using intercept term
cols <- names(Filter(is.numeric, chr1_gm12878_f))
chr1_gm12878_f[,cols] <- scale(chr1_gm12878_f[,cols], center = TRUE, scale = TRUE)


#Use smote to create balanced data 
chr1_gm12878_f_smote <- SMOTE(y ~ ., 
                              data=chr1_gm12878_f, 
                              perc.over = 100, 
                              perc.under = 200)

# Performing stepwise selection



#Using cross validation (10 fold)


#Both
k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_gm12878_f_smote), replace=TRUE)
cv.preds.both=matrix(NA, nrow=ncol(chr1_gm12878_f_smote)-1,ncol=k)
auc.model.both <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_gm12878_f_smote[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_gm12878_f_smote[folds!=j,], family = binomial)
  
  best.fit.both = step(glm.null, 
                      list(lower=formula(glm.null),
                           upper=formula(glm.full)),
                      direction="both",
                      trace=0)
  
  numpreds <- length(names(best.fit.both$coefficients)[-1])
  cv.preds.both[(1:numpreds),j] <- names(best.fit.both$coefficients)[-1]
  cols <- names(best.fit.both$model)
  model <- glm(y ~ . , data = chr1_gm12878_f_smote[folds==j,cols], family = binomial)
  pred.model <- predict(model, newdata=chr1_gm12878_f_smote[folds==j,cols], type="response")
  roc.model <- roc(chr1_gm12878_f_smote[folds==j,"y"], pred.model)
  auc.model.both[j] <- pROC::auc(roc.model)
  
}

saveRDS(cv.preds.both, "cv.preds.both.sm.rds")
saveRDS(auc.model.both, "auc.model.both.sm.rds")

