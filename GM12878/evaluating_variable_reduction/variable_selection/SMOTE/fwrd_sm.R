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

#forward
k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_gm12878_f_smote), replace=TRUE)
cv.preds.fwd=matrix(NA, nrow=ncol(chr1_gm12878_f_smote)-1,ncol=k)
auc.model.fwd <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_gm12878_f_smote[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_gm12878_f_smote[folds!=j,], family = binomial)
  
  best.fit.fwd = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="forward",
                      trace=0)
  
    numpreds <- length(names(best.fit.fwd$coefficients)[-1])
    cv.preds.fwd[(1:numpreds),j] <- names(best.fit.fwd$coefficients)[-1]
    cols <- names(best.fit.fwd$model)
    model <- glm(y ~ . , data = chr1_gm12878_f_smote[folds==j,cols], family = binomial)
    pred.model <- predict(model, newdata=chr1_gm12878_f_smote[folds==j,cols], type="response")
    roc.model <- roc(chr1_gm12878_f_smote[folds==j,"y"], pred.model)
    auc.model.fwd[j] <- pROC::auc(roc.model)
    
}

saveRDS(cv.preds.fwd, "cv.preds.fwd.sm.rds")
saveRDS(auc.model.fwd, "auc.model.fwd.sm.rds")

auc.model.fwd <- readRDS("auc.model.fwd.sm.rds")
cv.preds.fwd <- readRDS("cv.preds.fwd.sm.rds")

vars.fwd <- na.omit(cv.preds.fwd[,which.max(auc.model.fwd)])
vars.fwd[grep("_dist",vars.fwd,invert = TRUE)] <- unlist(lapply(vars.fwd[grep("_dist",vars.fwd,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))

chr1_gm12878_fwd_sm <- chr1_gm12878_f[,which((names(chr1_gm12878_f) %in% vars.fwd) | names(chr1_gm12878_f)=="y")]

dim(chr1_gm12878_fwd_sm)
#247632     41

saveRDS(chr1_gm12878_fwd_sm, "chr1_gm12878_fwd_sm.rds")
