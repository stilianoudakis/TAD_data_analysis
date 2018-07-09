# Variable selection using forward selection for k562 cell line

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

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

#Chromosome 1
chr1_k562_f <- readRDS("chr1_k562_f.rds")

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_k562_f$y=="No")
samps <- sample(which(chr1_k562_f$y=="No"),length(which(chr1_k562_f$y=="Yes")))
chr1_k562_f <- rbind.data.frame(chr1_k562_f[samps,],
                                chr1_k562_f[which(chr1_k562_f$y=="Yes"),])


# Performing stepwise selection

#center and scaling data to avoid using intercept term
cols <- names(Filter(is.numeric, chr1_k562_f))
chr1_k562_f[,cols] <- scale(chr1_k562_f[,cols], center = TRUE, scale = TRUE)



#Using cross validation (10 fold)

#forward
k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_k562_f), replace=TRUE)
cv.preds.fwd=matrix(NA, nrow=ncol(chr1_k562_f)-1,ncol=k)
auc.model.fwd <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_k562_f[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_k562_f[folds!=j,], family = binomial)
  
  best.fit.fwd = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="forward",
                      trace=0)
  
  numpreds <- length(names(best.fit.fwd$coefficients)[-1])
  cv.preds.fwd[(1:numpreds),j] <- names(best.fit.fwd$coefficients)[-1]
  cols <- names(best.fit.fwd$model)
  model <- glm(y ~ . , data = chr1_k562_f[folds==j,cols], family = binomial)
  pred.model <- predict(model, newdata=chr1_k562_f[folds==j,cols], type="response")
  roc.model <- roc(chr1_k562_f[folds==j,"y"], pred.model)
  auc.model.fwd[j] <- pROC::auc(roc.model)
  
}

saveRDS(cv.preds.fwd, "cv.preds.fwd.rds")
saveRDS(auc.model.fwd, "auc.model.fwd.rds")

auc.model.fwd <- readRDS("auc.model.fwd.rds")
cv.preds.fwd <- readRDS("cv.preds.fwd.rds")

vars.fwd <- na.omit(cv.preds.fwd[,which.max(auc.model.fwd)])
vars.fwd[grep("_dist",vars.fwd,invert = TRUE)] <- unlist(lapply(vars.fwd[grep("_dist",vars.fwd,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))

chr1_k562_fwd <- chr1_k562_f[,which((names(chr1_k562_f) %in% vars.fwd) | names(chr1_k562_f)=="y")]

dim(chr1_k562_fwd)
#247672     15

saveRDS(chr1_k562_fwd, "chr1_k562_fwd.rds")
