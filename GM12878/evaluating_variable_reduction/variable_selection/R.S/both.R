#Evaluating Variable Reduction Techniques using variable selection 

#loading packages 

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)
library(leaps)

#setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")
setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/variable_selection/R.S./")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")

chr1_gm12878_f$A <- as.numeric(chr1_gm12878_f$A)
chr1_gm12878_f$B <- as.numeric(chr1_gm12878_f$B)

#randomly sample to reduce dataset
set.seed(123)
zclass <- which(chr1_gm12878_f$y=="No")
samps <- sample(which(chr1_gm12878_f$y=="No"),length(which(chr1_gm12878_f$y=="Yes")))
chr1_gm12878_f <- rbind.data.frame(chr1_gm12878_f[samps,],
                                     chr1_gm12878_f[which(chr1_gm12878_f$y=="Yes"),])


# Performing stepwise selection

#center and scaling data to avoid using intercept term
cols <- names(Filter(is.numeric, chr1_gm12878_f))
chr1_gm12878_f[,cols] <- scale(chr1_gm12878_f[,cols], center = TRUE, scale = TRUE)



#Using cross validation (10 fold)


#Both
k = 10
set.seed(789)
folds = sample(1:k,nrow(chr1_gm12878_f), replace=TRUE)
cv.preds.both=matrix(NA, nrow=ncol(chr1_gm12878_f)-1,ncol=k)
auc.model.both <- numeric(k)

for(j in 1:k){
  #null model
  glm.null <- glm(y ~ 1, data = chr1_gm12878_f[folds!=j,], family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = chr1_gm12878_f[folds!=j,], family = binomial)
  
  best.fit.both = step(glm.null, 
                      list(lower=formula(glm.null),
                           upper=formula(glm.full)),
                      direction="both",
                      trace=0)
  
  numpreds <- length(names(best.fit.both$coefficients)[-1])
  cv.preds.both[(1:numpreds),j] <- names(best.fit.both$coefficients)[-1]
  cols <- names(best.fit.both$model)
  model <- glm(y ~ . , data = chr1_gm12878_f[folds==j,cols], family = binomial)
  pred.model <- predict(model, newdata=chr1_gm12878_f[folds==j,cols], type="response")
  roc.model <- roc(chr1_gm12878_f[folds==j,"y"], pred.model)
  auc.model.both[j] <- pROC::auc(roc.model)
  
}

saveRDS(cv.preds.both, "cv.preds.both.rds")
saveRDS(auc.model.both, "auc.model.both.rds")


auc.model.both <- readRDS("auc.model.both.rds")
cv.preds.both <- readRDS("cv.preds.both.rds")

vars.both <- na.omit(cv.preds.both[,which(order(auc.model.both)==1)])

chr1_gm12878_both <- chr1_gm12878_f[,which((names(chr1_gm12878_f) %in% vars.both) | 
                                            names(chr1_gm12878_f)=="y" | 
                                            names(chr1_gm12878_f)=="A" |
                                            names(chr1_gm12878_f)=="B")]

saveRDS(chr1_gm12878_both, "chr1_gm12878_both.rds")

