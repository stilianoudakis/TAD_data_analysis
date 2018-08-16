#Evaluating SMOTE

library(caret)
library(DMwR)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

chr1data_f <- readRDS("chr1data.rds")

#subsetting the data for chromosome 1
#randomly choosing 1000 rows (50 yes/950 no)
#randomly choosing 10 columns
zeroclass <- sample(which(chr1data_f$y=="No"),950)
oneclass <- sample(which(chr1data_f$y=="Yes"),50)
index <- c(oneclass,zeroclass)
index <- sample(index)
chr1data_f2 <- chr1data_f[index, c(1,sample(c(2:63),10))]
prop.table(table(chr1data_f2$y))


# Splitting the data
set.seed(5228)
inTrainingSet <- createDataPartition(chr1data_f2$y,p=.7,list=FALSE)
train <- chr1data_f2[inTrainingSet,]
test <- chr1data_f2[-inTrainingSet,]

#Using Smote

#oversampling n(100)% - n neighbors from the 5 nearest neighbors are 
#chosen and a new observation is created in the direction of each
#That is, for each case in the orginal data set belonging to the 
#minority class, perc.over/100 new examples of that class will be created
#undersampling n(100)% - the majority class is randomly sampled by ransomly
#removing observations until the minority class becomes some specified
#percentage (n) of the majority class

table(train$y)
#No Yes 
#665  35
prop.table(table(train$y))
#No  Yes 
#0.95 0.05 

#200/100
train_smote1 <- SMOTE(y ~ ., data=train, perc.over = 200, perc.under = 100)
table(train_smote1$y)
#No Yes 
#70 105

#Here, SMOTE created 70 more synthetic observations from the minority class
#35+70=105
#this is because the perc.over = 200 corresponds to creating 1 new observation
#for every existing minority observation
#Likewise the perc.under = 100 corresponds to randomly sampling n=70 observations
#from the majority class

#200/200
train_smote2 <- SMOTE(y ~ ., data=train, perc.over = 200, perc.under = 200)
table(train_smote2$y)
#No Yes 
#140 105

#This time perc.under = 200 corresponds to randomly sampling double 
#the number of created synthetic observations (70*2=140)

#100/100
train_smote3 <- SMOTE(y ~ ., data=train, perc.over = 100, perc.under = 100)
table(train_smote3$y)
#No Yes 
#35  70

#Here, perc.over = 100 corresponds to only creating 1 new synthetic observation 
#per original minority observation
#Since perc.under = 100, this means that only n=35 (number of new observations created)
#are randomly sampled from the majority class

#100/200
train_smote4 <- SMOTE(y ~ ., data=train, perc.over = 100, perc.under = 200)
table(train_smote4$y)
#No Yes 
#70  70

#This time since perc.under = 200, SMOTE has randomly sampled twice as many majority
#samples as the number of minority samples created (35*2=70)

#Note: the combination of perc.over = 100 and perc.under = 200 will always result 
#in perfectly balanced data
