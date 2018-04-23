
setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data_analysis")

#Performing logistic regression
tad_data <- data.frame(class=y, A=X_A, B=X_B)
dim(tad_data)
#703900      3

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
