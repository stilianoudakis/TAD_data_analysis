#Variables selected

#Full
full.vars <- rownames(enetlst.full[[4]])
length(full.vars)
fwd.vars <- rownames(enetlst.fwd[[4]])
length(fwd.vars)
bwd.vars <- rownames(enetlst.bwd[[4]])
length(bwd.vars)
fwd.sm.vars <- rownames(enetlst.fwd.sm[[4]])
length(fwd.sm.vars)
bwd.sm.vars <- rownames(enetlst.bwd.sm[[4]])
length(bwd.sm.vars)
b.vars <- rownames(enetlst.b[[4]])
length(b.vars)
b.r.vars <- rownames(enetlst.b.r[[4]])
length(b.r.vars)



#Model performance


auc.fwd <- mean(enetlst.fwd[[3]])
round(auc.fwd,3)

auc.bwd <- mean(enetlst.bwd[[3]])
round(auc.bwd,3)

auc.fwd.sm <- mean(enetlst.fwd.sm[[3]])
round(auc.fwd.sm,3)

auc.bwd.sm <- mean(enetlst.bwd.sm[[3]])
round(auc.bwd.sm,3)

auc.b <- mean(enetlst.b[[3]])
round(auc.b,3)

auc.b.r <- mean(enetlst.b.r[[3]])
round(auc.b.r,3)

auc.full <- mean(enetlst.full[[3]])
round(auc.full,3)

test.auc <- data.frame("Variable Reduction"=c("Forward R.S.", 
                                       "Backward R.S.",
                                       "Forward SMOTE",
                                       "Backward SMOTE",
                                       "Boruta Full",
                                       "None"),
                       auc=c(auc.fwd,
                             auc.bwd,
                             auc.fwd.sm,
                             auc.bwd.sm,
                             auc.b,
                             auc.full))

test.auc <- test.auc[order(test.auc$auc, decreasing=TRUE),]

test.auc$Variable.Reduction <-factor(test.auc$Variable.Reduction, 
                                     levels=test.auc$Variable.Reduction)

p<-ggplot(data=test.auc, aes(x=Variable.Reduction, y=auc)) + 
  xlab("Variable Reduction") + ylab("AUC") +
  geom_bar(stat="identity", fill="steelblue") + ylim(0,1) +
  geom_text(aes(label=c(21,22,40,20,65,36)), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p


plot(rowMeans(enetlst.fwd[[2]]),
     rowMeans(enetlst.fwd[[1]]), 
     type="l", col="red",
     xlab="1-Specificity",
     ylab="Sensitivity")
lines(rowMeans(enetlst.bwd[[2]]),
      rowMeans(enetlst.bwd[[1]]),
      type="l", col="blue")
lines(rowMeans(enetlst.fwd.sm[[2]]),
      rowMeans(enetlst.fwd.sm[[1]]),
      type="l", col="green")
lines(rowMeans(enetlst.bwd.sm[[2]]),
      rowMeans(enetlst.bwd.sm[[1]]),
      type="l", col="orange")
lines(rowMeans(enetlst.b[[2]]),
      rowMeans(enetlst.b[[1]]),
      type="l", col="pink")
lines(rowMeans(enetlst.full[[2]]),
      rowMeans(enetlst.full[[1]]),
      type="l", col="black")
abline(a=0, b=1)
legend("bottomright", legend = c("Forward R.S.", 
                                 "Backward R.S.",
                                 "Forward SMOTE",
                                 "Backward SMOTE",
                                 "Boruta",
                                 "None"),
       fill=c("red","blue","green","orange","yellow","black"),
       cex=.75)


#comparing results

#foward
varimp.enet.fwd <- as.vector(rowMeans(enetlst.fwd[[4]]))
Labels <- rownames(enetlst.fwd[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.enet.df.fwd <- data.frame(Feature=Labels,
                                 Importance=varimp.enet.fwd)
varimp.enet.df.fwd <- varimp.enet.df.fwd[order(varimp.enet.df.fwd$Importance),]
varimp.enet.df.fwd$Feature <- factor(varimp.enet.df.fwd$Feature,
                                     levels=varimp.enet.df.fwd$Feature)
enetp.fwd <- ggplot(varimp.enet.df.fwd, aes(x=Feature, 
                                            y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="red") +
  coord_flip()


#backward
varimp.enet.bwd <- as.vector(rowMeans(enetlst.bwd[[4]]))
Labels <- rownames(enetlst.bwd[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.enet.df.bwd <- data.frame(Feature=Labels,
                                 Importance=varimp.enet.bwd)
varimp.enet.df.bwd <- varimp.enet.df.bwd[order(varimp.enet.df.bwd$Importance),]
varimp.enet.df.bwd$Feature <- factor(varimp.enet.df.bwd$Feature,
                                     levels=varimp.enet.df.bwd$Feature)
enetp.bwd <- ggplot(varimp.enet.df.bwd, aes(x=Feature, 
                                            y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="blue") +
  coord_flip()

#forward smote
varimp.enet.fwd.sm <- as.vector(rowMeans(enetlst.fwd.sm[[4]]))
Labels <- rownames(enetlst.fwd.sm[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.enet.df.fwd.sm <- data.frame(Feature=Labels,
                                 Importance=varimp.enet.fwd.sm)
varimp.enet.df.fwd.sm <- varimp.enet.df.fwd.sm[order(varimp.enet.df.fwd.sm$Importance),]
varimp.enet.df.fwd.sm$Feature <- factor(varimp.enet.df.fwd.sm$Feature,
                                     levels=varimp.enet.df.fwd.sm$Feature)
numvarenet <- dim(varimp.enet.df.fwd.sm)[1]
varimp.enet.df.fwd.sm <- varimp.enet.df.fwd.sm[(numvarenet-19):numvarenet,]
enetp.fwd.sm <- ggplot(varimp.enet.df.fwd.sm, aes(x=Feature, 
                                            y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="green") +
  coord_flip()


#backward smote
varimp.enet.bwd.sm <- as.vector(rowMeans(enetlst.bwd.sm[[4]]))
Labels <- rownames(enetlst.bwd.sm[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.enet.df.bwd.sm <- data.frame(Feature=Labels,
                                    Importance=varimp.enet.bwd.sm)
varimp.enet.df.bwd.sm <- varimp.enet.df.bwd.sm[order(varimp.enet.df.bwd.sm$Importance),]
varimp.enet.df.bwd.sm$Feature <- factor(varimp.enet.df.bwd.sm$Feature,
                                        levels=varimp.enet.df.bwd.sm$Feature)
numvarenet <- dim(varimp.enet.df.bwd.sm)[1]
varimp.enet.df.bwd.sm <- varimp.enet.df.bwd.sm[(numvarenet-19):numvarenet,]
enetp.bwd.sm <- ggplot(varimp.enet.df.bwd.sm, aes(x=Feature, 
                                                  y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="orange") +
  coord_flip()


#boruta full
varimp.enet.b <- as.vector(rowMeans(enetlst.b[[4]]))
Labels <- rownames(enetlst.b[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.enet.df.b <- data.frame(Feature=Labels,
                                    Importance=varimp.enet.b)
varimp.enet.df.b <- varimp.enet.df.b[order(varimp.enet.df.b$Importance),]
varimp.enet.df.b$Feature <- factor(varimp.enet.df.b$Feature,
                                        levels=varimp.enet.df.b$Feature)
numvarenet <- dim(varimp.enet.df.b)[1]
varimp.enet.df.b <- varimp.enet.df.b[(numvarenet-19):numvarenet,]
enetp.b <- ggplot(varimp.enet.df.b, aes(x=Feature, 
                                                  y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="yellow") +
  coord_flip()

#boruta reduced
#varimp.enet.b.r <- as.vector(rowMeans(enetlst.b.r[[4]]))
#Labels <- rownames(enetlst.b.r[[4]])
#Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
#varimp.enet.df.b.r <- data.frame(Feature=Labels,
#                               Importance=varimp.enet.b.r)
#varimp.enet.df.b.r <- varimp.enet.df.b.r[order(varimp.enet.df.b.r$Importance),]
#varimp.enet.df.b.r$Feature <- factor(varimp.enet.df.b.r$Feature,
#                                   levels=varimp.enet.df.b.r$Feature)
#numvarenet <- dim(varimp.enet.df.b.r)[1]
#varimp.enet.df.b.r <- varimp.enet.df.b.r[(numvarenet-19):numvarenet,]
#enetp.b.r <- ggplot(varimp.enet.df.b.r, aes(x=Feature, 
#                                        y=Importance)) +
#  xlab("Predictors") +
#  ylab("Importance") +
#  #ggtitle("Importance Plot for Gradient Boosting Machine") +
#  geom_bar(stat="identity", 
#           width=.5, 
#           position="dodge",
#           fill="pink") +
#  coord_flip()


#no variable reduction
varimp.enet.full <- as.vector(rowMeans(enetlst.full[[4]]))
Labels <- rownames(enetlst.full[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.enet.df.full <- data.frame(Feature=Labels,
                                 Importance=varimp.enet.full)
varimp.enet.df.full <- varimp.enet.df.full[order(varimp.enet.df.full$Importance),]
varimp.enet.df.full$Feature <- factor(varimp.enet.df.full$Feature,
                                     levels=varimp.enet.df.full$Feature)
numvarenet <- dim(varimp.enet.df.full)[1]
varimp.enet.df.full <- varimp.enet.df.full[(numvarenet-19):numvarenet,]
enetp.full <- ggplot(varimp.enet.df.full, aes(x=Feature, 
                                            y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="black") +
  coord_flip()


grid.arrange(enetp.fwd,enetp.bwd,
             enetp.fwd.sm,enetp.bwd.sm,
             enetp.b,enetp.full,ncol=2)





