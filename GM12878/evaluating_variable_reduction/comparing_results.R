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
                                       "Boruta Reduced",
                                       "None"),
                       auc=c(auc.fwd,
                             auc.bwd,
                             auc.fwd.sm,
                             auc.bwd.sm,
                             auc.b,
                             auc.b.r,
                             auc.full))

test.auc <- test.auc[order(test.auc$auc, decreasing=TRUE),]

test.auc$Variable.Reduction <-factor(test.auc$Variable.Reduction, 
                                     levels=test.auc$Variable.Reduction)

p<-ggplot(data=test.auc, aes(x=Variable.Reduction, y=auc)) + 
  xlab("Variable Reduction") + ylab("AUC") +
  geom_bar(stat="identity", fill="steelblue") + ylim(0,1) +
  theme_minimal()
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
lines(rowMeans(enetlst.bwd[[2]]),
      rowMeans(enetlst.bwd[[1]]),
      type="l", col="yellow")
lines(rowMeans(enetlst.b[[2]]),
      rowMeans(enetlst.b[[1]]),
      type="l", col="pink")
lines(rowMeans(enetlst.b.r[[2]]),
      rowMeans(enetlst.b.r[[1]]),
      type="l", col="black")
abline(a=0, b=1)
legend("bottomright", legend = c("Forward R.S.", 
                                 "Backward R.S.",
                                 "Forward SMOTE",
                                 "Backward SMOTE",
                                 "Boruta Full",
                                 "Boruta Reduced",
                                 "None"),
       fill=c("red","blue","green","orange","yellow","pink","black"),
       cex=.75)


#comparing results

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




grid.arrange(enetp.fwd,enetp.bwd,ncol=2)



