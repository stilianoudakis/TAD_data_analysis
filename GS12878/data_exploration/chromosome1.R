#Data Exploration

library(ggplot2)
library(DT)

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/RData")

#Chromosome 1
chr1_gm12878_f <- readRDS("chr1_gm12878_f.rds")

#cols <- c(grep("dist",colnames(chr1_gm12878_f)))
#chr1_gm12878_f[,cols] <- apply(chr1_gm12878_f[,cols], 2, function(x){2^x})

#Full data
#gm12878_f <- readRDS("gm12878_f.rds")

#looking into summaries of the continous variables

#distance variables
chr1_gm12878_f_d <- chr1_gm12878_f[,c(1,grep("dist",colnames(chr1_gm12878_f)))]

#randomly sample from the majority class
set.seed(123)
sampids <- sample(which(chr1_gm12878_f_d$y=="No"),
                  length(which(chr1_gm12878_f_d$y=="Yes")))

chr1_gm12878_f_d_r <- rbind.data.frame(chr1_gm12878_f_d[sampids,],
                                       chr1_gm12878_f_d[which(chr1_gm12878_f_d$y=="Yes"),])




for(i in 2:dim(chr1_gm12878_f_d_r)[2]){
p <- ggplot(chr1_gm12878_f_d_r,aes(x=chr1_gm12878_f_d_r[,i],group=y,fill=y))+
  xlab(colnames(chr1_gm12878_f_d_r)[i])+
  geom_histogram(position="identity",alpha=0.5,binwidth=0.25)+
  theme_bw()
print(p)
}

#summaries
meanmat <- matrix(nrow=dim(chr1_gm12878_f_d_r)[2]-1,ncol = 2)
medianmat <- matrix(nrow=dim(chr1_gm12878_f_d_r)[2]-1,ncol = 2)
rangemat <- matrix(nrow=dim(chr1_gm12878_f_d_r)[2]-1,ncol = 2)
pvalmat <- matrix(nrow=dim(chr1_gm12878_f_d_r)[2]-1,ncol = 1)

for(i in 1:(dim(chr1_gm12878_f_d_r)[2]-1)){
meanmat[i,] <- round(tapply(chr1_gm12878_f_d_r[,i+1],chr1_gm12878_f_d_r$y,mean),2)
medianmat[i,] <- round(tapply(chr1_gm12878_f_d_r[,i+1],chr1_gm12878_f_d_r$y,median),2)
rangemat[i,] <- round(tapply(chr1_gm12878_f_d_r[,i+1],chr1_gm12878_f_d_r$y,function(x) max(x)-min(x)),2)
pvalmat[i,] <- wilcox.test(chr1_gm12878_f_d_r[which(chr1_gm12878_f_d_r$y=="Yes"),i+1],
                  chr1_gm12878_f_d_r[which(chr1_gm12878_f_d_r$y=="No"),i+1])$p.value
}

summarydf <- cbind.data.frame(Feature=colnames(chr1_gm12878_f_d_r)[-1],
                              Mean=meanmat,
                              Median=medianmat,
                              Range=rangemat,
                              PValue=pvalmat)


summarydf <- summarydf[order(summarydf$PValue,decreasing = FALSE),]

datatable(summarydf)

############################################################################3

#binary variables
chr1_gm12878_f_b <- chr1_gm12878_f[,c(grep("dist",colnames(chr1_gm12878_f),invert = TRUE))]

#randomly sample from the majority class
set.seed(123)
sampids <- sample(which(chr1_gm12878_f_b$y=="No"),
                  length(which(chr1_gm12878_f_b$y=="Yes")))

chr1_gm12878_f_b_r <- rbind.data.frame(chr1_gm12878_f_b[sampids,],
                                       chr1_gm12878_f_b[which(chr1_gm12878_f_b$y=="Yes"),])




for(i in 2:dim(chr1_gm12878_f_d_r)[2]){
  p <- ggplot(chr1_gm12878_f_b_r,aes(x=chr1_gm12878_f_b_r[,i],group=y,fill=y))+
         xlab(colnames(chr1_gm12878_f_b_r)[i])+
         geom_bar(position="fill")+
         theme_bw()
  print(p)
}

tab <- list()

for(i in 1:(dim(chr1_gm12878_f_b_r)[2]-1)){
tab[[i]] <- round(prop.table(table(chr1_gm12878_f_b_r[,i+1],chr1_gm12878_f_b_r$y),margin = 2),2)
}
tab 

tab <- rbind(tab[[1]],tab[[2]],tab[[3]],tab[[4]],tab[[5]],tab[[6]],tab[[7]],tab[[8]],tab[[9]],tab[[10]],tab[[11]],tab[[12]],tab[[13]],tab[[14]],tab[[15]],tab[[16]],tab[[17]],tab[[18]])

pval = numeric()
for(i in 1:(dim(chr1_gm12878_f_b_r)[2]-1)){
pval[i] <- chisq.test(table(chr1_gm12878_f_b_r[,i+1],chr1_gm12878_f_b_r$y))$p.value
}
pvaldf <- data.frame(Feature = colnames(chr1_gm12878_f_b_r)[-1], Pvalue=pval)
pvaldf <- pvaldf[order(pvaldf$Pvalue, decreasing = FALSE),]

datatable(pvaldf)

