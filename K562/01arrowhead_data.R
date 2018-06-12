#Arrowhead Data

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/data")

arrowhead_k562 <- read.table("GSE63525_K562_Arrowhead_domainlist.txt", header=T)
dim(arrowhead_k562)
#5975   12

arrowhead_k562 <- arrowhead_k562[,1:3]

arrowhead_k562$Chromosome <- paste("chr",arrowhead_k562$Chromosome,sep = "")
table(arrowhead_k562$Chromosome)

arrowhead_k562num <- arrowhead_k562[-which(arrowhead_k562$Chromosome=="chrX"),]
arrowhead_k562x <- arrowhead_k562[which(arrowhead_k562$Chromosome=="chrX"),]

arrowhead_k562num2 <- arrowhead_k562num[order(as.numeric(substr(arrowhead_k562num$Chromosome,4,5)), arrowhead_k562num$Start),]
arrowhead_k562x2 <- arrowhead_k562x[order(arrowhead_k562x$Start),]

arrowhead_k562_data <- rbind.data.frame(arrowhead_k562num2,arrowhead_k562x2)

dim(arrowhead_k562_data)
#5975    3

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data_analysis/data")

write.table(arrowhead_k562_data, file = "arrowhead_k562_data.txt",sep="\t",row.names=FALSE)

