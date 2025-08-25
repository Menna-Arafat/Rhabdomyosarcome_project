#AUther: MENNa Arafat

getwd()
setwd( "C:/Users/USER/Documents/rhabdo")

#load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(biomaRt)
library(expm)
library(gplots)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library("MetaboAnalystR")

#load data
list.files(getwd())
data= read.csv("input/ALL_PROT-ACC-zero.csv" )[-1,]
#remove duplicte proteins
data= data[!duplicated(data$NAME),]
sum(duplicated(data$X))
#set row names
row.names(data)= NULL
data = data %>% column_to_rownames(var = "NAME")

#change type to numeric
data[]= lapply(data, as.numeric)
metadata= read.csv("input/metadata_binarized.csv")
#normality check 
#normality check  for the whole data
m <- rowMeans(data, na.rm = T) #get mean for every metabolite (row)
ggpubr::ggdensity(m)
#shapiro test
shapiro.test(m)
#normality check for every group separately
alv = rowMeans(data[,grepl("Alv",names(data))],na.rm = T)
emb = rowMeans(data[,grepl("Emb",names(data))],na.rm = T)
ctrl = rowMeans(data[,grepl("Ctrl",names(data))],na.rm = T)
df <- data.frame(alv = alv, emb = emb, ctrl=ctrl )
normalityDf <- reshape2::melt(df, measure.vars = c("alv", "emb", "ctrl"))
#plot the data distribution before transformation
ggpubr::ggdensity(normalityDf, 'value', color= "variable") +xlim(1.615e+06,1.165e+08)
summary(normalityDf$value)

#normality check  for the whole data after transformation
m <- rowMeans(scaledDf, na.rm = T) #get mean for every metabolite (row)
ggpubr::ggdensity(m)
summary(m)
#Downstream analysis
list.files(paste0(getwd(), "/input"))
#setwd("New folder/")
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet,   "input/Emb_ctrl.csv"   , "colu", "disc")
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "AutoNorm", "NULL", ratio=FALSE)

mSet<-PlotNormSummary(mSet, "metabolites_norm", "png", 300, width=NA)
mSet<-PlotSampleNormSummary(mSet, "sample_norm", "png", 300, width=NA)


# normality test
x <- mSet$dataSet$norm
#write.csv(t(x), "input/normalized_rawdata.csv", row.names = TRUE)
means <- data.frame(apply(x, 2, mean))
shapiro.test(means$apply.x..2..mean.)

png("densityplot_after_norm.png")
plot(density(means$apply.x..2..mean.) , main = "After_normalization")
dev.off()



library(ggpubr)
plot4= ggqqplot(means$apply.x..2..mean.)
ggsave("qq_plot_after_pareto_norm.png",plot4, dpi = 300, width = 9, height = 5.5)

library(nortest)
sf.test(means$apply.x..2..mean.)  #Shapiro-Francia test for normality

# 1- Fold Change Analysis
mSet<-FC.Anal(mSet, 2, 0, FALSE) #Condition/control
######## soudy
x <- mSet$analSet$fc
y <- data.frame(fc_log = x$fc.log)
z= data.frame(fc_log = x$fc.all)
yz= cbind(z,y)
colnames(yz)[1]= "FC"
write.csv(yz, "fc.csv", row.names = TRUE)


# 2- T Tests
mSet<-Ttests.Anal(mSet, T, 0.05, FALSE, TRUE, "fdr", FALSE)
tt= mSet$analSet$tt$sig.mat
write.csv(tt , "output/Mann_Whitny_Emb_ctrl.csv" , row.names = TRUE)
x <- read.csv("t_test_pareto.csv", row.names = 1)


x$IDs <- rownames(x)
rownames(x) <-NULL
x$idu <- as.numeric(row.names(x))
xall <- x[c(1:23),]

