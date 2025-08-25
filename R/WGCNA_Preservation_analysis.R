#Auther: Menna Arafat
# This code is based on the tutorial "Corrected R code from chapter 12 of the book"
# by Steve Horvath, available at https://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html
#set WD
getwd()
setwd("C:/Users/USER/Documents/rhabdo/")
dir.create("preservation_analysis")


#load packages
library(pROC)
library(ape)
library(ggdendro)
library(WGCNA)
library(stats)
library(flashClust)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(tidyverse)
library(gridExtra)
library(gplots)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
allowWGCNAThreads()          # allow multi-threading (optional)

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

#keep only good samples and genes
gsg = goodSamplesGenes(t(data), verbose = 3);
data= data[gsg$goodGenes, gsg$goodSamples]


#metadata
metadata= read.csv("input/metadata_binarized.csv")
#-----------------------------------------
#module preservation
#calculate module preservation between two expression sets, that must have equal number of samples and genes    
# number of networks used in the consensus network analysis:
nSets = 2
# Vector with descriptive names of the two sets.
setLabels = c("ctrl","emb")
shortLabels = c("C","B")
#genes must be in columns
table(metadata$condition)
ctrl= data[, grepl("Ctrl", names(data))] %>% dplyr::select(-"Ctrl18_3") %>% t(.)  %>% as.data.frame()
alv= data[ ,grepl("Alv", names(data))]  %>% t(.) %>% as.data.frame()
emb= data[ ,grepl("Emb", names(data))]%>% dplyr::select(-"Emb14_2") %>%  t(.)  %>% as.data.frame()


#data contain genes or samples
# with zero variance or excessive counts of missing entries.
# Please use the function goodSamplesGenes on each set to remove them
gsg = goodSamplesGenes(ctrl, verbose = 3);
ctrl= ctrl[ gsg$goodSamples, gsg$goodGenes]
#alv= alv[ gsg$goodSamples, gsg$goodGenes]
emb= emb[ gsg$goodSamples, gsg$goodGenes]
#calculate module preservation between two expression sets, that must have equal number of samples and genes    
#same genes between two data sets
ctrl= ctrl[, names(ctrl) %in% names(emb)]
alv= alv[, names(alv) %in% names(emb)]
emb= emb[, names(emb) %in% names(ctrl)]
dim(ctrl)
dim(alv)
dim(emb)


#build reference network
# Choose a set of soft thresholding powers
powers = c(1:10)  # in practice this should include powers up to 20.
# choose power based on SFT criterion
sft = pickSoftThreshold(ctrl, powerVector = powers)
# Plot the results:
dir.create("preservation_analysis/plots")
png("plots/dendrogram.png", width = 8000, height = 6000, res= 600) 
par(mfrow = c(1, 2))
# SFT index as a function of different powers
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", 
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.8, col = "red")
# Mean connectivity as a function of different powers
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

dev.off()
#automatic module detetion
mergingThresh = 0.25
net = blockwiseModules(ctrl, corType = "pearson", maxBlockSize = 5000, 
                       networkType = "signed", power = 4, minModuleSize =13, mergeCutHeight = .25, 
                       numericLabels = F, saveTOMs = TRUE, pamRespectsDendro = FALSE, saveTOMFileBase = "alveolarTOM")

moduleColorsref = net$colors #vector that assign every gene to its module
unique(moduleColorsref)
#in case of numeric labels is set TRUE
# Convert labels to colors for plotting
#moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)

# A data frame with module eigengenes can be obtained as follows
MEs = net$MEs

# Define a list whose components contain the data
# We now set up the multi-set expression data, genes in columns, samples in rows
# and corresponding module colors:
multiExpr=list(ctrl=list(data=ctrl),
               emb=list(data=emb))

# match label between the initial modules for whole data and modules calculated for reference group only
#so that the group of proteins under certain module for example brown module for whole data correspond to the same proteins in the brown module forreference group
moduleColorsAlldata= read.csv("module_colors_assign.csv")
moduleColorsAutomatic= moduleColorsAlldata[,2] %>% as.vector()
names(moduleColorsAutomatic)= moduleColorsAlldata[,1]
moduleColorsAutomatic
moduleColorsref 
moduleColorsAutomatic_subset=  moduleColorsAutomatic[names(moduleColorsAutomatic) %in% names(moduleColorsref)]
modLabCons = matchLabels(moduleColorsref, moduleColorsAutomatic_subset) %>% as.character()
names(modLabCons) = names( moduleColorsAutomatic_subset)
unique(modLabCons) 
#check that works out
names(modLabCons[modLabCons== "brown"])  %in%  names(moduleColorsAutomatic_subset[moduleColorsAutomatic_subset== "brown"])

#extract new corresponding modules
brown <- names(modLabCons[modLabCons== "brown"])
turquoise<- names(modLabCons[modLabCons== "turquoise"])
blue <- names(modLabCons[modLabCons== "blue"])
grey <- names(modLabCons[modLabCons== "grey"])

# set corresponding module colors for reference group  
multicolor= list(ctrl= modLabCons)
#multiColor=list(alv=moduleColorsAutomatic[ names(moduleColorsAutomatic)%in% names(as.data.frame(alv))  ])
# Check that the data has the correct format:
exprSize = checkSets(multiExpr)
exprSize
#OR
# multiExpr = vector(mode = "list", length = nSets)
# multiExpr[[1]] = list(data = before)
# names(multiExpr[[1]]$data) = names(before)
# rownames(multiExpr[[1]]$data) = dimnames(before)[[1]]
# multiExpr[[2]] = list(data = after2)
# names(multiExpr[[2]]$data) = names(after2)
# rownames(multiExpr[[2]]$data) = dimnames(after2)[[1]]
#Run module preservation function
set.seed(1)
system.time({
  mp = modulePreservation(multiExpr, multicolor,
                          referenceNetworks = 1, nPermutations = 10,
                          randomSeed = 1, quickCor = 0, verbose = 3)
})
# Save the results of the module preservation analysis
save(mp, file = "preservation_analysis/modulePreservation_ctrl_emb.RData")
# If needed, reload the data:
load(file = "preservation_analysis/modulePreservation_ctrl_emb.RData")

# specify the reference and the test networks
ref=1; test = 2

Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats= mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats
# Z statistics from the permutation test analysis
Z.PreservationStats



# plot median rank
png("preservation_analysis/plots/module_preserv_ctrl_emb.png", width = 3500, height = 3000, res= 600)
modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label = modColors[selectModules]

#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres

par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSize[selectModules],medianRank[selectModules],col=1,
     bg=modColors[selectModules],pch = 21,main="medianRank Preservation",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x", ylim= c(0,4.5))
#labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)

# plot Zsummary versus module size
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,main="Zsummary Preservation",
     cex=2,ylab ="Zsummary", xlab = "Module size", log = "x", ylim= c(0,10.5))
#labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
dev.off()
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#(consensus modules)
# the genes in both expression sets must be in alignment with each other, ie equal number and same identity
alv= alv[ ,names(as.data.frame(alv)) %in% names(as.data.frame(emb))]
emb= emb[ ,names(as.data.frame(emb)) %in% names(as.data.frame(alv))]
multiExpr=list(alv=list(data=alv),
               emb=list(data=emb))

netConsensus = blockwiseConsensusModules(multiExpr, maxBlockSize = 5000, power = 5, 
                                         minModuleSize = 50, deepSplit = 2, pamRespectsDendro = FALSE, mergeCutHeight = 0.25, 
                                         numericLabels = TRUE, minKMEtoStay = 0, saveTOMs = TRUE, checkMinModuleSize = FALSE)


# list of consensus module eigengenes
consMEs = netConsensus$multiMEs
# module labels
modLabCons0 = netConsensus$colors
# Relabel the consensus modules so that their labels match those from the
# automatic analysis in before group only (control group)
help(matchLabels)
moduleColorsAutomatic= moduleColorsAutomatic[names(moduleColorsAutomatic) %in% names(modLabCons0)]
modLabCons = matchLabels(moduleColorsAutomatic, modLabCons0)
# check the agreement between consensus- and females-only modules
mean(modLabCons == moduleLabelsAutomatic)
# Convert the numeric labels to color labels
moduleColorsConsensus = labels2colors(modLabCons)

blocknumber = 1
consTree = netConsensus$dendrograms[[blocknumber]]
# plot the consensus dendrogram and module colors
datColors = data.frame(moduleColorsConsensus, moduleColorsAutomatic)[netConsensus$blockGenes[[blocknumber]], 
]

plotDendroAndColors(consTree, datColors, c("Cons module", "alv_mods."), 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Consensus gene dendrogram and module colors")


