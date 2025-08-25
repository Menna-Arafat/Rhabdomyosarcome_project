#Auther: Menna Arafat
# This code is based on the tutorial "Corrected R code from chapter 12 of the book"
# by Steve Horvath, available at https://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html


#-------------------------------------------------
#set workind directory
getwd()
setwd("C:/Users/USER/Documents/rhabdo_metabolomics")

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
list.files(paste0(getwd(), "/input"))
data= read.csv("input/ALL profile metabolomics-IDwithout18alv-19-20.csv" )[-1,] 

#remove duplicte proteins
data= data[!duplicated(data[1]),]
sum(duplicated(data$X))
#set row names
row.names(data)= NULL
data = data %>% column_to_rownames(var = names(data[1]))
datExpr= t(data)
goods <- goodSamplesGenes(datExpr, verbose = 3)
datExpr= datExpr[goods$goodSamples== TRUE, goods$goodGenes == TRUE ]

#change type to numeric
data[]= lapply(data, as.numeric)
#metadata
metadata= data.frame(Sample= colnames(data) ,
                     condition= gsub("[0-9_]+", "", colnames(data))) %>%
  column_to_rownames(var = "Sample")

metadata$condition= factor(metadata$condition, levels = c("Ctrl", "Alv", "Emb"))
design= model.matrix(~ 0+condition , metadata)
design


#metadata= read.csv("input/metadata_binarized.csv" )
#-------------------------------------------
#tree and modules and heatmap of associated traits/ phenotypes
#Hierarchical clustering of samples, detect outlier samples,and association of sample with certain trait
#with heatmap of such trait where red indicate high value

#Build adjacency matrix for samples
A = adjacency(data, type = "distance")
# this calculates the whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized connectivity
Z.k = scale(k)
# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -5  # often -2.5

# the color vector indicates outlyingness (red)
outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")

# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1 - A), method = "average")
# Convert traits to a color representation: where red indicates high
# values
traitColors = data.frame(numbers2colors(as.numeric(metadata$cond_binary), signed = TRUE))
#dimnames(traitColors)[[2]] = "Inflammation_lvl"
datColors = data.frame(outlier_Samples = outlierColor, Condition= traitColors)
colnames(datColors)[2]= "Condition"
# Plot the sample dendrogram and the colors underneath.
png("plots/WGCNA_dendrogram.png", width = 8000, height = 6000, res= 600)
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors, cex.rowText = 5,
                    main = "Sample dendrogram and Homogeneity of samples heatmap")

#grid.arrange(a1, a2, nrow = 2)
dev.off()

#-----------------------------------------------------------------
# Plot the dendrogram and the module colors before and after merging underneath
png("plots/dendrogram_merged_unmerged_modules.png", width = 2200, height = 2500, res= 600)
plotDendroAndColors(net$dendrograms[[1]], cbind(net$unmergedColors, net$colors),
                    c(paste0("Unmerged\n","Modules"), paste0("Merged\n", "Modules")),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    cex.colorLabels = 0.6,
                    guideHang = 0.05)
dev.off()
#-----------------------------------------------------------
#  TOM plot/ heatmap of modules for all proteins
datExpr= t(data)
dissTOM= 1 - TOMsimilarityFromExpr(datExpr, power= 5)
dendro= net$dendrograms[[1]]
moduleColorsAutomatic= net$colors

#visualizations
png("plots/TOM_PLOT_module_heatmap_proteins_all.png", width = 800, height = 600)
#myheatcol = colorpanel(250,'gold',"orange",'darkred')
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
# Transform dissTOM with a power to enhance visibility
TOMplot(dissTOM, dendro, moduleColorsAutomatic,col= myheatcol, 
        main = "Module Heatmap Plot, All Proteins")
dev.off()
#-------------------------------------------------------------
#  TOM plot/ heatmap of modules for selected proteins
#calculate fold change for proteins 
#Downstream analysis
library("MetaboAnalystR")
list.files(getwd())

#setwd("New folder/")
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "ALL_PROT-ACC-zero.csv" , "colu", "disc")
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE)
# normality test
x <- mSet$dataSet$norm
write.csv(x, "normalized_rawdata.csv", row.names = TRUE)
means <- data.frame(apply(x, 2, mean))
shapiro.test(means$apply.x..2..mean.)

# 1- Fold Change Analysis
mSet<-FC.Anal(mSet, 2, 0, FALSE) #condition/control
######## soudy
x <- mSet$analSet$fc
y <- data.frame(fc_log = x$fc.log)
z= data.frame(FC = x$fc.all)
yz= cbind(z,y)
write.csv(yz, "fc_alv_emb.csv", row.names = TRUE)
#----------
mSet<-Ttests.Anal(mSet, T, 0.05, FALSE, TRUE, "fdr", FALSE)
tt= mSet$analSet$tt$sig.mat
write.csv(tt , "output/wilcox_alv_emb.csv" , row.names = TRUE)


#subset fold change based on threshold 1.5, -1.5
fc_alv= read.csv("fc_alv_ctrl.csv") %>% filter(FC >= 1.5 | FC<= -1.5)
fc_emb= read.csv("fc_emb_ctrl.csv") %>% filter(FC >= 1.5 | FC<= -1.5)
alv_emb= read.csv("fc_alv_emb.csv") %>% filter(FC >= 1.5 | FC<= -1.5)

proteins= c(fc_alv$X , fc_emb$x, alv_emb$X) %>% unique()
# subset data to have only selected proteins
datExpr_subset= datExpr[,colnames(datExpr) %in% proteins]
dissTOM_subset= 1 - TOMsimilarityFromExpr(datExpr_subset, power= 5)
dendro_subset = hclust(as.dist(dissTOM_subset), method = "average")
moduleColors_subset= moduleColorsAutomatic[names(moduleColorsAutomatic)%in% proteins]

png("plots/module_heatmap_TOM_PLOT_selected.png", width = 2800, height = 3300, res= 600)
#myheatcol = colorpanel(250,'gold',"orange",'darkred')
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
# Transform dissTOM with a power to enhance visibility
TOMplot(dissTOM_subset, 
        dendro_subset, 
        moduleColors_subset,col= myheatcol,
        main = "Module Heatmap Plot, Selected Proteins")
dev.off()

#--------------------------------------------------------------------
#Module trait correlation
# Next use a single trait/variable or the whole metadata binarized to define the module significance 
#what module associated to what phenotype
traits= design[,c(2,3,1)]
head(traits)
# Define numbers of genes and samples
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)

module_eigengenes= read.csv("output/module_eigengenes.csv") %>% 
  column_to_rownames("X")

module.trait.corr <- WGCNA::cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

#module_trait heatmap of WGCNA package
# correlations and their p-values
png("plots/heatmap_module_trait_WGCNA_new.png", width = 3500, height = 4000, res= 600) 
textMatrix = paste(signif(module.trait.corr, 2), "\n(", signif(module.trait.corr.pvals, 1), ")", 
                   sep = "")
dim(textMatrix) = dim(module.trait.corr)
par(mar = c(6, 6, 4, 6))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.corr, xLabels = c( "Alv", "Emb", "Ctrl"), 
               yLabels = names(module_eigengenes), 
               ySymbols = names(module_eigengenes), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = T, cex.text = 0.8,
               zlim = c(-1, 1),xColorWidth = 1 * strheight("M"),
               yColorWidth = 1.5 * strwidth("M"),xColorOffset = strheight("M")/6, 
               yColorOffset = strwidth("M")/6, font.lab.x = 2, cex.legendLabel = 2,
               font.lab.y = 2, xLabelsAngle = 60, main = paste("Module-Condition Relationship"), plotLegend= TRUE)

dev.off()
#------------------------------------------------------------------
#MDS plot
png("plots/MDS_plot.png",  width = 2800, height = 3300, res= 600) 
dissTOM= 1 - TOMsimilarityFromExpr(datExpr, power= 5)
cmd1=cmdscale(as.dist(dissTOM),2)
par(mfrow=c(1,1))
plot(cmd1,col=moduleColorsAutomatic,main="MDS plot",
     xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")
dev.off()

#-----------------------------------------------------------
#Bar plot for phenotype discrimination by module eigengene
list.files(paste0(getwd(), "/output"))
module_eigengenes= read.csv("output/module_eigengenes.csv") %>% 
                        column_to_rownames("X") %>% lapply(., as.vector)
metadata= read.csv("input/metadata_binarized.csv" )
#Module eigengenes of WGCNA to statistically distinguish between the phenotypes of my data via Kruskal-Wallis test
shapiro.test(module_eigengenes$MEturquoise)

mods= names(module_eigengenes)
pvalues1= numeric(length = length(mods))

#check the differentiation between 2 groups (cancer/control) (non parametric) >> mann whitny
for (i in seq_along(mods)){
  wilcox_res1 <- wilcox.test(module_eigengenes[[i]] ~ metadata$cancer_status , paired= F)
  pvalues1[i] <-  wilcox_res1$p.value
}
fdr1 <- p.adjust(pvalues1, method = "BH") # Benjamini
names(fdr1)= mods
fdr1
res1= data.frame(pvalue= pvalues1,
                    fdr= fdr1)

res1$sig= ifelse(res1$fdr <= .05, "***", "")
res1 = res1%>% rownames_to_column("modules")
res1$modules= gsub("ME", "",res1$modules )
write.csv(res1, "output/module_sig_distinguish_cancer_control.csv")
#check the differentiation between 2 groups alveolar/ embryonal  (non parametric) >> wilcoxin
module_eigengenes= read.csv("output/module_eigengenes.csv") %>% 
                               filter(grepl("Alv|Emb" ,X))  %>%
                               column_to_rownames("X")  %>%
                               lapply(., as.vector)
mods= names(module_eigengenes)
pvalues2= numeric(length = length(mods))
for (i in seq_along(mods)){
  wilcox_res <- wilcox.test(module_eigengenes[[i]] ~ metadata[1:35, "con"], paired= F) #2 numeric vectors of moduleigengene for two groups
  pvalues2[i] <- wilcox_res$p.value
}
fdr2 <- p.adjust(pvalues2, method = "BH") # Benjamini
names(fdr2)= mods
fdr2
res2= data.frame(pvalue= pvalues2,
                    fdr= fdr2)
res2
res2$sig= ifelse(res2$fdr <= .05, "***", "")
res2 = res2%>% rownames_to_column("modules")
res2$modules= gsub("ME", "",res2$modules )
 write.csv(res2, "output/module_sig_distinguish_Alv_Emb.csv")


#visualization
p1=ggplot(res1) +
  geom_bar(aes(x = modules, y = -log10(fdr), fill = modules), stat = "identity") +
  coord_flip() + # Flip coordinates to make the bar plot horizontal
  scale_fill_identity() + # Use actual colors specified in the data frame
  geom_text(aes(x =  modules  , y =-log10(fdr) , label = sig), 
            position = position_dodge(width = 0.9), vjust = -0.5, size=.5 ,color = "black") +
  theme_minimal() + 
  labs(y = "-log10(FDR)", x = "Modules") +
  ggtitle("RMS/Ctrl Discrimination by Module Eigengene") + # Label axes and provide a title
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size=12),
        axis.title=element_text(size=13) ) + 
  ylim(0,12) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd= 1)

p1
-log10(res1$fdr)
#+
# annotate("text", x = .98, y = Inf, label = paste0("-------", "\n", "Sig. Cut-off"),
#   hjust = 1.05, vjust = 0, color = "black", size = 4.5, angle = 0)
p2=ggplot(res2) +
  geom_bar(aes(x = modules, y = -log10(fdr), fill = modules), stat = "identity") +
  coord_flip() + # Flip coordinates to make the bar plot horizontal
  scale_fill_identity() + # Use actual colors specified in the data frame
  geom_text(aes(x =  modules  , y =-log10(fdr) , label = sig), 
            position = position_dodge(width = 0.9), vjust = -0.5, size=.5, color = "black") +
  theme_minimal() + # Use a minimal theme for the plot
  labs(y = "-log10(FDR)", x = "Modules") +
  ggtitle("Alv/Emb RMS Discrimination by Module Eigengene") + # Label axes and provide a title
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size=12),
        axis.title=element_text(size=13) ) + 
  ylim(0, 10) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd= 1)

p2
p= grid.arrange(p1, p2, ncol = 2)

ggsave("plots/modules_discrimination_new.png", p,width= 10, height = 7, dpi=900 )
#dev.off()
#---------------------------------------------------------
#chord plot of key drivers
list.files(paste0(getwd(), "/output"))
hubs_alv= read.csv("output/hubproteins_gene_names_Alveolar.csv" ) %>% unlist(.) %>% .[. !="" & !is.na(.)] 
hubs_emb= read.csv("output/hubproteins_gene_names_Embrynoal.csv" ) %>% unlist(.) %>% .[. !="" & !is.na(.)] 
hubs_ctrl= read.csv("output/hubproteins_gene_names_Ctrl.csv" ) %>% unlist(.) %>% .[. !="" & !is.na(.)]
keydrivers= c(hubs_alv, hubs_emb, hubs_ctrl) %>% unique(.)


mtx= matrix(nrow= 3, ncol = length( keydrivers)) #number of all key drivers 
row.names(mtx)= c( paste(" Alveolar" ,"\n", "RMS"), paste("  Embryonal", "\n", "RMS"), paste(" Control"))
colnames(mtx)= paste0(keydrivers)

#build similarity matrix 
#colnames of matrix included in keydrivers specified for certain module/ phenotype(rows) then put in 1
mod= list(alv = hubs_alv ,emb=hubs_emb, Ctrl=hubs_ctrl  )
mod[[1]]
for (i in seq_along(mod)){
  for (j in 1:ncol(mtx)) {
    if ( colnames(mtx)[j] %in% mod[[i]] ) {
      mtx[i, j] <- 1
    } else {
      mtx[i, j] <- 0
    }
  }
}

library(circlize)

png("plots/chord_plot_hubproteins.png", width = 4700, height = 4700, res= 600)

par(cex = .8, mar = c(0, 0, 0, 0))
circos.par(gap.degree = 1, 
           track.margin = c(0.05, 0.05), 
           points.overflow.warning = FALSE
) 

chordDiagram(mtx, 
             annotationTrack = "grid",
             transparency = 0.5)

# List of labels to add an asterisk to "CETP"  "ZMYM6" "SBSN" 
labels_to_asterisk <- c("CETP", "ZMYM6", "SBSN")
# Labels to color red
labels_red <- c("SBSN")

# Customize the labels to be perpendicular, add asterisks, and color specific labels red
circos.track(track.index = 1, panel.fun = function(x, y) {
  label <- CELL_META$sector.index
  # Append asterisk to all specified labels
  modified_label <-ifelse(label %in% labels_to_asterisk , paste0(label, " ***"), label )
  # Check if the label should also be colored red
  label_color <- ifelse(label %in% labels_red, "red", "black")
  
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2]*3.5, # Adjust position as needed
              modified_label, col = label_color, facing = "outside", 
              niceFacing = TRUE, adj = c(0.5, 0))
}, bg.border = NA)

circos.clear()
dev.off()

#------------------------------------------------------------------



