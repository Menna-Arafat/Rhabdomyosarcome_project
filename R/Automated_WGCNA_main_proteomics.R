#Code Description:

# This script automates the process of conducting Weighted Gene Co-expression Network Analysis (WGCNA) to
# identify modules of highly correlated genes and relate them to external traits

#project name: proteomic and metabolomic profiling for RhabdoMyosarcoma subtypes (Alveolar and Embryonal)

#Author: Menna Arafat

#Date: 3/8/2024

#---------------------------------------------------------------------------------------------------------
#pre-requiste libraries
# install.packages("WGCNA")
# install.packages("flashClust")
# install.packages("ggplot2")
# install.packages("data.table")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")

setwd("C:/Users/USER/Documents/rhabdo_proteomics")

#create directtories
dir.create("input")
dir.create("output")
dir.create("plots")

#load libraries
library(WGCNA)
library(stats)
library(flashClust)
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
library(readxl)
library(openxlsx)
library(caret)
allowWGCNAThreads()          


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
#metadata
metadata= data.frame(Sample= colnames(data) ,
                     condition= gsub("[0-9_]+", "", colnames(data))) %>%
                     column_to_rownames(var = "Sample")

metadata$condition= factor(metadata$condition, levels = c("Ctrl", "Alv", "Emb"))
design= model.matrix(~ 0+condition , metadata)
design
#-------------------------------------------------------------------------------
##Applying WGCNA pipeline
#choose the power that the adjacency matrix raised for that allow for scale free network
# Choose a set of soft thresholding powers
powers = c(1:15)  # in practice this should include powers up to 20.
# choose power based on SFT criterion
sft = pickSoftThreshold(t(data), powerVector = powers, networkType = "signed")
# Plot the results:
png("plots/elbow_plot_signed.png", width = 1700, height = 1000, res=300)
par(mfrow = c(1, 2))
# SFT index as a function of different powers
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "red", cex = .7)
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.8, col = "red")
# Mean connectivity as a function of different powers
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red", cex = .7)

dev.off()

#-------------------------------------------------------------------------------
#set your parameters
power= 10  # we raise the adjacency matrix to this power to get rid of weak edges, hence emphasizing only strong correlations and ensure that the network degree follow a scale free distribution, recommended to be less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks  
minModuleSize = 35 #minimum module size is the minimum size/ number of genes per module 
networkType = "signed" # as for unsigned network, negative correlations are treated the same as positive correlations. While for signed network, it gives more weight to positive correlations and consider negative correlations irrelevant.
metadata_binary= design # metadata is constructed in a way that samples are encoded as 0, 1, for example, 0 in case of control and 1 in case of disease condition.
hub_threshold= 0.8 # reflect intra-modular connectivity of the node
GS_threshold= 0.2 #(.2) # Gene significance threshold reflect the correlation between the gene expression profile of a hub gene and the phenotype of interest
mapping_file= read.table("input/idmapping.tsv", header = T) #needs to be 2 columns that are the query id and the matched id of interest


#------------------------------------------------------------------------------  
Run_WGCNA= function(data, metadata_binary, power,minModuleSize, hub_threshold ,GS_threshold, networkType,
                     mapping_file=NULL){
  

  data= as.data.frame(data)
  data[]= lapply(data, as.numeric)
  datExpr= t(data) #so that samples become in row
  #This function checks data for missing entries, entries with weights below a threshold, and zero-variance genes, 
  goods <- goodSamplesGenes(datExpr, verbose = 3)
  datExpr= datExpr[goods$goodSamples== TRUE, goods$goodGenes == TRUE ]
  
  net = blockwiseModules(datExpr, corType = "pearson", maxBlockSize = 5000, 
                         networkType = networkType , power = power, minModuleSize =minModuleSize,
                         mergeCutHeight = 0.25, 
                         numericLabels = F, saveTOMs = TRUE, 
                         pamRespectsDendro = FALSE, saveTOMFileBase = "TOM")
  #table(net$colors)
  print("module detection is done")
  
  # A data frame with module eigengenes 
  module_eigengenes <- net$MEs
  write.csv(module_eigengenes, "output/module_eigengenes.csv")
  names(module_eigengenes)= gsub("ME", "", names(module_eigengenes) )
  print("module_eigengenes file is saved")
  
  #gene module assignment
  module.gene.assign= net$colors #vector of colors that assign each gene to its corresponding module
  write.csv(module.gene.assign, "output/gene_module_assignment.csv")
  print("gene_module_assignment file is saved")
  
  #Extract proteins for each modules
  module_names= names(module_eigengenes)
  module_lists= list()
  
  for (i in module_names){
    module_lists[[i]]=  module.gene.assign[module.gene.assign== i] %>% names()
  }
  
  df <- ldply(module_lists, rbind) %>% t()
  df[is.na(df)] <- ""
  df= as.data.frame(df)
  names(df) = df[1,] %>% as.vector()
  modules= df[-1,]
  write.csv(modules, paste0("output/", "modules_gene_id.csv"), row.names = FALSE)
  print("modules file is saved!")
  
  
  #calculate gene significance (a measure that show to what degree the gene is correlated with the phenotype of interest based on pearson correlation, and the authors consider GS > .2 is good enough) 
  nSamples <- nrow(datExpr)
  nGenes <- ncol(datExpr)
  gene.signf.corr <- cor(datExpr, metadata_binary, use = 'p')
  gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
  write.csv(gene.signf.corr, "output/gene.trait.corr.csv")
  gene.signf.corr= gene.signf.corr  %>% as.data.frame()
  #Module membership (intra-modular connectivity)
  #For a node i. its module membership kME(q) specifies how close node i is to module q. 
  #KME (module membership) defined as correlation between gene expression profile and the module eigengene (the first principal component or the mean expression profile of a module)
  #For example, if the module membership of a gene i (MMblue.i) is close to 1 or -1, this indicate that the gene has a positive or a negative relationship with the module eigengene, i.e., this gene is highly connected to the blue module genes.

  datKME = signedKME(datExpr, module_eigengenes)
  write.csv(datKME, "output/module_membership.csv")
  
  
  #extract hub genes
  #Extract hubs for each modules
  module_names= names(module_eigengenes)
  hub_lists <- vector("list", length(module_names))
  names(hub_lists)= names(module_eigengenes)
  threshold= hub_threshold
  
  for(j in 1:length(module_names)){
    hub_mod= row.names(datKME)[datKME[,j] >= threshold]
    hub_lists[[j]]= intersect(unlist(modules[j]),  hub_mod)
  }
  
  df= ldply(hub_lists, rbind) %>% t() 
  df[is.na(df)] <- ""
  df= as.data.frame(df)
  names(df) = paste0(unlist(df[1,]), "_hub") 
  hubs= df[-1,]
  write.csv(hubs, paste0("output/", "module_hubs.csv"), row.names = FALSE)
  print("hubs file is saved!")
  
  #Extract hubs correlated with phenotype
  hub_pheno <- vector("list", length= ncol(gene.signf.corr))
  names(hub_pheno)= names(gene.signf.corr)
  GS_threshold= GS_threshold
  
  for(j in 1:ncol(gene.signf.corr)){
    
    all_hubs= hub_lists %>% unlist() %>% .[. !="" & !is.na(.)]
    hub_pheno[[j]]= row.names(gene.signf.corr)[gene.signf.corr[,j] > GS_threshold] %>%  
                    .[. %in% all_hubs]
  }
  df= ldply(hub_pheno, rbind) %>% t() 
  df[is.na(df)] <- ""
  df= as.data.frame(df)
  names(df) = paste0(unlist(df[1,]), "_hub") 
  hubs= df[-1,]
  write.csv(hubs, paste0("output/", "phenotype_cor_hubs.csv"), row.names = FALSE)
  print("hubs_pheno file is saved!")
  
  #MAPPING 
  if(!is.null(mapping_file)){
  #mapping of hubs for modules and for phenotype
  mapping= mapping_file
  names(mapping)= c("query", "name")
  
  modules= module_lists
  modules_mapped= vector("list", length= length(module_lists))
  names(modules_mapped) = names(module_lists)
  
  for(i in seq_along(modules)){
    modules_mapped[[i]]= mapping$name[match(modules[[i]], mapping$query) ]
  }
  df= ldply( modules_mapped, rbind) %>% t() %>% as.data.frame()
  df[is.na(df)]= ""
  names(df)= as.character(df[1,])
  write.csv(df[-1,], "output/modules_gene_name.csv", row.names = F)
  
  
  combined_hubs= list(mod= hub_lists, pheno= hub_pheno)
  
  mapped_lists = list(mod= setNames(vector("list", length = length(hub_lists)), names(hub_lists)),
                     pheno= setNames(vector("list", length = length(hub_pheno)), names(hub_pheno)))
  
  unmapped_lists = list(mod= setNames(vector("list", length = length(hub_lists)), names(hub_lists)),
                       pheno= setNames(vector("list", length = length(hub_pheno)), names(hub_pheno)))
  
  for(i in seq_along(combined_hubs)){
  for (j in seq_along(combined_hubs[[i]])) {
      mapped_lists[[i]][[j]] <- mapping$name[match(combined_hubs[[i]][[j]], mapping$query)] 
      unmapped_lists[[i]][[j]] = combined_hubs[[i]][[j]][is.na(match(combined_hubs[[i]][[j]], mapping$query)) ]
  }
    df_mapped= ldply(mapped_lists[[i]], rbind) %>% t() %>%  as.data.frame()
    df_mapped[is.na( df_mapped)]= ""
    names(df_mapped)=  as.vector(df_mapped[1,])
    mapped_hubs=  df_mapped[-1,]
    write.csv(mapped_hubs, paste0("output/", names(combined_hubs)[i],"_hubs_names.csv"), row.names = FALSE)
    print("hubs_name file is saved!")
    
    df_unmapped= ldply(unmapped_lists[[i]], rbind) %>% t() %>%  as.data.frame()
    df_unmapped[is.na( df_unmapped)] <- ""
    names(df_unmapped)=  as.vector(df_unmapped[1,])
    unmapped_hubs=  df_unmapped[-1,]
    write.csv(unmapped_hubs, paste0("output/", names(combined_hubs)[i],"_hubs_unmapped.csv"), row.names = FALSE)
    print("hubs_unmapped file is saved!")
    
  }
  }
  else 
    { print("Supply a mapping id file from uniprot database! and Set mapping_source= \"uniprot\"")}
  
  #export networks for cytoscape
  module_names= names(module_eigengenes)
  for( i in   module_names) {
    datexpr_mod = datExpr[,module.gene.assign == i]
    TOM_mod = TOMsimilarityFromExpr(datexpr_mod, power = power, networkType = "signed", TOMType="signed Nowick");
    genes = colnames(datexpr_mod)
    dimnames(TOM_mod) = list(genes, genes)
    # nTop = 30;
    # IMConn = softConnectivity(datExpr[, probes]);
    # top = (rank(-IMConn) <= nTop)
    cyt = exportNetworkToCytoscape(TOM_mod,
                                   edgeFile = paste("output/CytoscapeInput-edges-.02tom", paste(i, collapse="-"), ".txt", sep=""),
                                   #nodeFile = paste("output/CytoscapeInput-nodes-", paste(i, collapse="-"), ".txt", sep=""),
                                   weighted = TRUE,
                                   threshold = 0.02, #threshold for including edges in the output, Default is 0.02
                                   nodeNames = genes,
                                   altNodeNames = genes);
  }
  print("networks exported!")
  print("ALL IS DONE!")
  return(net)
  
}
res= Run_WGCNA(data, metadata_binary=metadata_binary, power=power, minModuleSize= minModuleSize,networkType = networkType,
               hub_threshold= hub_threshold, GS_threshold=GS_threshold, mapping_file=mapping_file )

