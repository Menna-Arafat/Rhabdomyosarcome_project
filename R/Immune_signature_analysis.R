# Author: Menna Arafat
# Description: This script applies enrichment analysis to gene signatures within modules identified by WGCNA (Weighted Gene Co-expression Network Analysis).
# It utilizes a hypergeometric test to assess the enrichment of specified gene signatures in the gene modules.
# Input:
#   sig_list : A list of gene sets representing signature genes for specific processes or pathways.
#   modules : A list of WGCNA modules, which are clusters of coexpressed genes that function together.

# Source of the immune signature:
#   "Immune Cell Gene Signatures for Profiling the Microenvironment of Solid Tumors" by Ajit J Nirmal et al.
#   "The Immune Landscape of Cancer" by Vésteinn Thorsson et al.


#=========================================================================================
#set WD
getwd()
setwd("C:/Users/USER/Documents/rhabdo")

#load packages
library(viridis)
library(RColorBrewer)
library(readxl)
library(ape)
library(ggdendro)
library(stats)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(tidyverse)
library(gridExtra)
library(gplots)
library(ggplot2)

#------------------------------------------------------------------------------
#load data
list.files("C:/Users/USER/Documents/rhabdo/output")

#get modules
mapping= read.delim("string_mapping.tsv"   )[,c(2,4)]
module_proteins= read.csv("output/module_proteins.csv") %>% as.list()
#get gene symboles for uniprot ids
modules= list()
for(i in seq_along(module_proteins)){
modules[[i]]= mapping$preferredName[match(module_proteins[[i]], mapping$queryItem)] %>% .[. != "" & !is.na(.)]
}

df= qpcR:::cbind.na(modules[[1]],modules[[2]], modules[[3]], modules[[4]])
df[is.na(df)] <- ""
df= df %>%  as.data.frame()
names(df)= c(names(module_proteins))



module= df$brown
# write.csv(df, "output/modules_genes.csv", row.names = F)

#risk genes associated with rhabdomyosarcoma
risk_genes= read_xlsx("input/risk_genes_Rhabdo_Jack F. Shern_2014.xlsx") %>% select(1) %>% pull() 
#----------------------------------------------------
# gene signatures of immune cells and immune activity and other important tumor activities
immune_sig= read_xlsx("immune_signature_Ajit J. Nirmal_2018.xlsx")[-1,] %>%
                                               set_names(.[1,]) %>% slice(-1)
unique(immune_sig$Signature)
sig_list = split(immune_sig, immune_sig$Signature)
sig_list[[1]]$`HGNC symbols`
#----------------------------------------------------
#module enrichment
total_genes <- 20000 #protein coding genes
results <- list()

for(i in seq_along(sig_list)){
  sig= sig_list[[i]]
  # Intersection between module genes and signature genes
  k <- length(intersect(module,  sig$`HGNC symbols`))
  # Total number of genes in the module
  m <- length(module)
  # Total Number of risk/signature genes in signature i
  n <- length(sig_list[[i]]$`HGNC symbols`)
  # Hypergeometric test
  p_value <- phyper(k - 1, m, total_genes - m, n, lower.tail = FALSE)
  
  # Combine p-value and signature, signature genes, intersected genes into a list
  results[[i]] <- list(
    signature_for= sig$Signature, 
    p_value = p_value,
    signature_genes = sig$`HGNC symbols`, 
    intersected_genes= intersect(module, sig$`HGNC symbols`)
  )
  
  # Assign name to the ith element of results
  names(results)[i] <- sig$Signature
  
}

#subset by Filter function, which is a higher-order function in R that allows you to filter elements of a list based on a predicate function
res_signif <- Filter(function(x) x$p_value < 0.06, results)
res_signif[[1]]
module_genes[[1]] %in% sig_list[[1]]
#----------------------------------------------------------------------------------
#data
#curated signatures from literature from immune subtype clustering paper
list.files("C:/Users/USER/Documents/Immune-Subtype-Clustering/shiny-app/Immune-Subtype-Clustering/data")
load("C:/Users/USER/Documents/Immune-Subtype-Clustering/shiny-app/Immune-Subtype-Clustering/data/comparative_immuneSigs_geneLists4.rda")

sigs1_2_eg2[[30]]$probes
modules= read.csv("output/modules_genes.csv") %>% as.list()
module= modules[["grey"]] %>% .[. != ""]
# Total number of protein coding genes in the genome
total_genes <- 20000 



results <- list()

for(i in seq_along(sigs1_2_eg2)){
  sig= sigs1_2_eg2[[i]]
  # Intersection between module genes and signature genes
  k <- length(intersect(module, sig$probes))
  # Total number of genes in the module
  m <- length(module)
  # Total Number of risk/signature genes in signature i
  n <- length(sig$probes)
  # Hypergeometric test
  p_value <- phyper(k - 1, m, total_genes - m, n, lower.tail = FALSE)
  
  # Combine p-value and probes into a list
  results[[i]] <- list(
    signature_for= sig$src, 
    p_value = p_value,
    signature_genes = sig$probes, 
    intersected_genes= intersect(module, sig$probes)
  )
  
  # Assign name to the ith element of results
  names(results)[i] <-sig$src

}

#subset by Filter function, which is a higher-order function in R that allows you to filter elements of a list based on a predicate function
res_signif <- Filter(function(x) x$p_value < 0.06, results)
res_signif[[1]]


# Convert filtered results to a dataframe
results_df <- do.call(rbind, lapply(res_signif, function(x) {
  data.frame(
    Signature_For = x$signature_for,
    PValue = x$p_value,
    Signature_Genes = paste(x$signature_genes, collapse = ", "),
    Intersected_Genes = paste(x$intersected_genes, collapse = ", "),
    stringsAsFactors = FALSE
  )
}))

write.csv(results_df, "output/immune_signature.csv", row.names = F)



