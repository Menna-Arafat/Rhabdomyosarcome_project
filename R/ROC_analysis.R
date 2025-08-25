#AUther: MENNa Arafat
getwd()
setwd("C:/Users/USER/Documents/rhabdo")
#load packages
library(viridis)
library(RColorBrewer)
library(pROC)
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
allowWGCNAThreads()          # allow multi-threading (optional)

#load data
list.files("C:/Users/USER/Documents/rhabdo/output")
data= read.csv("input/ALL_PROT-ACC-zero.csv" )[-1,]
data[]= lapply(data, as.numeric)
mapping= read.delim("string_mapping.tsv"   )[,c(2,4)]
hubs_alv= read.csv("output/hubproteins_gene_names_Alv.csv")
hubs_emb= read.csv("output/hubproteins_gene_names_Emb.csv")
hubs_alv= hubs_alv$brown %>% as.vector() %>% .[. != "" & !is.na(.)]
hubs_emb= hubs_emb$brown %>% as.vector() %>% .[. != "" & !is.na(.)]
hubs= c(hubs_alv, hubs_emb) %>% unique()

data= read.csv("input/rhabdo_gene_symbols.csv") %>% column_to_rownames("Gene")
hubs= read.csv("output/hubproteins_gene_names.csv")[,1]

# find area under the curve ,Example for one gene for example "P11597"
data_d= data %>% dplyr::select(1:35) %>% 
  filter(row.names(.)== "P11597" ) %>% t(.)%>%
  as.data.frame() %>% mutate(class = ifelse(grepl("Alv", row.names(.)), 1, 2))

#Input for roc function is expression profile for the gene and the corresponding label of the group and it should be numeric
roc_result <- roc(response = data_d$class, predictor = data_d[,1])
auc(roc_result)

#loop over vector of hub genes, for every gene, extract expression profile (row),
#then transpose it to create an opposing column thatindicate the patient's group/ class
res= list()
auc_values <- numeric(length(hubs))
for(i in 1:length(hubs)){
  protein= hubs[i]
  data_d= data %>% select(matches("Alv|Emb")) %>% filter(row.names(.)==protein ) %>% t(.)%>%
          as.data.frame() %>% mutate(class = ifelse(grepl("Alv", row.names(.)), 1, 2))
  roc_result = roc(response = data_d$class, predictor = data_d[,1])
  res[[i]] = roc_result
  auc_values[i]= auc(roc_result)
}
names(res) = hubs
auc_values
diag= cbind(hubs, auc_values) %>% as.data.frame(.) %>% arrange(desc(auc_values))
#diag$marker_for= c(rep(paste0("Alveolar"),6), rep(paste0("Embryonal"),11))
#diag$gene_name= mapping$preferredName[match(diag$hubs, mapping$queryItem)]
write.csv(diag, "output/ROC_results.csv",row.names = F)


#plots
#plot AUC for single gene
#roc_result is the outputof roc() function
plot(roc_result)
# Plot the ROC curves combined
names(res)
brewer.pal(6, "Set3")
col= viridis(6, alpha = 1, begin = 0, end = .9, direction = 1, option = "inferno")

#plot
plot.roc(res[["CETP"]], col="#000004FF", main="ROC Curves for Top Hub Genes",
         legacy.axes= T,xlim= c(.8,0))
lines(res[["SBSN"]], col="#390963FF" )
lines(smooth(res[["ZMYM6"]]), col="#84206BFF")
lines(smooth(res[["F11"]]), col="#C9404AFF")
lines(smooth(res[["FGA"]]), col="#F67F13FF")
lines(smooth(res[["CFH"]]), col= "#F6D645FF")
legend("bottomright", legend=c("CETP:  1", "SBSN:  1", "ZMYM6: 0.85", "F11:  0.71", "FGA:  0.69", "CFH:  0.68"), 
       col= col, lwd=2, pch=19, title= "AUC Scores", title.font= 2)

#----------------------------------------------------
dev.off()
