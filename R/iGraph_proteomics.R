setwd("C:/Users/USER/Documents/rhabdo_proteomics/KD")
list.files(getwd())


#BiocManager::install("estrogen")
# library(affy)
# library(estrogen)
# library(vsn)
# library(genefilter)
library(igraph)
library(tidyr)
library(readxl)
library(openxlsx)
library(org.Hs.eg.db)
library(scales)
library(ggplot2)
#load data
KD= read.xlsx( "input_ppi_for_embryonal.xlsx", colNames =T)
KD= c(KD$hubs_cor_pheno, KD$hubs_blue)
#mapping of KD
mapping= read.delim("idmapping.tsv")
names(mapping)= c("From", "To")
KD= mapping$To[match(KD, mapping$From)]
KD

net= read.delim("ppi_embryonal.tsv")[,1:2]

#build igraph object from edge list
graph =graph_from_edgelist(net %>%as.matrix(), directed = FALSE)

#centrality measures
# Degree Centrality, edges cconnected to a node
degree_centrality <- degree(graph, mode = "all")  %>% sort(.,decreasing = T) 
print(names(degree_centrality )[1:5])
DC= names(degree_centrality)[1:5] %>% as.vector()


# Betweenness Centrality, The number of shortest paths that pass through a vertex
betweenness_centrality <- betweenness(graph, normalized = TRUE)  %>% sort(.,decreasing = T) 
print(names(betweenness_centrality)[1:5])
BC= names(betweenness_centrality)[1:5] %>% as.vector()

# Eigenvector Centrality, The importance of a node based on the importance of its neighbors.
eigenvector_centrality <- eigen_centrality(graph)$vector %>% sort(.,decreasing = T) 
print(eigenvector_centrality)
eigC= names(eigenvector_centrality)[1:5] %>% as.vector()

#closeness centrality
closeness_centrality= closeness(graph) %>% sort(.,decreasing = T)
closeness_centrality

# PageRank,Nodes are considered important if they are linked to by other important nodes, with a damping factor accounting for random jumps.
pagerank <- page_rank(graph)$vector %>% sort(.,decreasing = T) 
print(pagerank)
prC= names(pagerank)[1:5] %>% as.vector()

all_drivers= c(KD,DC, BC, eigC, prC ) %>% unique() %>% na.omit(.) %>% as.character()
all_drivers
#visualization
# Create a color vector, defaulting to black
node_colors <- rep("gray80", vcount(graph))

# Set the color of specific nodes to red
node_colors[V(graph)$name %in% all_drivers] <-  "orange"
node_colors
#change label position
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
n= vcount(graph)
lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)


png("igraph_ppi_emb.png", width=17000, height=18000, res= 1200)
plot(graph, edge.arrow.size=.5, vertex.color= node_colors, vertex.size=5, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.8, vertex.label.dist= 1, vertex.label.degree=lab.locs,
     layout=layout_in_circle, #layout_with_fr
     main="PPI for Embryonal RMS with KEy Drivers") 

dev.off()

#WGCNA networks
net_turquoise= read.delim("output/networks/CytoscapeInput-edges-turquoise.txt")[,1:2]
net_brown= read.delim("output/networks/CytoscapeInput-edges-brown.txt")[,1:2]
hubs= read.csv("output/hubproteins_gene_names.csv")

mapping= read.delim("string_mapping.tsv" )[,c(2,4)]
#mapping
networks= list(net_brown, net_turquoise)
network_names = c("brown",  "turquoise")
mapped_networks= list()

for (i in seq_along(networks)){
  net= networks[[i]]
  mapped1= mapping$preferredName[match(net[[1]]  , mapping$queryItem)]
  mapped2= mapping$preferredName[match(net[[2]] , mapping$queryItem)]
  mapped= cbind(mapped1, mapped2) %>% as.data.frame(.)
  mapped = na.omit(mapped)
  mapped_networks[[network_names[i]]]= mapped
  
}
#save
for (name in network_names) {
  
  file_name = paste0("output/",name, "_wgcna_network_mapped" ,".txt")
  
  write.table(mapped_networks[[name]], file_name, row.names = FALSE, sep= "\t")
}

brown= mapped_networks[["brown"]]
turquoise= mapped_networks[["turquoise"]]
hub_brown= hubs$brown
hub_turquoise= hubs$turquoise %>% .[. != ""]

graph= graph_from_edgelist(turquoise %>%as.matrix(), directed = FALSE)

nodes_to_label <- hub_turquoise

# Create a vector for labels, defaulting to empty strings
labels <- setNames(rep("", vcount(graph)), V(graph)$name)

# Set labels for selected nodes
labels[nodes_to_label] <- nodes_to_label
labels
V(graph)$color <- ifelse(V(graph)$name %in% nodes_to_label,"darkgrey", "grey")

# Define node sizes (replace with your own data, perhaps normalized gene expression levels)
V(graph)$size <- ifelse(V(graph)$name %in% nodes_to_label, 2, 1)

# Define edge widths based on co-expression values (replace with your actual co-expression data)
E(graph)$width <- runif(ecount(graph), 0.5, 2)

# Define edge colors if needed
E(graph)$color <- "grey"

# Customize plot parameters
png("plots/igraph_turquoise_network.png", width=19000, height=19000, res= 1200)
plot(graph, layout=layout_with_fr(graph), # Use Fruchterman-Reingold layout
     vertex.label=labels, 
     vertex.label.color="black",
     vertex.label.cex=0.8, 
     vertex.label.dist=0.5, # Set label distance from nodes
     vertex.frame.color="white", # Remove the node borders
     vertex.color=V(graph)$color,
     vertex.size=V(graph)$size, # Use the sizes defined earlier
     edge.width=E(graph)$width, # Use the edge widths defined earlier
     edge.color="grey", # Use the edge colors defined earlier
     main="Interaction Network of Turquoise Module") # Add a title to the plot

dev.off()

#-------------------------------------------------------------------------------
#adjacency
graph[]
get.adjacency(net, attr="weight", sparse=F)
#all vertices
V(graph)$name
#all edges
E(graph)
#get neighbours
n= neighborhood(graph, order=1, V(graph)$name %in% c( "APCS", "SERPINA3","A1BG" ))

names(unlist(n))    

#get edges between certain node and another node, and their count/length         
E(graph)[ V(graph)[name=="SERPINA3" ] %--% V(graph)[name=="A1BG"] ] %>% length()   
