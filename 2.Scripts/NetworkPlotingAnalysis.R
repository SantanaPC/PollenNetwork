library("ggplot2")
library("tidyverse")
library("igraph")

#loading the matrix of interaction
pollen_net <- read.csv("PollenNetwork/1.Data/Network_data.csv", header = T, sep="," , row.names = 1)
pollen_one_net <- as.one.mode(pollen_net)
graph<-graph_from_biadjacency_matrix(pollen_net, weighted = T)
E(graph)$weight#looking at the interactions

net_matrix<-as_adjacency_matrix(graph, sparse = FALSE)
net_matrix[which(is.na(net_matrix) == T)]<-0

#loading the data about the species
net_info <- read_csv("PollenNetwork/1.Data/Species_data.csv")

# General metrics
plants<-net_info|> filter(resource==1)|> nrow()
animals <- net_info|> filter(resource==0)|> nrow()
size <- plants*animals
richness <- plants+animals
links <- net_matrix|> nrow()
conn <- round(links/size, 3)#network not very connected

isResource <- net_info$resource |> as.logical() # 0/1 converted to FALSE/TRUE
isHive <- net_info$color |> as.factor()

V(graph)$type <- !(isResource) 
E(graph)$arrow.mode = "-"
V(graph)$hive <- isHive

colrs <- c("#714925","#eec04b", "#CF6B13","#F69F22")
status<- c("#714925","#eec04b", "#d7191c", "#CF6B13","#F69F22")

#V(graph)$name <- net_info$species
V(graph)$group <- net_info$hive
V(graph)$colHive <- colrs[V(graph)$hive]
V(graph)$status <- as.factor(net_info$cons_status)
V(graph)$colSta <- status[V(graph)$status]

modulos<-cluster_spinglass(graph)

layout <- layout_as_bipartite(graph)

#png("bipartite_graph.png", width = 800, height = 600)

#length(degree(graph))

plot(graph, 
     edge.width=0.3*igraph::degree(graph), 
     vertex.color=vertex_attr(graph)$colSta, 
     vertex.size=0.6*igraph::degree(graph), 
     vertex.label=NA,
     layout=layout_as_bipartite,
     #vertex.label=vertex_attr(graph)$group,
     vertex.frame.color="#777777", 
     vertex.label.color="black", 
     vertex.label.cex=0.8, 
     #vertex.label.dist=2, 
     edge.curved=0.1,
     #mark.groups= modulos,
     asp=0.4
)

legend(x=-1, y=-0.2, c("Nannotrigona testaceicornis", "Native Plants", "Non-native Plants", "Paratrigona opaca", "Tetragonisca angustula"), pch=21, col="#777777", pt.bg=c("#714925","#eec04b", "#d7191c", "#CF6B13","#F69F22"), pt.cex=2, cex=.8, bty="n", ncol=1)

#heatmap from the interaction matrix to visualise the species which are more connected
pollen_net<-as.matrix(pollen_net)

heatmap(t(pollen_net), Rowv=NA, Colv=NA, col=colorRampPalette(c("#f5f5f5","#80cdc1", "#35978f", "#003c30"))(100), scale="none", margins=c(10,5))

##############
#############

#network with the pooled data (hives combined)
pool_net <- read.csv("PollenNetwork/1.Data/Network_pooled_data.csv", header = T, sep="," , row.names = 1)
pool<-graph_from_biadjacency_matrix(pool_net, weighted = TRUE)
E(pool)$weight
length(V(pool)$name)

#selecting the info for the pooled data
net_pool_info<-net_info[c(1:50,53,56),]#always check if the order is right

bees <- c("#c7eae5","#8c510a", "#bf812d","#d8b365")
status<- c("darkgreen", "#d7191c")

isResource <- net_pool_info$resource %>% as.logical() # 0/1 converted to FALSE/TRUE

# General metrics
plants<-pool_net|> nrow()
animals <- pool_net|> ncol()
size <- plants*animals
richness <- plants+animals
# Calculando métricas de rede, incluindo conectância
networklevel(pool_net, index = "connectance")
networklevel(pool_net, index = "H2")
networklevel(pool_net, index = "NODF")
grouplevel(pool_net, index = "cluster coefficient")
grouplevel(pool_net, index = "niche overlap")
moduls<-computeModules(pool_net)
moduls@likelihood
plotModuleWeb(moduls)

# Exibindo o resultado
result

V(pool)$type <- !(isResource) 
V(pool)$name <- net_pool_info$species
V(pool)$hive <- factor(net_pool_info$color, levels = c("P", "N", "T", "PO"))
V(pool)$colBee <- bees[V(pool)$hive]
V(pool)$status <- as.factor(net_pool_info$cons_status)
V(pool)$colSta <- status[V(pool)$status]

layout <- layout_as_bipartite(pool)

#naming species
vertex_labels<-net_pool_info$group

#plot of the pooled data
plot(pool, 
     edge.width=0.4*igraph::degree(pool),
     vertex.color=vertex_attr(pool)$colSta, 
     vertex.size=8,
     #vertex.label=NA, 
     edge.color="grey50",
     vertex.frame.color="#777777", 
     vertex.label.color="black", 
     vertex.label.cex=1, 
     #main="Pool network",
     vertex.label = vertex_labels,  
     vertex.label.font=1,
     vertex.label.dist=2,
     vertex.label.degree=pi/2,
     asp=0.5,
     #edge.curved=0.2,
     #mark.groups= cluster_spinglass(pool),
     layout=layout_as_bipartite
     #layout = vertical_layout,
) 

# Add legend
legend(x=0.8, y=-0.5, c("Plants", "Nannotrigona testaceicornis", "Paratrigona opaca", "Tetragonisca angustula"), pch=21, col="#777777", pt.bg=colrs <- c("#c7eae5","#8c510a", "#bf812d","#d8b365"), pt.cex=1, cex=.6, bty="n", ncol=1)

legend(x=0.8, y=-0.5, c("Native", "Non-native"), pch=21, col="#777777", pt.bg=c("#2c7bb6", "#d7191c"), pt.cex=2, cex=.8, bty="n", ncol=1)

