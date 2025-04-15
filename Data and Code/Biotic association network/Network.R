rm(list=ls())


####networks
library(WGCNA)
library(igraph)
library(ggplot2)

#Correlation matrix calculation
#It is necessary to remove the taxa with more 0 values in the group before analysis
occor<-corAndPvalue(t(otu),method="pearson",use="p")
occor.r = occor$cor 
occor.p = occor$p 
occor.r[occor.p>0.01|abs(occor.r)<0.8] = 0
diag(occor.r)<-0
occor.r[is.na(occor.r)]<-0
#edge
sum(abs(occor.r)>0)/2 
#node
sum(colSums(abs(occor.r))>0)
#Export the data and plot using Gephi
write.csv(occor.r,file="network.csv")
##### What is obtained here is only the co-occurrence network, and two steps of filtering are needed to obtain the microbial biological association network.







#### keystone nodes 
library(reshape2)
edge<-read.table("clipboard",header = T) #Import the edge file obtained by the biological association network
#The edge file is converted into an association matrix
occor.r<-dcast(edge,edge$Source~edge$Target,mean)
row.names(occor.r)<-occor.r$`edge$Source`
occor.r[is.na(occor.r)]=0
occor.r<-occor.r[-1]
occor.r[abs(occor.r)>0]=1
adjacency_unweight<-occor.r

# Obtain an undirected network with no weights
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)

#computational node
V(igraph)$degree <- degree(igraph)

set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#integrated data
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list) 

#Calculate within-module connectivity (Zi) and among-module connectivity (Pi)
source('hub_node.r')

row.names(nodes_list)<-nodes_list$nodes_id
nodes_list<-nodes_list[-1]

zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

zi_pi <- na.omit(zi_pi)   # remove NA
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type,shape=type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)
