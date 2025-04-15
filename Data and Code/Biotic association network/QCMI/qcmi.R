library("phyloseq")
library("ggClusterNet")
library("tidyverse")
library("WGCNA")
library("igraph")
library("Matrix")
library("vegan")
library("magrittr")
library("reticulate")
library("SpiecEasi")
library("ggstatsplot")
library("qcmi")
library("dplyr")



################The ecological association of diffusion limitation is classified based on geographical distance.
####result_dl= assigned_process(link_table_row=edgelist, OTUabd=otu_table, p=0.05 , data= geo,cutoff=0, method="dl")


################Filtered environmental ecological association classification, based on selected environmental differences


env<-read.table("env.txt",header = T)
otu<-read.table("ASV.txt",header = T)
edge<-read.table("edge.txt",header = T)
geo<-read.table("geo.txt",header = T)

result_dl= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= geo,cutoff=0, method="dl")

dl_link<-row.names(result_dl)

dledge=edge[dl_link,]

######Detecting Environmental Redundancy
library(Hmisc)
plot(varclus(as.matrix(env) ))


result_AP= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['AP'],cutoff=0, method="ef")
result_TP= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['TP'],cutoff=0, method="ef")
result_NO= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['NO3'],cutoff=0, method="ef")
result_NH= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['NH4'],cutoff=0, method="ef")
result_Salt= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['Salt'],cutoff=0, method="ef")
result_pH= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['pH'],cutoff=0, method="ef")
result_SWC= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['SWC'],cutoff=0, method="ef")
result_TN= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['TN'],cutoff=0, method="ef")
result_SOC= assigned_process(link_table_row=edge, OTUabd=otu, p=0.05 , data= env['SOC'],cutoff=0, method="ef")


length(row.names(result_SOC))
length(row.names(result_TN))

total_link=row.names(edge)
ef_link=union(as.character(row.names(result_AP)),as.character(row.names(result_TP)))
ef_link=union(ef_link,as.character(row.names(result_NO)))
ef_link=union(ef_link,as.character(row.names(result_NH)))
ef_link=union(ef_link,as.character(row.names(result_Salt)))
ef_link=union(ef_link,as.character(row.names(result_SWC)))
ef_link=union(ef_link,as.character(row.names(result_pH)))
ef_link=union(ef_link,as.character(row.names(result_TN)))
ef_link=union(ef_link,as.character(row.names(result_SOC)))


bi_link=setdiff(total_link,ef_link)

bi_link=setdiff(bi_link,dl_link)


efedge=edge[ef_link,]
write.csv(efedge,"S5_ef_edge.csv")
length(row.names(efedge))

biedge=edge[bi_link,]
write.csv(biedge,"S5_bi_edge.csv")
length(row.names(biedge))

ig.bi=graph_from_edgelist(as.matrix(biedge[,1:2]) , directed = FALSE)
ig.bi=set_edge_attr(ig.bi, 'weight', index = E(ig.bi),  as.numeric(biedge[,3]) )

##result_bi= qcmi(igraph= ig.bi, OTU= otu_table, pers.cutoff=0)

result_bi= qcmi(igraph= ig.bi, OTU= otu, pers.cutoff=0)


int<-data.frame(result_bi$interaction.neg,result_bi$interaction.pos)

write.csv(int,"int_S5.csv")

