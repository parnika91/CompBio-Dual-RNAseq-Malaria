# Script to find modules in the network; 
# locate immune cell marker genes and
# calculate the enrichment of these markers in the network


library(igraph)
library(grid)
library(limma)
library(EnsDb.Hsapiens.v79)

options(echo=TRUE)
args <- commandArgs(TRUE)
dataset <- args[1] #, eg blood_all_bipartite for overall and blood_core_network for core network

############################# communities
network <- dataset
data <- readRDS(paste0(network, ".rds"), collapse = '')

ig <- graph_from_data_frame(data, directed = F)
vnames = c(unique(as.character(data[,1])), unique(as.character(data[,2])))
hnum = length(unique(as.character(data[,1])))
pnum = length(unique(as.character(data[,2])))
E(ig)$weight = data[,3]

# clustering using edge betweenness algorithm
mst <- mst(ig, algorithm="prim")
mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)

commSummary <- data.frame(
  mst.communities$names,
  mst.communities$membership)
colnames(commSummary) <- c("Gene", "Community")

# check how the number of modules affect modularity of the network
m <- vector()
for (s in 0:200)
{
  memb <- cut_at(mst.communities, no=s)
  m <- c(m, modularity (ig, memb, weights=NULL))
} 

png("modularity_betwenness.png", width = 30, height = 10, unit = "cm", res = 300)
plot(0:200, m, col="blue",xlab="Modules",ylab="Modularity")
dev.off()

#### communities with fast greedy algorithm
comm.fg <- cluster_fast_greedy(ig, merges = TRUE, modularity = TRUE,
  membership = TRUE, weights = abs(E(ig)$weight))

m <- vector()

for (s in 0:200)
{
  memb <- cut_at(comm.fg, no=s)
  m <- c(m, modularity(ig, memb, weights=NULL))
} 

png("modularity_fast_greedy.png", width = 30, height = 10, unit = "cm", res = 300)
plot(0:200, m, col="blue",xlab="Modules",ylab="Modularity")
dev.off()


###### communities with leading eigenvector ############

comm.eig <- cluster_leading_eigen(ig, steps = -1, weights = E(ig)$weight,
  start = NULL, options = arpack_defaults, callback = NULL,
  extra = NULL, env = parent.frame())

m <- vector()
for (s in 0:200)
{
  memb <- cut_at(comm.eig, no=s)
  m <- c(m, modularity (ig, memb, weights=NULL))
} 
png("modularity_leading_eigen.png", width = 30, height = 10, unit = "cm", res = 300)
plot(20:200, m, col="blue",xlab="Modules",ylab="Modularity")
dev.off()


##################### convert gene IDs of the marker gene names to ensembl IDs for humans.

# read marker genes file
# Immune_cell_marker_for_R.csv is the supplementary file 4 from Vallejo et al., 2018
markers <- list()
markers[["B_cells"]] <- read.columns(file = "Immune_cell_marker_for_R.csv", sep = ',', required.col = 1)
markers[["monocytes"]] <- read.columns(file = "Immune_cell_marker_for_R.csv", sep = ',', required.col = 2)
markers[["myeloid_DC"]] <- read.columns(file = "Immune_cell_marker_for_R.csv", sep = ',', required.col = 3)
markers[["neutrophils"]] <- read.columns(file = "Immune_cell_marker_for_R.csv", sep = ',', required.col = 4)
markers[["NK_cells"]] <- read.columns(file = "Immune_cell_marker_for_R.csv", sep = ',', required.col = 5)
markers[["T_cells"]] <- read.columns(file = "Immune_cell_marker_for_R.csv", sep = ',', required.col = 6)

# remove empty rows
for(i in names(markers))
	markers[[i]] <- markers[[i]][!apply(markers[[i]] == "", 1, all),]

ensIDs <- list()
for(j in names(markers))
{
	# convert IDs
	ensIDs[[j]] <- ensembldb::select(EnsDb.Hsapiens.v79, keys= markers[[j]], 
		keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
	# remove "LRG" IDs
	ensIDs[[j]] <- ensIDs[[j]][-grep(pattern = "LRG_", ensIDs[[j]]$GENEID),]
}

# convert commSummary to human ens IDs
pOG <- read.delim("parasite_orthogroups.txt", stringsAsFactors=FALSE)
hOG <- read.delim("host_orthogroups.txt", stringsAsFactors=FALSE)

ens <- c()
for(i in 1:nrow(commSummary))
{
	if(grepl(pattern = "h_OG", commSummary$Gene[i]))
		ens[i] <- hOG[grep(pattern = commSummary$Gene[i], hOG$Orthogroup), "human"]
	if(grepl(pattern = "p_OG", commSummary$Gene[i]))
		ens[i] <- pOG[grep(pattern = commSummary$Gene[i], pOG$Orthogroup), "Pberghei"]
}

ens_commSummary <- data.frame(ens = ens, community = commSummary$Community)

# make a vector that registers if any of these community genes are one of the marker genes.
# If yes, put the name of the cell that this is a marker

clusters <- ens_commSummary
cell_marker <- data.frame(Marker = rep("", nrow(clusters)))

for(i in 1:nrow(clusters))
{
	for(j in names(ensIDs))
	{
		if(any(ensIDs[[j]]$GENEID %in% clusters[i,1]))
			cell_marker[i,1] <- j

	}
}

clusters <- data.frame(Gene = ens_commSummary$ens, Membership = ens_commSummary$community, Marker = cell_marker)
# save to use in network in Cytoscape
write.csv(clusters, paste0(network, "_betweenness_communities.csv", collapse = ''), row.names = F, quote = F)

# calculate the enrichment of each a cell type
all_host_genes <- hOG[, "human"]
network_genes <- hOG[ hOG$Orthogroup %in% data[,1], "human"]

for(i in names(ensIDs))
{
  print(names(ensIDs)[j])
  print(table(In_network = all_host_genes %in% network_genes, Cell_type = all_host_genes %in% ensIDs[[i]]$GENEID))
  fisher.test(table(In_network = all_host_genes %in% network_genes, Cell_type = all_host_genes %in% ensIDs[[i]]$GENEID))
}