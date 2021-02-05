# Script to find out which bipartite interactions are present across multiple datasets
# and to create the core network

library(tidyverse) # glue package of tidyverse is required, installs by invoking tidyverse
library(igraph)

# collection of studies that constitute the blood overall dataset
studies <- c("blood_all",
  # human studies
  "DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int", #hpv_bl is the collection of human - P.vivax blood runs
  # mouse studies
  "ERP110375_str", "ERP004598_int", "m_bl_int", # m_bl is a collection of mouse - Plasmodium studies
  # monkey studies
  "SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int")

for(i in studies)
  readRDS(paste0(substr(i, 1, 9), "/cor/", i, "_bipartite.rds", collapse = "")) %>%
  unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
  select(hp) %>%
  mutate("{substr(i, 1, 9)}_c" := rep(1, length(.))) -> study_list[[i]]


# put a 1 to indicate the presence of the interaction (rows) in each study (columns)
# it is important that the overall network (or blood_all) has to be the first one in study_list
# Then left_join joins based on the interactions in blood_all, thus, maintaining our requirement 
# of using blood_all as the scaffold to make the core network
edge_counter <- list(study_list) %>% 
  reduce(left_join, by = "hp")

# adds up the number of datasets an interaction is found in
edge_rowsums <- rowSums(edge_counter[,2:ncol(edge_counter)], na.rm = T)
edge_counter$edge_rowsums <- edge_rowsums

#####################################################################################################

# finding edges present in multiple datasets
e_1 <- edge_counter[edge_counter$edge_rowsums == 1,] #present only in the overall network
e_4 <- edge_counter[edge_counter$edge_rowsums == 4,] #present only in the overall+3 networks
e_5 <- edge_counter[edge_counter$edge_rowsums == 5,] #present only in the overall+4 networks
e_6 <- edge_counter[edge_counter$edge_rowsums == 6,] #present only in the overall+5 networks
e_7 <- edge_counter[edge_counter$edge_rowsums == 7,] #present only in the overall+6 networks

# a sampling of the edges present, where possible, to make a representative network
e <- rbind(e_1[sample(200, 200),], e_5[sample(200, 200),], e_6[sample(200, 200),], e_7)

# separate the names of the genes in the interaction pair joined by underscores
hp <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', e$hp), ' ')

h <- sapply(hp, function(x) x[[1]])
p <- sapply(hp, function(x) x[[2]])
edges_df<- data.frame(h = h, p = p, e[,2:ncol(e)])
ig <- graph_from_data_frame(edges_df[,c(1,2,ncol(edges_df))], directed = F)

E(ig)$weight <- e$edge_rowsums
V(ig)$color = c(rep("grey", length(unique(edges_df$h))), rep("tomato", length(unique(edges_df$p))))

E(ig)$color[E(ig)$weight == 1] <- 'yellow'
E(ig)$color[E(ig)$weight == 5] <- 'green'
E(ig)$color[E(ig)$weight == 6] <- 'navy'
E(ig)$color[E(ig)$weight == 7] <- 'red'

svg("blood_edge_counter_sampled.svg", width = 50, height = 50)
plot(ig, vertex.color = V(ig)$color,
     edge.color = E(ig)$color,
     vertex.size = 2, vertex.label=NA,
     layout = layout_with_graphopt,
     edge.width = 8,
     edge.curved=TRUE,)
dev.off()

#####################################################################################################

############ Making the core network ################
# Core network definition: the edge has to be present in the overall network, 
# in at least 1 human network and in at least 1 network from another host organism
# We ensure the presence in the overall network in the construction of "edge_counter": 
# all the edges there are definitely present in the overall network

# to make sure that an edge in the core network is present in at least one human dataset
human_edges <- edge_counter[(!is.na(edge_counter$ERP106451_c) | !is.na(edge_counter$DRP000987_c) | 
  !is.na(edge_counter$hpv_c) | !is.na(edge_counter$SRP032775_c) | !is.na(edge_counter$SRP233153_c) | !is.na(edge_counter$ERP023982_c)),]

# to make sure that an edge in the core network is present in at least one other host organism dataset - 
# so that the edge is not just present in the human dataset
core_edges <- human_edges[(!is.na(human_edges$ERP110375_c) | !is.na(human_edges$ERP004598_c) | !is.na(human_edges$mbl_c) | 
!is.na(human_edges$SRP118827_c) | !is.na(human_edges$SRP116793_c) | !is.na(human_edges$SRP116593_c) | !is.na(human_edges$SRP118996_c) | 
!is.na(human_edges$SRP108356_c) | !is.na(human_edges$SRP118503_c)),]

# separate the names of the genes in the interaction pair joined by underscores
hp <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', core_edges$hp), ' ')

h <- sapply(hp, function(x) x[[1]])
p <- sapply(hp, function(x) x[[2]])
blood_core_network <- data.frame(host = h, parasite = p, weight = core_edges$edge_rowsums)

saveRDS(blood_core_network, file = "blood_core_network.rds") # with orthogroups

#### convert orthogroups to ensembl and PlasmoDB IDs
pOG <- read.delim("parasite_orthogroups.txt", stringsAsFactors=FALSE)
hOG <- read.delim("host_orthogroups.txt", stringsAsFactors=FALSE)

host <- sapply(h, function(x) hOG[hOG$Orthogroup %in% x, "human"])
parasite <- sapply(p, function(x) pOG[pOG$Orthogroup %in% x, "Pberghei"])

blood_core_network_IDs <- data.frame(host = host, parasite = parasite, weight = core_edges$edge_rowsums)
saveRDS(blood_core_network_IDs, file = "blood_core_network_IDs.rds") # with gene IDs
write.csv(blood_core_network_IDs, "blood_core_network_IDs.csv", row.names = F, quote = F)