# # Script to find out which bipartite interactions are present across multiple datasets
# # and to create the core network

# library(tidyverse)
# library(igraph)

# loadRData <- function(fileName){
#   #loads an RData file, and returns it
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }

# # overall network
# all <- loadRData("overall_addblood/cor/blood_all_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(all_c = rep(1, length(.)))

# # human studies
# DRP000987 <- loadRData("DRP000987/cor/DRP000987_str_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(DRP000987_c = rep(1, length(.)))
# ERP106451 <- loadRData("ERP106451/cor/ERP106451_all_bipartite.RData")%>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(ERP106451_c = rep(1, length(.)))
# SRP032775 <- loadRData("SRP032775/cor/SRP032775_str_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP032775_c = rep(1, length(.)))
# SRP233153 <- loadRData("SRP233153/cor/SRP233153_int_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP233153_c = rep(1, length(.)))
# ERP023982 <- loadRData("ERP023982/cor/ERP023982_int_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(ERP023982_c = rep(1, length(.)))
# hpv <- loadRData("hpv_bl/cor/hpv_bl_int_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(hpv_c = rep(1, length(.)))

# # mouse studies
# ERP110375 <- loadRData("ERP110375/cor/ERP110375_str_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(ERP110375_c = rep(1, length(.)))
# ERP004598 <- loadRData("ERP004598/cor/ERP004598_int_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(ERP004598_c = rep(1, length(.)))
# mbl <- loadRData("m_bl/cor/m_bl_int_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(mbl_c = rep(1, length(.)))


# # monkey studies
# SRP118827 <- loadRData("SRP118827/cor/SRP118827_all_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP118827_c = rep(1, length(.)))
# SRP116793 <- loadRData("SRP116793/cor/SRP116793_all_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP116793_c = rep(1, length(.)))
# SRP118996 <- loadRData("SRP118996/cor/SRP118996_all_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP118996_c = rep(1, length(.)))
# SRP116593 <- loadRData("SRP116593/cor/SRP116593_int_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP116593_c = rep(1, length(.)))
# SRP108356 <- loadRData("SRP108356/cor/SRP108356_str_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP108356_c = rep(1, length(.)))
# SRP118503 <- loadRData("SRP118503/cor/SRP118503_int_bipartite.RData") %>%
#   unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
#   select(hp) %>%
#   mutate(SRP118503_c = rep(1, length(.)))

# # put a 1 to indicate the presence of the interaction (rows) in each study (columns)
# edge_counter <- list(all, DRP000987, ERP106451, SRP233153, SRP032775, ERP023982, # human studies
#   ERP110375, ERP004598, # mouse studies
#   SRP116793, SRP118827, SRP118996, SRP116593, SRP108356, SRP118503, # monkey studies
#   hpv, mbl) %>% reduce(left_join, by = "hp")

# # adds up the number of datasets an interaction is found in
# edge_rowsums <- rowSums(edge_counter_para[,2:ncol(edge_counter_para)], na.rm = T)
# edge_counter$edge_rowsums <- edge_rowsums

# #####################################################################################################

# ############ Making the core network ################
# # Core network definition: the edge has to be present in the overall network, in at least 1 human network and in at least 1 network from another host organism
# # We ensure the presence in the overall network in the construction of "edge_counter": all the edges there are definitely present in the overall network

# # to make sure that an edge in the core network is present in at least one human dataset
# human_edges <- edge_counter_para[(!is.na(edge_counter_para$ERP106451_c) | !is.na(edge_counter_para$DRP000987_c) | 
#   !is.na(edge_counter_para$hpv_c) | !is.na(edge_counter_para$SRP032775_c) | !is.na(edge_counter_para$SRP233153_c) | !is.na(edge_counter_para$ERP023982_c)),]

# # to make sure that an edge in the core network is present in at least one other host organism dataset - so that the edge is not just present in the human dataset
# core_edges <- human_edges[(!is.na(human_edges$ERP110375_c) | !is.na(human_edges$ERP004598_c) |
# !is.na(human_edges$mbl_c) | 
# !is.na(human_edges$SRP118827_c) | !is.na(human_edges$SRP116793_c) | !is.na(human_edges$SRP116593_c) | !is.na(human_edges$SRP118996_c) | !is.na(human_edges$SRP108356_c) | !is.na(human_edges$SRP118503_c)),]

# blood_core_edges <- core_edges

# # finding edges present in multiple datasets
# e_1 <- blood_edge_counter[blood_edge_counter$edge_rowsums == 1,] #present only in the overall network
# e_4 <- blood_edge_counter[blood_edge_counter$edge_rowsums == 4,] #present only in the overall+3 networks
# e_5 <- blood_edge_counter[blood_edge_counter$edge_rowsums == 5,] #present only in the overall+4 networks
# e_6 <- blood_edge_counter[blood_edge_counter$edge_rowsums == 6,] #present only in the overall+5 networks
# e_7 <- blood_edge_counter[blood_edge_counter$edge_rowsums == 7,] #present only in the overall+6 networks

# # a sampling of the edges present, where possible, to make a representative network
# e <- rbind(e_1[sample(200, 200),], e_5[sample(200, 200),], e_6[sample(200, 200),], e_7)

# hp <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', e$hp), ' ')

# h <- sapply(hp, function(x) x[[1]])
# p <- sapply(hp, function(x) x[[2]])
# edges_df<- data.frame(h = h, p = p, e[,2:ncol(e)])
# ig <- graph_from_data_frame(edges_df[,c(1,2,ncol(edges_df))], directed = F)

# E(ig)$weight <- e$edge_rowsums
# V(ig)$color = c(rep("grey", length(unique(edges_df$h))), rep("tomato", length(unique(edges_df$p))))

# E(ig)$color[E(ig)$weight == 1] <- 'yellow'
# E(ig)$color[E(ig)$weight == 5] <- 'green'
# E(ig)$color[E(ig)$weight == 6] <- 'navy'
# E(ig)$color[E(ig)$weight == 7] <- 'red'

# svg("blood_edge_counter_sampled.svg", width = 50, height = 50)
# plot(ig, vertex.color = V(ig)$color,
#      edge.color = E(ig)$color,
#      vertex.size = 2, vertex.label=NA,
#      layout = layout_with_graphopt,
#      edge.width = 8,
#      edge.curved=TRUE,)
# dev.off()