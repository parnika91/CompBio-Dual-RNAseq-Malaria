# Script to construct beta regression models with the properties of gene co-expression


library(dplyr)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)
library(stringr)
library(stargazer)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

options(echo=TRUE)
args <- commandArgs(TRUE)
dataset <- args[1] # enter overall or core

# function to get the core network of parasite-parasite edges 
blood_core_para_edges <- function()
{
  studies <- c("blood_all",
  # human studies
  "DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int", #hpv_bl is the collection of human - P.vivax blood runs
  # mouse studies
  "ERP110375_str", "ERP004598_int", "m_bl_int", # m_bl is a collection of mouse - Plasmodium studies
  # monkey studies
  "SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int")

  for(i in studies)
    readRDS(paste0(substr(i, 1, 9), "/cor/", i, "_para.rds", collapse = "")) %>%
    unite("pp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
    select(pp) %>%
    mutate("{substr(i, 1, 9)}_c" := rep(1, length(.))) -> study_list[[i]]


  # put a 1 to indicate the presence of the interaction (rows) in each study (columns)
  # it is important that the overall network (or blood_all) has to be the first one in study_list
  # Then left_join joins based on the interactions in blood_all, thus, maintaining our requirement
  # of using blood_all as the scaffold to make the core network
  edge_counter <- list(study_list) %>%
    reduce(left_join, by = "pp")

  # adds up the number of datasets an interaction is found in
  edge_rowsums <- rowSums(edge_counter[,2:ncol(edge_counter)], na.rm = T)
  edge_counter$edge_rowsums <- edge_rowsums

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
  blood_core_para <- data.frame(host = h, parasite = p, weight = core_edges$edge_rowsums)

  return(blood_core_para) # with orthogroups
}

if(dataset =="overall")
{
  data <- readRDS("blood_all_para.rds")
}else
{
  #data <- blood_core_para_edges()
  data <- readRDS("blood_core_para.rds")
}


d <- data.frame(gene1 = as.character(data[,1]), gene2 = as.character(data[,2]))
ig <- graph_from_data_frame(d, directed = F)

dg <- degree(ig, v = V(ig), loops = F, normalized = F)
bw <- betweenness(ig, v = V(ig), directed = FALSE)
cl <- closeness(ig, vids = V(ig))
ec <- eigen_centrality(ig, directed = FALSE)
kc <- coreness(ig, mode = c("all"))

dg_p <- dg[grep(pattern = "p_OG", names(dg))]
dg_df <- as.data.frame(dg_p) %>%
  tibble::rownames_to_column("Orthogroup")

bw_p <- bw[grep(pattern = "p_OG", names(bw))]
bw_df <- as.data.frame(bw_p) %>%
  tibble::rownames_to_column("Orthogroup")

cl_p <- cl[grep(pattern = "p_OG", names(cl))]
cl_df <- as.data.frame(cl_p) %>%
  tibble::rownames_to_column("Orthogroup")

ec_p <- ec$vector[grep(pattern = "p_OG", names(ec$vector))]
ec_df <- as.data.frame(ec_p) %>%
  tibble::rownames_to_column("Orthogroup")

kc_p <- kc[grep(pattern = "p_OG", names(kc))]
kc_df <- as.data.frame(kc_p) %>%
  tibble::rownames_to_column("Orthogroup")

join_df <-plyr::join_all(list(dg_df, bw_df, cl_df,
                              ec_df, kc_df), by = "Orthogroup", type = "full")

join_df[is.na(join_df)] <- 0

colnames(join_df) <- c("Orthogroup",paste0(dataset, "_dg", collapse = ""),
                       paste0(dataset, "_bw", collapse = ""),
                       paste0(dataset, "_cl", collapse = ""),
                       paste0(dataset, "_ec", collapse = ""),
                       paste0(dataset, "_kc", collapse = ""))



RGR_MIS_phenotype = readRDS("RGR_MIS_phenotype.rds")
RGR_MIS_phenotype <- left_join(RGR_MIS_phenotype, join_df)
saveRDS(RGR_MIS_phenotype, file = "RGR_MIS_phenotype.rds")

### models ###

#RGR
dg_rgr <- betareg(data = RGR_MIS_phenotype, RGR ~ RGR_MIS_phenotype[,c(paste0(dataset, "_dg", collapse = ""))])
bw_rgr <- betareg(data = RGR_MIS_phenotype, RGR ~ RGR_MIS_phenotype[,c(paste0(dataset, "_bw", collapse = ""))])
ec_rgr <- betareg(data = RGR_MIS_phenotype, RGR ~ RGR_MIS_phenotype[,c(paste0(dataset, "_ec", collapse = ""))])

dg_bw_rgr <- betareg(data = RGR_MIS_phenotype, RGR ~ RGR_MIS_phenotype[,c(paste0(dataset, "_dg", collapse = ""))] + RGR_MIS_phenotype[,c(paste0(dataset, "_bw", collapse = ""))])
dg_ec_rgr <- betareg(data = RGR_MIS_phenotype, RGR ~ RGR_MIS_phenotype[,c(paste0(dataset, "_dg", collapse = ""))] + RGR_MIS_phenotype[,c(paste0(dataset, "_ec", collapse = ""))])

#MIS
dg_mis <- betareg(data = RGR_MIS_phenotype, MIS ~ RGR_MIS_phenotype[,c(paste0(dataset, "_dg", collapse = ""))])
ec_mis <- betareg(data = RGR_MIS_phenotype, MIS ~ RGR_MIS_phenotype[,c(paste0(dataset, "_ec", collapse = ""))])
bw_mis <- betareg(data = RGR_MIS_phenotype, MIS ~ RGR_MIS_phenotype[,c(paste0(dataset, "_bw", collapse = ""))])

dg_bw_mis <- betareg(data = RGR_MIS_phenotype, MIS ~ RGR_MIS_phenotype[,c(paste0(dataset, "_dg", collapse = ""))] + RGR_MIS_phenotype[,c(paste0(dataset, "_bw", collapse = ""))])
dg_ec_mis <- betareg(data = RGR_MIS_phenotype, MIS ~ RGR_MIS_phenotype[,c(paste0(dataset, "_dg", collapse = ""))] + RGR_MIS_phenotype[,c(paste0(dataset, "_ec", collapse = ""))])


models_list <- list(dg_rgr, ec_rgr, bw_rgr, dg_bw_rgr, dg_ec_rgr,
                    dg_mis, ec_mis, bw_mis, dg_bw_mis, dg_ec_mis)
names(models_list) <- c(paste0(dataset, "_dg_rgr"), paste0(dataset, "_ec_rgr"), paste0(dataset, "_bw_rgr"), paste0(dataset, "_dg_bw_rgr"), paste0(dataset, "_dg_ec_rgr"),
                        paste0(dataset, "_dg_mis"), paste0(dataset, "_ec_mis"), paste0(dataset, "_bw_mis"), paste0(dataset, "_dg_bw_mis"), paste0(dataset, "_dg_ec_mis"))


# making a table with model properties from the current models
model_table <- function(models)
{
  tab <- data.frame(model_name = rep("mn", 10),
    dg = rep("dg", 10),
    bw = rep("bw", 10),
    cl = rep("cl", 10),
    ec = rep("ec", 10),
    kc = rep("kc", 10),
    prsq = rep("prsq", 10))
  tab[] <- lapply(tab, as.character)

  for(i in 1:length(models))
  {
    tab[i,1] <- names(models[i])
    m <- models[i]
    coefs <- coef(summary(m[[1]]))$mean
    
    for(j in 2:nrow(coefs))
    {
      pred <- str_sub(rownames(coefs)[j], 41, 42)
      effect_size <- format(exp(coefs[j, 1]), nsmall = 3, scientific = T, digits = 5)
      pval <- format(coef(summary(m[[1]]))$mean[j,4], nsmall = 3, scientific = T, digits = 5)
      eff_pv <- paste(effect_size, pval, sep = "; ")

      tab[i,grep(pattern = pred, colnames(tab))] = eff_pv
    }

    tab[i,7] <- format(m[[1]][["pseudo.r.squared"]], nsmall = 3, scientific = T, digits = 5)
  }
  return(tab)
}

models_table <- model_table(models_list)

write.table(models_table, paste0(dataset, "_models_table_exp.txt", collapse = ""), sep = '\t', row.names = F)


## RGR plots ##
ggplot(data=RGR_MIS_phenotype, aes(x=RGR_MIS_phenotype[,c(paste0(dataset, "_ec", collapse = ""))], y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  geom_line(aes(y = predict(models_list[[paste0(dataset, "_ec_mis")]], RGR_MIS_phenotype))) +
  ggtitle("Relative growth rate vs eigen centrality") +
  xlab("Eigenvector centrality") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_EC.png")

## MIS plots ##
ggplot(data=RGR_MIS_phenotype, aes(x=RGR_MIS_phenotype[,c(paste0(dataset, "_ec", collapse = ""))], y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  geom_line(aes(y = predict(models_list[[paste0(dataset, "_ec_mis")]], RGR_MIS_phenotype))) +
  ggtitle("Mutagenesis index score vs eigen centrality") +
  xlab("Eigenvector centrality") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_EC.png")

stargazer(models_list, title="Results", align=TRUE, type = "text")
