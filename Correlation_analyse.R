# function to assign an RData object to a name
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

library(WGCNA)
library(reshape2)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(grid)
library(stringr)

########### Main operations ###########

studyID <- "blood_overall"

# make sum of files with uncorrected pvalues (permutation scores) of all the permutation tests for blood overall dataset
# to get total permutation score of a gene pair for further analysis
sum <- matrix(rep(0, 17991*4005), nrow = 17991) #for hpxpp analysis - bipartite and parasite-parasite edges

files <- grep(pattern = "outer_blood_overall_", list.files())
file.names <- list.files()[files]

for(i in 1:length(file.names))
{
  print(i)
  if(grepl(pattern = ".RData", file.names[i]))
  {
    load(file.names[i])
    sum <- sum + outer 
  }
}

hp <- t(loadRData("blood.ortho.data.RData"))
p <- hp[13987:17991,]

ori_cor <- cor(t(hp), t(p), use = "pairwise.complete.obs") # original correlation coefficients - test statistic

# bipartite part of the correlation matrix
hp_ori_cor <- ori_cor[1:13986,] 

#  correlation matrix to data.frame
cor_melt <- melt(hpori_cor)
colnames(cor_melt) <- c("gene1", "gene2", "cor")

# uncoreected pvalue (permutation scores) matrix to data.frame
pval <- sum[1:13986,]
pval_melt <- melt(pval)
colnames(pval_melt) <- c("gene1", "gene2", "permute_score")

# merging correlation data frame with permutation score data frame
pval_cor <- cbind(cor_melt, pval_melt$permute_score)
colnames(pval_cor)[4] <- "permute_score"

# remove possible NAs from the data frame
pval_cor_na.omit <- na.omit(pval_cor)
                
# get gene pairs with permutation score (uncorrected pvalues) 0 - these are the pairs we use for further analysis
pval0 <- pval_cor_na.omit[pval_cor_na.omit$permute_score==0,]

# find number of unuique host genes in the result
unique_hg <- unique(as.character(pval0[grep(pattern = "h_OG", pval0$gene2),1]))

# find number of unuique parasite genes in the result
unique_pg <- unique(as.character(pval0[grep(pattern = "p_OG", pval0$gene2),2]))

# data frame of all bipartite edges - should be the same as pval0
bipartite <- pval0

colnames(bipartite) <- sapply(colnames(bipartite), function(x) paste0(studyID,"_",type,"_",x, collapse = ''))
saveRDS(bipartite, file = paste0(studyID,"_", type, "_bipartite.rds", collapse = ''))

################### parasite-parasite part of the correlation matrix
p_ori_cor <- ori_cor[13987:17991,]
pcor_melt <- melt(p_ori_cor)
colnames(pcor_melt) <- c("gene1", "gene2", "cor")

ppval <- sum[13987:17991,]
ppval_melt <- melt(ppval)
colnames(ppval_melt) <- c("gene1", "gene2", "permute_score")

ppval_cor <- cbind(pcor_melt, ppval_melt$permute_score)
colnames(ppval_cor)[4] <- "permute_score"

ppval_cor_na.omit <- na.omit(ppval_cor)

ppval0 <- ppval_cor_na.omit[ppval_cor_na.omit$permute_score==0,]

para <- ppval0
colnames(para) <- sapply(colnames(para), function(x) paste0(studyID,"_",type,"_",x, collapse = ''))
saveRDS(para, file =  paste0(studyID,"_", type, "_para.rds", collapse = ''))

##############################################################################################################

########### Permutation score (uncorrected pvalue) vs Correlation coefficients #########

# there are about 56 million data points (13000 * 4005) pertaining to this relationship
# because there are a total of 56 million gene pairs that we are comparing.
# We group the permutation scores into two groups: Group 0 and Group non-zero (greater then zero)
# we plot the density of these groups across the range of correlation coefficients

# to plot from one study:
colnames(bipartite) <- c("gene1", "gene2", "Correlation", "Score")

cor_coefs <- rbind(bipartite[,c(3,4)])
colnames(cor_coefs) <- c("Correlation", "Score")

cor_coefs$Group <- ifelse(cor_coefs$Score == 0, "0", "1")
save(cor_coefs, file = "Cor_distribution.RData")

p<-ggplot(cor_coefs, aes(x=Correlation, fill=Group)) +
  geom_density(alpha=0.4) +
  theme_bw()+
  xlab("Correlation coefficient [-1,1]") +
  theme(text = element_text(size = 15))
ggsave(plot = p, "Cor_distribution.pdf")


# for the plot in the Mukherjee et al. 2020 manuscript:

#bl <- readRDS("blood_overall_all_bioartite.rds")
#mo <- readRDS("SRP118827_all_bioartite.rds")
#m <- readRDS("ERP004598_int_bioartite.rds")
#h <- readRDS("DRP000987_str_bioartite.rds")
#
#colnames(bl) <- c("gene1", "gene2", "Correlation", "Score")
#colnames(mo) <- c("gene1", "gene2", "Correlation", "Score")
#colnames(m) <- c("gene1", "gene2", "Correlation", "Score")
#colnames(h) <- c("gene1", "gene2", "Correlation", "Score")
#
#cor_coefs <- rbind(bl[,c(3,4)], mo[,c(3,4)], m[,c(3,4)], h[,c(3,4)])
#colnames(cor_coefs) <- c("Correlation", "Score")
#cor_coefs$Study <- c(rep("Overall", nrow(bl)), rep("SRP118827_all", nrow(mo)), rep("ERP004598_int", nrow(m)), rep("DRP000987_str", nrow(h)))
#
#cor_coefs$Group <- ifelse(cor_coefs$Score == 0, "0", "1")
#save(cor_coefs, file = "Cor_distribution.RData")
#
#
#p<-ggplot(cor_coefs, aes(x=Correlation, fill=Group)) +
#  geom_density(alpha=0.4) +
#  theme_bw() + 
#  facet_wrap(Study ~ ., scales = "free_y") + 
#  xlab("Correlation coefficient [-1,1]") +
#  theme(text = element_text(size = 15))
#ggsave(plot = p, "Cor_distribution.pdf")

################################################################################################################
 
############################ Reduction in the number of bipartite edges with permutation score (uncorrected pvalue) 0 with increase in number of permutations

hp <- t(loadRData("blood.ortho.data.RData"))
p <- hp[13987:17991,]

ori_cor <- cor(t(hp), t(p), use = "pairwise.complete.obs")
cor_melt <- melt(ori_cor)
pval_cor <- cbind(cor_melt, sum_melt$permute_score)
colnames(pval_cor)[4] <- "permute_score"

      
bipartite_edges_in_correlation_range <- function(pcc)
{
  # sum all permutation score files to get permutation scores of all gene pairs
  sum <- matrix(rep(0, 17991*4005), nrow = 17991)
  df <- data.frame(Perms = seq(10000, 100000, 10000))
  a = 1
  
  # after every 10000 permutations (max score from each scoring file) the number of edges are counted 
  for(i in 1:length(file.names))
  {
    print(i)
    if(grepl(pattern = ".RData", file.names[i]))
    {
      load(file.names[i])
      sum <- outer + sum
      sum_melt <- melt(sum)
      
      colnames(sum_melt) <- c("gene1", "gene2", "permute_score")
      
      pval_cor_na.omit <- na.omit(pval_cor)
      pval0 <- pval_cor_na.omit[pval_cor_na.omit$permute_score==0,]
      colnames(pval0) <- c("gene1", "gene2", "cor", "permute_score")
      
      if(pcc >= 0.9){ r1 = 0.9; r2 = 1.0 }
      if(pcc >= 0.8 & pcc < 0.9){ r1 = 0.8; r2 = 0.9 }
      if(pcc >= 0.7 & pcc < 0.8){ r1 = 0.7; r2 = 0.8 }
      if(pcc >= 0.6 & pcc < 0.7){ r1 = 0.6; r2 = 0.7 }
      if(pcc >= 0.5 & pcc < 0.6){ r1 = 0.5; r2 = 0.6 }
      if(pcc >= 0.4 & pcc < 0.5){ r1 = 0.4; r2 = 0.5 }
      if(pcc >= 0.3 & pcc < 0.4){ r1 = 0.3; r2 = 0.4 }
      if(pcc >= 0.2 & pcc < 0.3){ r1 = 0.2; r2 = 0.3 }
      if(pcc >= 0.1 & pcc < 0.2){ r1 = 0.1; r2 = 0.2 }
      all.edges <- pval0[(abs(pval0$cor) >= r1 & abs(pval0$cor) < r2),]
      hp.edges <- all.edges[((grepl(pattern = "h_OG", as.character(all.edges$gene1)) & grepl(pattern = "p_OG", as.character(all.edges$gene2))) |
       (grepl(pattern = "p_OG", as.character(all.edges$gene1)) & grepl(pattern = "h_OG", as.character(all.edges$gene2)))),]
      hp.edges.number <- nrow(hp.edges)
      
      df[a,2] <- pcc
      df[a,3] <- hp.edges.number
      
      a = a+1
    }
  }
  return(df)
}
bipartite_edges_cor_range <- rbind(bipartite_edges_in_correlation_range(0.9), 
                         bipartite_edges_in_correlation_range(0.8), 
                         bipartite_edges_in_correlation_range(0.7),
                         bipartite_edges_in_correlation_range(0.6),
                         bipartite_edges_in_correlation_range(0.5),
                         bipartite_edges_in_correlation_range(0.4),
                         bipartite_edges_in_correlation_range(0.3),
                         bipartite_edges_in_correlation_range(0.2),
                         bipartite_edges_in_correlation_range(0.1))

colnames(bipartite_edges_cor_range) <- c("Perms", "cor", "hp.edges")
df_melt <- melt(bipartite_edges_cor_range, id = c("Perms", "cor"))
colnames(df_melt) <- c("Perms", "cor", "edge.type", "count")

png("blood_overall_hp_edges_cor_range.png", width = 800, height = 500, units = "px")
edge.plot <- ggplot(df_melt, aes(x = Perms, y = count, colour = factor(cor))) +
  geom_line() + scale_x_continuous(breaks = seq(100000, 1000000, 100000), limits = c(100000, 1000000))+
  theme_bw() +
  theme(axis.text.x=element_blank()) +
  facet_wrap(. ~ edge.type, scales="free_y")
dev.off()

#######################################################################################################################

########### script to make annotation table for heatmaps #########

# allHPexp <- read.delim("allHPexp.txt", sep = ',', stringsAsFactors = F)

# datasets <- c("DRP000987.ortho.data.all", "DRP000987.ortho.data.int", "DRP000987.ortho.data.str", 
#               "ERP106451.ortho.data.all", "ERP106451.ortho.data.int", "ERP106451.ortho.data.str",
#               "ERP023982.ortho.data.all", "ERP023982.ortho.data.int",
#               "ERP004598.ortho.data.all", "ERP004598.ortho.data.int", "ERP004598.ortho.data.str",
#               "ERP110375.ortho.data.all", "ERP110375.ortho.data.int", "ERP110375.ortho.data.str",
#               "SRP118996.ortho.data.all", "SRP118996.ortho.data.int", "SRP118996.ortho.data.str",
#               "SRP118827.ortho.data.all", "SRP118827.ortho.data.int", "SRP118827.ortho.data.str",
#               "SRP116593.ortho.data.all", "SRP116593.ortho.data.int", "SRP116593.ortho.data.str",
#               "SRP116793.ortho.data.all", "SRP116793.ortho.data.int", "SRP116793.ortho.data.str",
#               "SRP032775.ortho.data.all", "SRP032775.ortho.data.int", "SRP032775.ortho.data.str",
#               "SRP108356.ortho.data.all", "SRP108356.ortho.data.int", "SRP108356.ortho.data.str",
#               "SRP118503.ortho.data.all", "SRP118503.ortho.data.int", "SRP233153.ortho.data.int",
#               "hpv_bl.ortho.data.int", "m_bl.ortho.data.int"
#               )

# rn <- c("DRP000987_all", "DRP000987_int", "DRP000987_str", 
#         "ERP106451_all", "ERP106451_int", "ERP106451_str",
#         "ERP023982_all", "ERP023982_int", 
#         "ERP004598_all", "ERP004598_int", "ERP004598_str",
#         "ERP110375_all", "ERP110375_int", "ERP110375_str",
#         "SRP118996_all", "SRP118996_int", "SRP118996_str",
#         "SRP118827_all", "SRP118827_int", "SRP118827_str",
#         "SRP116593_all", "SRP116593_int", "SRP116593_str",
#         "SRP116793_all", "SRP116793_int", "SRP116793_str",
#         "SRP032775_all", "SRP032775_int", "SRP032775_str",
#         "SRP108356_all", "SRP108356_int", "SRP108356_str",
#         "SRP118503_all", "SRP118503_int", "SRP233153_int",
#         "hpv_bl_int", "m_bl_int"
#         )

# anno <- data.frame()

# for(i in 1:length(datasets))
# {
#   ds <- loadRData(paste0(str_sub(rn[i], 1, -5), "/cor/", datasets[i],".RData", collapse = ))
#   runs <- sapply(colnames(ds), function(x) strsplit(x, split = "_")[[1]][1])
  
#   pp <- sapply(runs, function(x) allHPexp[which(allHPexp$RunID==x),"ParaGenomeProp"])
#   pp_median <- median(pp)
#   hh <- sapply(runs, function(x) allHPexp[which(allHPexp$RunID==x),"HostGenomeProp"])
#   hh_median <- median(hh)

#   host <- as.character(allHPexp[which(allHPexp$RunID==runs[1]),"Host"])
#   parasite <- as.character(allHPexp[which(allHPexp$RunID==runs[1]),"Parasite"])
  
#   anno[i,1] <- datasets[i]
#   anno[i,2] <- host
#   anno[i,3] <- parasite
#   anno[i,4] <- pp_median
#   anno[i,5] <- hh_median
#   anno[i,6] <- rn[i]
# }

# rownames(anno) <- anno[,6]
# anno <- anno[,c(2,3,4,5)]
# colnames(anno) <- c("Host", "Parasite", "Median_ParaTransProp", "Median_HostTransProp")

# write.table(anno, "anno.txt", sep = '\t', row.names = T)
# saveRDS(anno, file = "anno.rds")


########### Function to produce comparison matrices #########

# studyID <- c("ERP106451_str", "ERP106451_int", "ERP106451_all",
#  "SRP118996_str", "SRP118996_int", "SRP118996_all",
#  "SRP118827_str", "SRP118827_int", "SRP118827_all",
#  "SRP116793_str", "SRP116793_int", "SRP116793_all",
#  "SRP116593_str", "SRP116593_int", "SRP116593_all",
#  "DRP000987_str", "DRP000987_int", "DRP000987_all",
#  "ERP023982_int", "ERP023982_all",
#  "ERP004598_str", "ERP004598_int", "ERP004598_all",
#  "ERP110375_str", "ERP110375_int", "ERP110375_all",
#  "SRP032775_all", "SRP032775_int", "SRP032775_str",
#  "SRP108356_all", "SRP108356_int", "SRP108356_str",
#  "SRP118503_all", "SRP118503_int", "SRP233153_int",
#  "hpv_bl_int", "m_bl_int")

# datasets <- c()
# list_ds <- list()
# l <- 0
# for(i in 1:length(studyID))
# {
#   ds <- loadRData(paste0(str_sub(studyID[i], 1, -5), "/cor/", studyID[i],"_bipartite.RData"))
#   colnames(ds) <- sapply(colnames(ds), function(x) paste0(studyID[i],"_", feature, "_", x))
    
#   list_ds[[paste0(studyID[i], "_bipartite_edges")]] <- ds
#   l <- length(list_ds)
# }

# data <- list_ds
# feature = "b_edges" # bipartite edges
# operations <- function(data, op, feature)
# {

#  # size of intersecting interactions (common interactions) between all dataset pairs
#  if(op == "overlap")
#   {
#         feature = feature
#         m <- length(data)
#         n <- ncol(data[[1]])
        
#         upset <- list()
#         mat <- matrix(rep(0, (length(data)^2)), nrow = length(data), ncol = length(data))
         
#         if(n>=3)
#         {
#           # finding intersecting edges
#           for(i in 1:length(data))
#           {
#             upset[[i]] <- data[[i]][which(data[[i]][,4] == 0),1:2] #here 4 because cor column was included
#             upset[[i]] <- paste(upset[[i]][,1], upset[[i]][,2], sep = '_')
#           }
#           names(upset) <- names(data)
#           for(i in 1:length(upset))
#             for(j in 1:length(upset))
#               mat[i,j] <- 2*length(intersect(upset[[i]], upset[[j]]))/(length(upset[[i]]) + length(upset[[j]]))
#           rownames(mat) <- colnames(mat) <- names(upset)
#           save(mat, file = paste0(feature, "b_overlap_matrix_all_datasets_all_liver.RData"))
#           save(upset, file =  paste0(feature,"b_overlap_upset_all_datasets_all_liver.RData"))
#           write.table(mat,  paste0(feature, "b_overlap_matrix_all_datasets_all_liver.txt"), sep = '\t', row.names = T)
          
#           mat1 <- matrix(rep(0, (length(data)^2)), nrow = length(data), ncol = length(data))
          
#           for(i in 1:length(upset))
#             for(j in 1:length(upset))
#               mat1[i,j] <- length(intersect(upset[[i]], upset[[j]]))
#           rownames(mat1) <- colnames(mat1) <- names(upset)
#           save(mat1, file = paste0(feature, "b_overlap_raw_numbers_matrix_all_datasets_all_liver.RData"))
#           write.table(mat1,  paste0(feature, "b_overlap_raw_numbers_matrix_all_datasets_all_liver.txt"), sep = '\t', row.names = T)
#         }
    
#     if(n==1)
#     {
#       # finding intersecting genes
#       upset <- data
#       names(upset) <- names(data)
#       for(i in 1:length(data))
#         for(j in 1:length(data))
#           mat[i,j] <- 2*length(intersect(data[[i]][,1], data[[j]][,1]))/(nrow(data[[i]]) + nrow(data[[j]]))
      
#       rownames(mat) <- colnames(mat) <- names(upset)
#       save(mat, file = paste0(feature, "_overlap_matrix_all_datasets_liver.RData"))
#       save(upset, file =  paste0(feature,"_overlap_upset_all_datasets_liver.RData"))
#       write.table(mat,  paste0(feature, "_overlap_matrix_all_datasets_liver.txt"), sep = '\t', row.names = T)
#      }
#     res <- list(upset, mat)
#   }
  
#   ###### significance of the size of intersections ###########

#   if(op == "sig")
#   {
#     feature = feature
#     m <- length(data)
#     n <- ncol(data[[1]])
    
#     sig <- data.frame()
#     a = 1
#     upset <- list()
    
#     # finding intersecting edges
#     for(i in 1:length(data))
#     {
#       upset[[i]] <- data[[i]][which(data[[i]][,4] == 0),1:2] # edges with ISIGEM score 0
#       upset[[i]] <- paste(upset[[i]][,1], upset[[i]][,2], sep = '_')
#     }
#     names(upset) <- names(data)
#     system.time(
#     for(i in 1:length(upset))
#     {
#       #U1 = nrow(data[[i]])#[which(data[[i]][,4] != 100000),])
#       for(j in 1:length(upset))
#       {
#         #U2 = nrow(data[[j]])#[which(data[[j]][,4] != 100000),])
#         #U = U1 + U2
#         U = 13981*4005
#         X = length(union(upset[[i]], upset[[j]]))
#         Y = length(setdiff(upset[[i]],upset[[j]]))
#         W = length(setdiff(upset[[j]],upset[[i]]))
#         Z = length(intersect(upset[[i]], upset[[j]]))
        
#         significance <- fisher.test(matrix(c(U-X, Y, W, Z), nrow=2), alternative = "greater")
        
#         sig[a,1] <- names(upset)[i]
#         sig[a,2] <- length(upset[[i]])
#         sig[a,3] <- names(upset)[j]
#         sig[a,4] <- length(upset[[j]])
#         sig[a,5] <- U
#         sig[a,6] <- Z
#         sig[a,7] <- significance$p.value
        
#         a = a+1
#       }
        
#     })
    
#     colnames(sig) <- c("Set1","Set1_edges#", "Set2", "Set2_edges#", "Universe_size", "Intersection_size", "pvalue")
#     padj <- p.adjust(sig[,7], method = "BH")
#     sig$p.adj <- padj
    
#     #cast the df into a matrix
#     sig.mat <- reshape(sig[,c(1,3,8)], timevar = "Set1", idvar = "Set2", direction = "wide")
#     rownames(sig.mat) <- sig.mat[,1]
#     sig.mat <- sig.mat[,-1]
#     colnames(sig.mat) <- str_sub(colnames(sig.mat),1, -9)
#     rownames(sig.mat) <- str_sub(rownames(sig.mat), 1, -9)
#     logt.sig.mat <- log10(1/(sig.mat + 1e-5))
    
#     save(sig, file = paste0(feature, "_overlap_signif_all_liver.RData"))
#     write.table(sig,  paste0(feature, "_overlap_signif_all_liver.txt"), sep = '\t', row.names = T)
#     save(sig.mat, file = paste0(feature, "_overlap_signif_all_liver_matrix.RData"))
#     write.table(sig.mat,  paste0(feature, "_overlap_signif_all_liver_matrix.txt"), sep = '\t', row.names = T)
#     save(logt.sig.mat, file = paste0(feature, "_logt_overlap_signif_all_liver_matrix.RData"))
    
#   }
#   return(res)
# }

# Cross_study_comparison(feature = "b_edges", op = "cor")

# colnames(mat1) <- sapply(colnames(mat1), function(x) str_sub(x, 1, -9))
# rownames(mat1) <- sapply(rownames(mat1), function(x) str_sub(x, 1, -9))

# pdf("pheatmap_b_edges_overlap.pdf" ) #, width = 20, height = 20, unit = "cm", res = 300
# pheatmap::pheatmap(log10(mat1+1), annotation_row = anno, 
#                    fontsize = 8, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
#                    main = "Intersection size of bipartite edges",
#                    cluster_cols = F, cluster_rows = T) #-log10(1/(sig.mat + 1e-06))
# dev.off()

# pdf("pheatmap_b_edges_significance.pdf")
# pheatmap::pheatmap(sig.mat, annotation_row = anno, 
#                    fontsize = 8, color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100),
#                    main = "Significance of intersection size of bipartite edges",
#                    cluster_cols = F, cluster_rows = F) #-log10(1/(sig.mat + 1e-06))
# dev.off()

# ##### Computing and plotting Jaccard index #############

# Jaccard <- function(edgelist)
# {
#   d <- length(edgelist)
#   jac.sim <- matrix(rep(0, d^2), nrow = d)
  
#   for(i in 1:d)
#   {
#     for(j in 1:d)
#     {
#       A <- edgelist[[i]]
#       B <- edgelist[[j]]

#       numerator <- length(intersect(A, B))
#       denominator <- length(union(A, B))

#       jac.ind <- numerator/denominator

#       jac.sim[i,j] <- jac.ind
#     }
#   }
#   colnames(jac.sim) <- names(edgelist)
#   rownames(jac.sim) <- names(edgelist)

#   return(jac.sim)
# }

# jaccard.index <- Jaccard(upset)

# jac <- -log10(jaccard.index) # or the -log10 of the jaccard indices
# colnames(jac) <- sapply(colnames(jac), function(x) str_sub(x, 1, -9))
# rownames(jac) <- sapply(rownames(jac), function(x) str_sub(x, 1, -9))

# is.na(jac)<-sapply(jac, is.infinite)
# jac[is.na(jac)]<-8

# pdf("pheatmap_b_edges_jaccard.pdf")
# pheatmap::pheatmap(jac, annotation_row = anno, 
#                    fontsize = 8, color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100),
#                    main = "Intersection size of bipartite edges - -log10(Jaccard Index)",
#                    cluster_cols = T, cluster_rows = T) #-log10(1/(sig.mat + 1e-06))
# dev.off()