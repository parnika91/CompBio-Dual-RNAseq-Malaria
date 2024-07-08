## Script to calculate network characteristics between blood and liver networks

# libraries
library(tidyverse)
library(igraph)
library(lsa)
library(bipartite)

blood_nw <- readRDS("blood_all_bipartite.rds") %>% 
  rename(gene1 = ends_with("gene1"), 
         gene2 = ends_with("gene2"), 
         cor = ends_with("cor")) %>% 
  select(gene1, gene2, cor)

liver_nw <- read.delim("liver_overall_bipartite.txt")%>% 
  rename(gene1 = Host, 
         gene2 = Parasite, 
         cor = Correlation_coef) %>% 
  select(gene1, gene2, cor)

blood_g <- graph_from_data_frame(blood_nw[,1:3], directed = F); blood_g$weights = blood_nw[,3]
liver_g <- graph_from_data_frame(liver_nw[,1:3], directed = F); liver_g$weights = liver_nw[,3]

bld_g <- blood_g; liv_g <- liver_g

#### degree and ec distribution by host and parasite

blood_host_degree <- data.frame(degree = igraph::degree(bld_g, v = V(bld_g)[1:12652], normalized = T, loops = F, mode = "all"))
blood_para_degree <- data.frame(degree = igraph::degree(bld_g, v = V(bld_g)[12653:16648], normalized = T, loops = F, mode = "all"))

liver_host_degree_vector <- unlist(liver_host_degree)
liver_para_degree_vector <- unlist(liver_para_degree)

blood_host_degree_vector <- unlist(blood_host_degree)
blood_para_degree_vector <- unlist(blood_para_degree)

wilcox.test(x = blood_host_degree_vector, y = liver_host_degree_vector, paired = F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  blood_host_degree_vector and liver_host_degree_vector
# W = 58207685, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(x = blood_para_degree_vector, y = liver_para_degree_vector, paired = F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  blood_para_degree_vector and liver_para_degree_vector
# W = 12038000, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

### eigen centrality

blood_host_eigen <- data.frame(EC = bld_eigen$EC[1:12652])
blood_para_eigen <- data.frame(EC = bld_eigen$EC[12653:16648])

# liver

liver_host_eigen <- data.frame(EC = liv_eigen$EC[1:5810])
liver_para_eigen <- data.frame(EC = liv_eigen$EC[5811:9506])

liver_host_eigen_vector <- unlist(liver_host_eigen)
liver_para_eigen_vector <- unlist(liver_para_eigen)

blood_host_eigen_vector <- unlist(blood_host_eigen)
blood_para_eigen_vector <- unlist(blood_para_eigen)

wilcox.test(x = blood_para_eigen_vector, y = liver_para_eigen_vector, paired = F)
# Wilcoxon rank sum test with continuity correction
# 
# data:  blood_para_eigen_vector and liver_para_eigen_vector
# W = 10932622, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(x = blood_host_eigen_vector, y = liver_host_eigen_vector, paired = F)

# Wilcoxon rank sum test with continuity correction
# 
# data:  blood_host_eigen_vector and liver_host_eigen_vector
# W = 58168886, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0



## plot same organism for different networks together
blood_host_degree$type = rep("blood_host_degree", nrow(blood_host_degree))
liver_host_degree$type = rep("liver_host_degree", nrow(liver_host_degree))

blood_liver_host_degree <- rbind(blood_host_degree, liver_host_degree)                                 
ggplot(blood_liver_host_degree, aes(degree, fill = type)) +
  geom_density(alpha = 0.5, size = 0.08) +
  scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  labs(fill = "Network") +
  xlab("Degree distribution") +
  ylab("Density") +
  theme(text = element_text(size = 12)) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("darkgoldenrod1", "darkolivegreen3")) +
  ggtitle("Degree distribution of host genes in blood and liver networks")
ggsave("blood_liver_host_degree_dist.png", width = 20, height = 10, units = "cm", dpi = 300)
#---
blood_para_degree$type = rep("blood_para_degree", nrow(blood_para_degree))
liver_para_degree$type = rep("liver_para_degree", nrow(liver_para_degree))

blood_liver_para_degree <- rbind(blood_para_degree, liver_para_degree)                                 
ggplot(blood_liver_para_degree, aes(degree, fill = type)) +
  geom_density(alpha = 0.5, size = 0.08) +
  scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  labs(fill = "Network") +
  xlab("Degree distribution") +
  ylab("Density") +
  theme(text = element_text(size = 12)) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("darkgoldenrod1", "darkolivegreen3")) +
  ggtitle("Degree distribution of parasite genes in blood and liver networks")
ggsave("blood_liver_para_degree_dist.png", width = 20, height = 10, units = "cm", dpi = 300)

#--eigen EC
blood_host_eigen$type = rep("blood_host_eigen", nrow(blood_host_eigen))
liver_host_eigen$type = rep("liver_host_eigen", nrow(liver_host_eigen))

blood_liver_host_eigen <- rbind(blood_host_eigen, liver_host_eigen)                                 
ggplot(blood_liver_host_eigen, aes(EC, fill = type)) +
  geom_density(alpha = 0.5, size = 0.08) +
  scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  labs(fill = "Network") +
  xlab("EC distribution") +
  ylab("Density") +
  theme(text = element_text(size = 12)) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("darkorchid3", "burlywood1")) +
  ggtitle("Eigenvector centrality (EC) distribution of host genes\nin blood and liver networks")
ggsave("blood_liver_host_eigen_dist.png", width = 20, height = 10, units = "cm", dpi = 300)
#---
blood_para_eigen$type = rep("blood_para_eigen", nrow(blood_para_eigen))
liver_para_eigen$type = rep("liver_para_eigen", nrow(liver_para_eigen))

blood_liver_para_eigen <- rbind(blood_para_eigen, liver_para_eigen)                                 
ggplot(blood_liver_para_eigen, aes(EC, fill = type)) +
  geom_density(alpha = 0.5, size = 0.08) +
  scale_x_log10() +
  #scale_y_log10() +
  theme_bw() +
  labs(fill = "Network") +
  xlab("EC distribution") +
  ylab("Density") +
  theme(text = element_text(size = 12)) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("darkorchid3", "burlywood1")) +
  ggtitle("Eigenvector centrality (EC) distribution of parasite genes\nin blood and liver networks")
ggsave("blood_liver_para_eigen_dist.png", width = 20, height = 10, units = "cm", dpi = 300)



### compute nw properties for blood and liver nws

liver_nw_wide <- liver_nw %>% 
  pivot_wider(names_from = gene2, values_from = cor) %>% 
  tibble::column_to_rownames("gene1") 

liver_wide <- data.frame(apply(liver_nw_wide, 2, as.numeric))
rownames(liver_wide) <- rownames(liver_nw_wide)

blood_nw_wide <- blood_nw %>% 
  pivot_wider(names_from = gene2, values_from = cor) %>% 
  tibble::column_to_rownames("gene1") 

blood_wide <- data.frame(apply(blood_nw_wide, 2, as.numeric))
rownames(blood_wide) <- rownames(blood_nw_wide)


net.level.liver <- networklevel(liver_wide, index=c("connectance", 
                                                  "cluster coefficient"))

net.level.blood <- networklevel(blood_wide, index=c("connectance", 
                                                  "cluster coefficient"))

grouplevel(blood_wide, index = "number of species")
grouplevel(blood_wide, index = "mean number of links")
grouplevel(blood_wide, index = "cluster coefficient")

grouplevel(liver_wide, index = "number of species")
grouplevel(liver_wide, index = "mean number of links")
grouplevel(liver_wide, index = "cluster coefficient")

