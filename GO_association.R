# script to find associations between genes
# From the parasite side:
# get all the enriched GO terms, eg. pathogenesis
# get all parasite genes in there, find out all the host interactors, do GO term analysis of host GO terms
# gives us what host GO terms pathogenesis is associated with
# repeat from the host side

library(topGO)
library(dplyr)
library(Rgraphviz)
library(biomaRt)
library(org.Mm.eg.db)
library(ggplot2)


library(dplyr)
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


options(echo=TRUE)
args <- commandArgs(TRUE)
dataset <- args[1] # enter overall or core

overall <- readRDS("blood_all_bipartite.rds")

if(dataset=="core")
{
  net <- readRDS("blood_core_network.rds"); colnames(net)[1] <- "host"; colnames(net)[2] <- "para"
  host_universe_for_core <- unique(as.character(overall[,1]))
  parasite_universe_for_core <- unique(as.character(overall[,2]))
}else
  net = overall; colnames(net)[1] <- "host"; colnames(net)[2] <- "para"


hostGOenr <- function(host_genes, GO)
{
	host_orthogroups <- read.delim("host_orthogroups.txt", stringsAsFactors=FALSE)
	host_genes <- data.frame(Orthogroup = host_genes)
  host_in <- inner_join(host_genes, host_orthogroups)
  host_in <- host_in[,c(1,3)]
  host_in <- unique(as.character(host_in[,2]))
  host_bg <- host_orthogroups[,3]
  h_geneList <- as.integer(host_bg %in% host_in)
  names(h_geneList) <- host_bg

  if(dataset == "core")
  {
  # universe for core blood 
  host_uni <- host_universe_for_core
  host_bg <- host_orthogroups$Orthogroup
  h_bg <- intersect(host_uni, host_bg)
  h_bg <- host_orthogroups[host_orthogroups$Orthogroup%in%h_bg,"mouse"]
  host_in <- host_genes[,1]
  host_in <- host_orthogroups[host_orthogroups$Orthogroup%in%host_in,"mouse"]
  
  h_geneList <- as.factor(as.integer(h_bg %in% host_in))
  names(h_geneList) <- h_bg
  } 
  
  topDiffGenes <- function(allScore) 
  {
    return(allScore == 1)
  }
  x <- topDiffGenes(h_geneList)
  
  hGOdata <- new("topGOdata",
                 ontology = "BP",
                 allGenes = h_geneList,
                 nodeSize = 10,
                 annotationFun = annFUN.org,
                 geneSelectionFun = topDiffGenes,
                 mapping = "org.Mm.eg",
                 ID = "ensembl")
  
  resultKS=runTest(hGOdata, algorithm='weight01', statistic='KS')
  allGO=usedGO(hGOdata)
  all_res=GenTable(hGOdata, KS=resultKS, orderBy="KS", topNodes=length(allGO))
  
  GenesForGOterm <- c()
  myterms = all_res$GO.ID
  mygenes <- genesInTerm(hGOdata, myterms)
  for (i in 1:length(myterms))
  {
    myterm <- mygenes[myterms[i]][[1]]
    mygenesforterm <- myterm[which(myterm %in% host_in == TRUE)]
    mygenesforterm <- paste(mygenesforterm, collapse=',')
    GenesForGOterm[i] <- mygenesforterm
  }
  
  all_res$GenesForGOterm <- GenesForGOterm
  return(all_res)
}

paraGOenr <- function(para_genes, GO)
{
	para_in <- as.matrix(para_genes)
	colnames(para_in)[1] <- "Orthogroup"
	pOG <- read.delim("parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
	geneID2GO <- readMappings(file = "p_OG_GOterms.txt") # 3659 genes
  all_res <- data.frame(GO.ID = "GO:000000", Term = "xyz", Annotated = 0, Significant = 0, Expected = 0, Fisher = 1, GenesForGOterm = "")
	# universe and genelist for overall:
  #para_in <- unique(as.character(para_genes))
  para_bg <- names(geneID2GO)
  geneList = factor(as.integer(para_bg %in% para_in))

  if(nlevels(factor(as.integer(para_bg %in% para_in))) == 2)
  {
  names(geneList)= para_bg

  if(dataset == "core")
  # universe for core network: 
  {
    para_uni <- parasite_universe_for_core
    para_in <- unique(as.character(para_genes))
    bg <- intersect(para_uni, names(geneID2GO))
    geneList = factor(as.integer(bg %in% para_in))
    names(geneList)= bg
  }
  # 

  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO,
                nodeSize = 1)
  resultKS=runTest(GOdata, algorithm='weight01', statistic='KS') 
  allGO=usedGO(GOdata)
  all_res=GenTable(GOdata, KS=resultKS, orderBy="KS", topNodes=length(allGO))
  
  GenesForGOterm <- c()
  myterms = all_res$GO.ID
  mygenes <- genesInTerm(GOdata, myterms)
  for (i in 1:length(myterms))
  {
    myterm <- mygenes[myterms[i]][[1]]
    mygenesforterm <- myterm[which(myterm %in% para_in == TRUE)]
    mygenesforterm <- sapply(mygenesforterm, 
                             function(x) pOG[grep(pattern = x, pOG$Orthogroup),"Pberghei"])
    mygenesforterm <- paste(mygenesforterm, collapse=',')
    GenesForGOterm[i] <- mygenesforterm
  }

  all_res$GenesForGOterm <- GenesForGOterm
  }
  return(all_res)
  
}

# get host and parasite prthogroups
pOG <- read.delim("parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
hOG <- read.delim("host_orthogroups.txt", stringsAsFactors=FALSE)

# get host and parasite BP GO terms
p_GO <- read.table(paste0("p_OG_topGO_BP_blood_", dataset, "_para_result.txt", collapse = ""), header = T, sep = '\t', stringsAsFactors=FALSE) %>%
        dplyr::filter(KS <= 0.05)
h_GO <- read.delim(paste0("Mmus_topGO_BP_blood_", dataset, "_host_result.txt", collapse = ''), header = T, sep = '\t', stringsAsFactors=FALSE) %>%
        dplyr::filter(KS <= 0.05)



GO_asso_from_parasite <- list()
for(i in 1:2)# nrow(p_GO)
{
	genes_in_GO_term <- as.character(p_GO[i,"GenesForGOterm"])
	if(genes_in_GO_term != "")
	{
	  genes_in_GO_term <- strsplit(genes_in_GO_term, split = ",")[[1]]

	# Find their interactors
  	hg_vector <- c()
  	for(j in 1:length(genes_in_GO_term))
  	{
  		pg <- genes_in_GO_term[j]
  		p_ortho <- pOG[grep(pattern = pg, pOG$Pberghei), "Orthogroup"]
  
  		# find its interactors
  		hg_vector <- c(hg_vector, as.character(net[grep(pattern = p_ortho, net$para), "host"]))
  	}
	

	# Go analysis of host genes
	allres <- hostGOenr(hg_vector, GO = p_GO[i,1])%>%
		dplyr::filter(KS <= 0.05)

	GO_asso_from_parasite[[i]] <- allres
	names(GO_asso_from_parasite)[i] <- paste0(p_GO[i,"GO.ID"], "_", p_GO[i,"Term"])
	}
}
saveRDS(GO_asso_from_parasite, file = "GO_asso_from_parasite.rds")


GO_asso_from_host <- list()
for(i in 1:2)# nrow(h_GO)
{
  print(i)
	hgenes_in_GO_term <- as.character(h_GO[i,"GenesForGOterm"])
	if(hgenes_in_GO_term != "")
	{ 
	hgenes_in_GO_term <- strsplit(hgenes_in_GO_term, split = ",")[[1]]

  # Find their interactors
	pg_vector <- c()
	for(j in 1:length(hgenes_in_GO_term))
	{
		hg <- as.character(hgenes_in_GO_term[j])
	 
  		h_ortho <- hOG[grep(pattern = hg, hOG$mouse), "Orthogroup"]
  
  		# find its interactors
  		pg_vector <- c(pg_vector, as.character(net[grep(pattern = h_ortho, net$host), "para"]))

  	}

  	# Go analysis of host genes
  	allres <- paraGOenr(para_genes = as.character(pg_vector), GO = h_GO[i,1])
  	allres <- allres[which(allres$KS <= 0.05),]
  		
  	GO_asso_from_host[[i]] <- allres
  	names(GO_asso_from_host)[i] <- paste0(h_GO[i,"GO.ID"], "_", h_GO[i,"Term"])
	}
}

saveRDS(GO_asso_from_host, file = "GO_asso_from_host.rds")

# make GO term edges to visualise as networks

GO_asso_from_parasite_edges <- data.frame()
a = 1
for(i in 1:length(GO_asso_from_parasite))
{
	rows <- nrow(GO_asso_from_parasite[[i]])
	if(length(rows) != 0)
	{GO_asso_from_parasite_edges[a:(a+rows-1),1] <- rep(strsplit(names(GO_asso_from_parasite)[i], split = "_")[[1]][2], rows)
	GO_asso_from_parasite_edges[a:(a+rows-1),2] <- GO_asso_from_parasite[[i]]$Term

	a = a+rows}
}

g = 1
GO_asso_list = list()
for(k in 1:length(GO_asso_from_host))
{
  if(!is.null(GO_asso_from_host[[k]]))
    {
      GO_asso_list[[g]] <- GO_asso_from_host[[k]]
      names(GO_asso_list)[[g]] <- names(GO_asso_from_host)[[k]]
      g = g+1
    }
}

GO_asso_from_host_edges <- data.frame()
a = 1
for(i in 1:length(GO_asso_list))
{
  if(nrow(GO_asso_list[[i]]) > 0)
  {
   rows <- nrow(GO_asso_list[[i]])
	if(!is.null(rows))
	{
	GO_asso_from_host_edges[a:(a+rows-1),1] <- rep(strsplit(names(GO_asso_list)[i], split = "_")[[1]][2], rows)
	GO_asso_from_host_edges[a:(a+rows-1),2] <- GO_asso_list[[i]]$Term

	a = a+rows
	}
}
}

colnames(GO_asso_from_host_edges) <- c("Host", "Parasite")
colnames(GO_asso_from_parasite_edges) <- c("Parasite", "Host")

# Put host and parasite GO associations together; find distinct rows
GO_asso <- full_join(GO_asso_from_host_edges, GO_asso_from_parasite_edges) %>%
  dplyr::distinct(.)

saveRDS(GO_asso, file = paste0("GO_asso_", dataset, ".rds", collapse = ''))
write.csv2(GO_asso, paste0("GO_asso_", dataset, ".csv", collapse = ''), row.names = F, quote = F)


# attach node properties to the GO asso network
GO_asso_nodes_host <- data.frame(Term = unique(c(unique(as.character(GO_asso_from_host_edges[,1])), unique(as.character(GO_asso_from_parasite_edges[,2])))))
GO_asso_nodes_parasite <- data.frame(Term = unique(c(unique(as.character(GO_asso_from_host_edges[,2])), unique(as.character(GO_asso_from_parasite_edges[,1])))))

GO_asso_nodes <- data.frame(Term = c(GO_asso_nodes_host[,1], GO_asso_nodes_parasite[,1]), 
  Organism = c(rep("Host", length(GO_asso_nodes_host[,1])), rep("Parasite", length(GO_asso_nodes_parasite[,1]))))

write.csv2(GO_asso_nodes, paste0("GO_asso_nodes_", dataset, ".csv", collapse = ''), row.names = F, quote = F)
