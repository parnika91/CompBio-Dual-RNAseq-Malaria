# To perform GO analysis with topGO package from Bioconductor
# to do this, we first need to do cexpression analysis, and get the genes in the network
# Here I have found enriched GO terms for the bipartite part of the network

# In addition, we also need the .gaf (gene annotation file) file for each parasite here (ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pberghei/)
# To maximaise the number of GO terms annotated for Plasmodium sp., I also downloaded the ensembl biomart Go annotations for each species
# These two resources were combined to obtain GO terms for a single species
# These terms were then merged with the orthogroup IDs to get GO annotations for the orthogroups
# I found GO annotations for P. falciparum, P. vivax and P. berghei. Their terms were combined to encompass as many species as possible.

# For hosts, this was not necessary, as I used the dataset from the org.Mm.eg.db package in Bioconductor, which contains genome wide annotation for mouse
# For parasites this was not prefered as the only dataset available was for P. falciparum. 
# First, P. falciparum is phylogenetically quite distant to the other Plasmodium species here, and
# second, the dataset is much smaller than the dataset created from using orthogroups of 3 species.

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
  
#BiocManager::install("topGO")
library(topGO)
library(dplyr)
library(Rgraphviz)
library(biomaRt)
library(org.Mm.eg.db)
library(ggplot2)

# function to assign an RData object to a name
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## Parasite gaf and biomart files are read in

# Pb_Gdb <- read.delim("~/topGO/Pberghei.gaf", header=FALSE, 
#                      stringsAsFactors=FALSE)
# Pb_bm <- read.delim("~/topGO/Pb_mart_export.txt", header=T,
#                    stringsAsFactors=FALSE)
# 
# Pf_Gdb <- read.delim("~/topGO/Pfalciparum.gaf", header=FALSE, 
#                      stringsAsFactors=FALSE)
# Pf_bm <- read.delim("~/topGO/Pf_mart_export.txt", header=T,
#                     stringsAsFactors=FALSE)
# 
# Pv_Gdb <- read.delim("~/topGO/PvivaxP01.gaf", header=FALSE, 
#                      stringsAsFactors=FALSE)
# Pv_bm <- read.delim("~/topGO/Pv_mart_export.txt", header=T,
#                     stringsAsFactors=FALSE)


# # treating Pf tables
# # remove rows with empty GO accession from both
# Pf_bm <- Pf_bm[which(Pf_bm$GO.term.accession != ""),]
# Pf_Gdb <- Pf_Gdb[which(Pf_Gdb$V5 != ""),]
# 
# # remove .1 from Pf_Gdb
# Pf_Gdb$V2 <- substr(Pf_Gdb$V2, 1, 13)
# 
# # collect rows wih gene IDs and GO term
# Pf_bm_red <- Pf_bm[,c("Gene.stable.ID", "GO.term.accession")]
# Pf_Gdb_red <- Pf_Gdb[,c("V2", "V5")]
# colnames(Pf_Gdb_red) <- colnames(Pf_bm_red)
# 
# # rbind bm and Gdb
# Pf <- rbind(Pf_bm_red, Pf_Gdb_red)
# Pf <- aggregate(Pf$GO.term.accession,Pf['Gene.stable.ID'],paste,collapse=',')
#  
# # keep unique GOterms
# Pf$GO <- sapply(Pf$x, function(y) paste(unique(strsplit(y, split = ",")[[1]]), collapse = ","))
# Pf <- Pf[, -c(2)]


# # treating Pb tables
# # remove rows with empty GO accession from both
# Pb_bm <- Pb_bm[which(Pb_bm$GO.term.accession != ""),]
# Pb_Gdb <- Pb_Gdb[which(Pb_Gdb$V5 != ""),]
# 
# # remove .1 from Pb_Gdb
# Pb_Gdb$V2 <- substr(Pb_Gdb$V2, 1, 13)
# 
# # collect rows wih gene IDs and GO term
# Pb_bm_red <- Pb_bm[,c("Gene.stable.ID", "GO.term.accession")]
# Pb_Gdb_red <- Pb_Gdb[,c("V2", "V5")]
# colnames(Pb_Gdb_red) <- colnames(Pb_bm_red)
# 
# # rbind bm and Gdb
# Pb <- rbind(Pb_bm_red, Pb_Gdb_red)
# Pb <- aggregate(Pb$GO.term.accession,Pb['Gene.stable.ID'],paste,collapse=',')
# 
# # keep unique GOterms
# Pb$GO <- sapply(Pb$x, function(y) paste(unique(strsplit(y, split = ",")[[1]]), collapse = ","))
# Pb <- Pb[, -c(2)]


# # treating Pv tables
# # remove rows with empty GO accession from both
# Pv_bm <- Pv_bm[which(Pv_bm$GO.term.accession != ""),]
# Pv_Gdb <- Pv_Gdb[which(Pv_Gdb$V5 != ""),]
# 
# # remove .1 from Pb_Gdb
# Pv_Gdb$V2 <- substr(Pv_Gdb$V2, 1, 13)
# 
# # collect rows wih gene IDs and GO term
# Pv_bm_red <- Pv_bm[,c("Gene.stable.ID", "GO.term.accession")]
# Pv_Gdb_red <- Pv_Gdb[,c("V2", "V5")]
# colnames(Pv_Gdb_red) <- colnames(Pv_bm_red)
# 
# # rbind bm and Gdb
# Pv <- rbind(Pv_bm_red, Pv_Gdb_red)
# Pv <- aggregate(Pv$GO.term.accession,Pv['Gene.stable.ID'],paste,collapse=',')
# 
# # keep unique GOterms
# Pv$GO <- sapply(Pv$x, function(y) paste(unique(strsplit(y, split = ",")[[1]]), collapse = ","))
# Pv <- Pv[, -c(2)]


# save(Pf, file = "Pf_annot.RData")
# save(Pb, file = "Pb_annot.RData")
# save(Pv, file = "Pv_annot.RData")


# # merge with orthogroups
# Pb_OG <- merge(Pb, pOG[,c("Orthogroup", "Pb_g")], by.x = "Gene.stable.ID", by.y = "Pb_g")
# Pf_OG <- merge(Pf, pOG[,c("Orthogroup", "Pf_g")], by.x = "Gene.stable.ID", by.y = "Pf_g")
# Pv_OG <- merge(Pv, pOG[,c("Orthogroup", "Pv_g")], by.x = "Gene.stable.ID", by.y = "Pv_g")
# 
# Pb_Pf <- full_join(Pb_OG, Pf_OG, by = "Orthogroup")
# Pb_Pf_Pv <- full_join(Pb_Pf, Pv_OG, by = "Orthogroup")
# 
# # keeping only orthogroups IDs
# GO <- data.frame(paste(as.character(Pb_Pf_Pv$GO.x), as.character(Pb_Pf_Pv$GO.y), as.character(Pb_Pf_Pv$GO), sep = ","))
# GO_ <- data.frame()
# for(i in 1:nrow(GO))
# {
#   t <- unique(strsplit(as.character(GO[i,]), split = ",")[[1]])
#   t <- t[!t %in% "NA"]
#   GO_[i,1] <- paste(t, collapse = ",")
# }
# 
# Pb_Pf_Pv_OG <- data.frame(Orthogroup = Pb_Pf_Pv$Orthogroup, GOterm = GO_)
# colnames(Pb_Pf_Pv_OG)[2] <- "GO"
# save(Pb_Pf_Pv_OG, file = "p_OG_GOterms.RData")
# colnames(Pb_Pf_Pv_OG) <- NULL


###################### This is the file that can be used here as universal geneset ####################
# write.table(Pb_Pf_Pv_OG, "p_OG_GOterms.txt", sep = '\t', row.names = F, quote = F)

### Accept as input the studyID and the network that we want to do the GO analysis for
# for example ,study DRP000987 and dataset DRP000987_str_bipartite.RData, which is the bipartite network in the "str" threshold of DRP000987
# such an Rdata file would have the first column as the host genes and the second column as the parasite genes that each host gene is associated (interacting) with
# Note please that these are still represented as orthogroups IDs and not really the gene names of the individual species
options(echo=TRUE)
args <- commandArgs(TRUE)
studyID <- args[1]
dataset <- args[2] #, eg DRP000987_str_bipartite.RData

net <- loadRData(paste0(studyID, "/", dataset, ".RData", collapse = ''))
para_genes <- unique(as.character(net[,2]))
host_genes <- unique(as.character(net[,1]))
write.table(para_genes, paste0(studyID, "/", studyID, "_para_genes.txt", collapse = ''), quote = F, row.names = F)
write.table(host_genes, paste0(studyID, "/", studyID, "_host_genes.txt", collapse = ''), quote = F, row.names = F)

pOG <- read.delim("~/Documents/Data/parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)

bip_studies <- studyID

# for example, if we want to do this aalysis for many studies at once, then the above step can be done first by commenting out the next step
# and then in bip_studies, include all the studyIDs as a vector as shown in the comments below this block and run the for loop below.
# The for loop below loops over each study and BP, CC, MF for parasite and for host
# I will make this process more elegant in the future

# bip_studies <- c("SRP250329", "ERP105548", 
#                  "SRP110282", "SRP096160",
#                  "liver.int.overall")
# bip_studies <- c("DRP000987","ERP106451",
#                  "ERP110375", "ERP004598",
#                  "SRP118827", "SRP116793"
#                  )


geneont <- c("BP", "CC", "MF")

for(m in 1:length(bip_studies))
{
  print(bip_studies[m])
  for(n in 1:length(geneont))
  {
    print(geneont[n])
    study <- bip_studies[m]
    GeneOnt <- geneont[n]
    
     
    geneID2GO <- readMappings(file = "topGO/p_OG_GOterms.txt") # 3659 genes
    
    para_genes <- read.delim(paste0("~/Documents/Data/", study ,"_para_genes.txt",
                                  collapse = ""), stringsAsFactors=FALSE, header = T)

    # parasite genes of interest
    para_in <- unique(as.character(para_genes[,1]))
    # universe of parasite genes
    para_bg <- names(geneID2GO)
    # to know which genes are interesting within the universe, we do %in% with background genes
    geneList = factor(as.integer(para_bg %in% para_in))
    # make it a named vectore, that's what topGO requires
    names(geneList)= para_bg

    GOdata <- new("topGOdata",
                  ontology = GeneOnt,
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
                  nodeSize = 1)
    # Expected: Under random chance, number of genes that would be expected 
    # to be significantly DE and annotated with that term
    # The column Expected represents the expected number of interesting genes mapped to the 
    # GO term if the interesting genes were randomly distributed over all GO terms.
    
    resultKS=runTest(GOdata, algorithm='weight01', statistic='Fisher') 
    allGO=usedGO(GOdata)
    all_res=GenTable(GOdata, Fisher=resultKS, orderBy="Fisher", topNodes=length(allGO))
    #par(cex = 0.4)
    ## Plotting results
    #showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo ='all')
    ##red: significant, yellow = not sig
    #printGraph(GOdata, resultKS, firstSigNodes = 5, fn.prefix = paste0("tGO_", study,"_parasite_", GeneOnt, collapse = ''), useInfo = "all", pdfSW = TRUE)
    
    # Get genes in a particular GO term:
    GenesForGOterm <- c()
    myterms = all_res$GO.ID
    mygenes <- genesInTerm(GOdata, myterms)
    for (i in 1:length(myterms))
    {
      myterm <- mygenes[myterms[i]][[1]]
      mygenesforterm <- myterm[which(myterm %in% para_in == TRUE)]
      mygenesforterm <- sapply(mygenesforterm, 
                               function(x) pOG[grep(pattern = x, pOG$Orthogroup),"Pb_g"])
      mygenesforterm <- paste(mygenesforterm, collapse=',')
      GenesForGOterm[i] <- mygenesforterm
    }

    all_res$GenesForGOterm <- GenesForGOterm
    write.table(all_res, paste0("p_OG_topGO_", GeneOnt,"_", study, "_para_result.txt", collapse = ''), sep = '\t', row.names = F)
    
    goEnrichment <- all_res[all_res$KS<0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","KS", "Significant")]
    goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
    goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
    goEnrichment$KS <- as.numeric(goEnrichment$KS)
    
    ggplot(goEnrichment, aes(x=Term, y=log10(Significant), fill = KS)) +
      stat_summary(geom = "bar", width = 0.5, fun = mean, position = "dodge") +
      xlab("Molecular Function") +
      ylab("Number of parasite genes in GO term") +
      ggtitle("Significant Plasmodium GOterms") +
      #scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
      theme_bw(base_size=20) +
      theme(
        #legend.position='none',
        #legend.background=element_rect(),
        plot.title=element_text(angle=0, size=20, face="bold", vjust=1),
        axis.text.x=element_text(angle=0, size=20, face="bold", hjust=1.10),
        axis.text.y=element_text(angle=0, size=20, face="bold", vjust=0.5),
        axis.title=element_text(size=20, face="bold"),
        #legend.key=element_blank(),     #removes the border
        #legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        #legend.text=element_text(size=8),  #Text size
        title=element_text(size=18)) +
      #guides(colour=guide_legend(override.aes=list(size=2.5))) +
      coord_flip()
      ggsave(paste0(study, "_", GeneOnt, "_p_OG_Enrichment.png"), height = 30, width = 45, units = "cm")
      


    # Host
    # Background genes
    ############################################## good code ################
    host_genes <- read.delim(paste0("~/Documents/Data/", study ,"_host_genes.txt",
                                  collapse = ""), header =T,
                           stringsAsFactors=FALSE)
    host_orthogroups <- read.delim("~/Downloads/host_orthogroups.txt", stringsAsFactors=FALSE)
    colnames(host_genes)[1] <- "Orthogroup"
    host_in <- inner_join(host_genes, host_orthogroups)
    host_in <- host_in[,c(1,3)]
    host_in <- unique(as.character(host_in[,2]))
    
    host_bg <- host_orthogroups[,3]
    ########################################################################

    ### The code here might work but I didn't end up using it
    
    # library(biomaRt)
    # getBM() wasn't returning results
    # ensembl <- useMart("ensembl")
    # datasets <- listDatasets(ensembl)
    # ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
    # go_ids = getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), 
    #               filters='external_gene_name', values=host_bg, mart=ensembl)
    
    # downloaded dataset from biomart webpage
    # Mmus_annot <- read.delim("~/Downloads/Mmus_annot.txt", stringsAsFactors=FALSE)
    # Mmus_annot <- Mmus_annot[-which(Mmus_annot[,3]==""),c(1,3)]
    # Mmus_annot <- aggregate(Mmus_annot$GO.term.accession, Mmus_annot['Gene.stable.ID'],paste,collapse=', ')
    # write.table(Mmus_annot, "Mmus_agg_annot.txt", row.names = F, sep = "\t", quote = F)
    
    ##################### good code ##################################
    h_geneList <- as.integer(host_bg %in% host_in)
    names(h_geneList) <- host_bg
    
    topDiffGenes <- function(allScore) 
    {
      return(allScore == 1)
    }
    x <- topDiffGenes(h_geneList)
    
    hGOdata <- new("topGOdata",
                   ontology = GeneOnt,
                   allGenes = h_geneList,
                   nodeSize = 10,
                   annotationFun = annFUN.org,
                   geneSelectionFun = topDiffGenes,
                   mapping = "org.Mm.eg",
                   ID = "ensembl")
    
    resultKS=runTest(hGOdata, algorithm='weight01', statistic='Fisher')
    allGO=usedGO(hGOdata)
    all_res=GenTable(hGOdata, Fisher=resultKS, orderBy="Fisher", topNodes=length(allGO))
    ################################################################
    
    # Plotting results
    #par(cex = 0.4)
    #showSigOfNodes(hGOdata, score(resultKS), firstSigNodes = 5, useInfo ='all')
    # # # red: significant, yellow = not sig
    #printGraph(hGOdata, resultKS, firstSigNodes = 5, fn.prefix = paste0("tGO_", study,"_host_", GeneOnt, collapse = ''), useInfo = "all", pdfSW = TRUE)
    ############################
    #
    # # Get genes in a particular GO term:
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
    write.table(all_res, paste0("Mmus_topGO_", GeneOnt,"_", study, "_host_result.txt", collapse = ''), sep = '\t', row.names = F)
    
    goEnrichment <- all_res[all_res$KS<0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","KS", "Significant")]
    goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
    goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
    goEnrichment$KS <- as.numeric(goEnrichment$KS)
    
    ggplot(goEnrichment, aes(x=Term, y=Significant, fill = KS)) +
      stat_summary(geom = "bar", width = 0.3, fun = mean, position = "dodge") +
      xlab(GeneOnt) +
      ylab("Number of host genes in GO term") +
      ggtitle(study) +
      #scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
      theme_bw(base_size=14) +
      theme(
        #legend.position='none',
        #legend.background=element_rect(),
        plot.title=element_text(angle=0, size=10, face="bold", vjust=1),
        axis.text.x=element_text(angle=0, size=10, face="bold", hjust=1.10),
        axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5),
        axis.title=element_text(size=10, face="bold"),
        #legend.key=element_blank(),     #removes the border
        #legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        #legend.text=element_text(size=8),  #Text size
        title=element_text(size=8)) +
      #guides(colour=guide_legend(override.aes=list(size=2.5))) +
      coord_flip()
    ggsave(paste0(study, "_", GeneOnt, "_host_Enrichment.png"), height = 20, width = 20, units = "cm")
    
  }
}

