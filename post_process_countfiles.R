## This script is to make the gene count table ready for coexpression analysis
# Take as input a .txt file like that of TranscriptomeAnalysisENA.sh or the file current_studies_for_analysis.txt made in the same script

### Concatenate all runs in a study
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(S4Vectors)
library(dplyr)

options(echo=TRUE)
args <- commandArgs(TRUE)
studyIDs <- args[1] # a file name
studyIDs <- read.table(paste0(studyIDs, ".txt", collapse = ''), header = T, sep = '\t')
studyIDs <- studyIDs[,1]
positive_experiments <- read.delim("positive_experiments.txt", sep = ",", header = F)

################# Step 1: Bring all runs together ###################

ConcatRunsToStudy <- data.frame()

for(j in 1:length(studyIDs))
{
  print(studyIDs[j])
  runIDs <- read.table(paste0(studyIDs[j], "/runs_",studyIDs[j],".txt", collapse = ''), header = F, sep = ',')
  number_of_runs <- nrow(runIDs)

  if(number_of_runs == 1)
  {
    firstRun <- runIDs[1,1]
    FirstCountfile <- read.table(paste0(studyIDs[j], "/countWithGFF3_",firstRun,".txt",collapse=''), header = T, sep = '\t')
    write.table(FirstCountfile, paste0(studyIDs[j],"/ConcatRunsToStudy_", studyIDs[j],".txt"), sep = '\t', row.names=F)
  }

  if(number_of_runs > 1)
  {
    firstRun <- runIDs[1,1]
    if(file.exists(paste0(studyIDs[j], "/countWithGFF3_",firstRun,".txt",collapse='')))
    {
      FirstCountfile <- read.table(paste0(studyIDs[j], "/countWithGFF3_",firstRun,".txt",collapse=''), header = T, sep = '\t')
      ConcatRunsToStudy <- FirstCountfile
      colnames(ConcatRunsToStudy)[6] <- paste0(as.character(firstRun), "_",as.character(studyIDs[j]),collapse='')
    }
    a = 2
    
    for(i in 2:number_of_runs)
    {
      runID <- runIDs[i,1]
      # get runID
      if(file.exists(paste0( studyIDs[j], "/countWithGFF3_",runID,".txt",collapse='')))
        {
          countfile <- read.table(paste0(studyIDs[j], "/countWithGFF3_",runID,".txt",collapse=''), header = T, sep = '\t')
  
          ConcatRunsToStudy[1:nrow(FirstCountfile),(a+5)] <- countfile[,6]
          #ConcatRunsToStudy <- merge(ConcatRunsToStudy, countfile, by = c("seqnames", "start", "end", "width", "strand"))
          colnames(ConcatRunsToStudy)[(a+5)] <- paste0(as.character(runID), "_",as.character(studyIDs[j]),collapse='')
          a = a+1
        }  
    }
    write.table(ConcatRunsToStudy, paste0(studyIDs[j], "/ConcatRunsToStudy_", studyIDs[j],".txt", collapse=''), sep = '\t', row.names=F)
  }
}


################################# Step 2: Get gene names for all reads #################

for( i in 1:length(studyIDs))
{
  print(studyIDs[i])
  study <- read.csv2(paste0(studyIDs[i], "/ConcatRunsToStudy_", studyIDs[i],".txt", collapse=''), sep = '\t', header=T)
  
  #get host and parasite
  
  host <- as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),2])
  para <- as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),3])
  
  genes <- import(paste0("Genomes/annotation/",host,para,".gtf", collapse=''), format = "gtf")
  genes <- genes[genes$type%in%"exon"]
  genes.df <- as.data.frame(genes)
  genes.df.gene_name <- genes.df[,c("seqnames", "start", "end", "width", "strand", "gene_id")]
  
  mergeStudy.genes.df.gene_name <- merge(study, genes.df.gene_name, by = c("seqnames", "start", "end", "width", "strand"))
  
  mergeStudy.genes.df.gene_name <- mergeStudy.genes.df.gene_name[,6:ncol(mergeStudy.genes.df.gene_name)]
  mergeStudy.genes.df.gene_name.combineGenes <- data.frame()
  mergeStudy.genes.df.gene_name.combineGenes <- aggregate(mergeStudy.genes.df.gene_name[,1] ~ gene_id, data = mergeStudy.genes.df.gene_name, sum)
  colnames(mergeStudy.genes.df.gene_name.combineGenes)[2] <- colnames(mergeStudy.genes.df.gene_name)[1]
  
  if(ncol(study) > 6)
  {
    for(k in 2:(ncol(mergeStudy.genes.df.gene_name)-1))
    {
      agg <- aggregate(mergeStudy.genes.df.gene_name[,k] ~ gene_id, data = mergeStudy.genes.df.gene_name, sum)
      mergeStudy.genes.df.gene_name.combineGenes <- merge(mergeStudy.genes.df.gene_name.combineGenes, agg, by = c("gene_id"))
      colnames(mergeStudy.genes.df.gene_name.combineGenes)[k+1] <- colnames(mergeStudy.genes.df.gene_name)[k]
    }
  }
  
  t.study <- t(mergeStudy.genes.df.gene_name.combineGenes)
  colnames(t.study) <- mergeStudy.genes.df.gene_name.combineGenes$gene_id
  t.study <- t.study[-1,]
  class(t.study) <- "numeric"
  write.table(t.study, paste0(studyIDs[i], "/", studyIDs[i],".txt", collapse=''), sep = '\t', row.names=T)
}


#################################### Step 3: Make a new countfile with only the coding genes ########################################

# Keep the one with pseudogenes in case it is required later on

################ host-parasite pairs ###############

for(i in 1:length(studyIDs))
{
  print(i)
  # get study.txt including all runs
  study <- as.data.frame(t(read.delim(paste0("Output/", studyIDs[i], "/", studyIDs[i], ".txt", collapse = ''), sep = '\t', header = T)))
  hp <- paste(as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),2]), 
    as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),3]), sep="")

  # Possible values of hp are:
  # humanPfalciparum
  # humanPvivax
  # humanPberghei
  # mousePchabaudi
  # mousePberghei
  # mousePyoelii
  # monkeyPcynomolgi
  # monkeyPcoatneyi
  
  coding = as.data.frame(import(paste0("Genomes/annotation/", hp, ".gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id)
  
  # filter study to keep only protein-coding genes
  study_coding_genes <- study %>%
    tibble::rownames_to_column('gene') %>%
    filter(rownames(study)%in%coding$gene_id) %>%
    tibble::column_to_rownames('gene')
   
  # write the table out
  write.table(study_coding_genes, paste0(studyIDs[i],"/", studyIDs[i], "_coding_genes.txt", collapse = ''), sep ='\t', row.names = T)
}

############################## Step 4: Get orthologous groups for each study ########################

# OrthoFinder was used to obtain the orthogroups for the hosts and parasites
# to run OrthoFinder: ./orthofinder -f [path to folder with all parasite protein sequences] -t  6
# parasite_orthogroups.txt contains only single-copy orthologs extracted from the output of OrthoFinder
# Its first column contains the orthogroup ID, such as p_OG000001, followed by the corresponding genes from each parasite species

parasite_orthogroups <- read.delim("OrthoFinder/parasite_orthogroups.txt", stringsAsFactors=FALSE) 
host_orthogroups <- read.delim("OrthoFinder/host_orthogroups.txt", stringsAsFactors=FALSE) 

for(i in 1:length(studyIDs))
{ print(i)
  #if the study_coding_genes exists, merge with orthogroups (join functions)
  filepath = paste0(studyIDs[i],"/",studyIDs[i],"_coding_genes.txt", collapse = "")
  if(file.exists(filepath))
  {
    # find out what host and parasite the study is
    host <- as.character(positive_experiments[grep(pattern = studyIDs[i], positive_experiments[,1]),2])
    para <- as.character(positive_experiments[grep(pattern = studyIDs[i], positive_experiments[,1]),3])

    # take the host and para orthogroups and make a df -> orthogroup | gene name 

    h.df <- data.frame(Orthogroup = host_orthogroups[,1], Org = host_orthogroups[,grep(pattern = host, colnames(host_orthogroups))])
    p.df <- data.frame(Orthogroup = parasite_orthogroups[,1], Org = parasite_orthogroups[,grep(pattern = para, colnames(parasite_orthogroups))])

    hp.df <- rbind(h.df, p.df)

    # read table
    file = read.delim(filepath, header = T) %>% tibble::rownames_to_column("Gene")

    ortho.table = merge(file, hp.df, by.x = "Gene", by.y = "Org")
    ortho.table <- data.frame(Gene = ortho.table$Gene, Orthogroup = ortho.table$Orthogroup, ortho.table[,2:(ncol(ortho.table)-1)])

    write.table(ortho.table, paste0(studyIDs[i],"/",studyIDs[i],"_orthogroups.txt", collapse = ""), sep = '\t', row.names = F)
  }

}