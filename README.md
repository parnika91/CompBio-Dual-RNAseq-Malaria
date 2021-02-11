# CompBio-Dual-RNAseq-Malaria
Computational analysis of host parasite interactions in malaria from existing data


Malaria is one of the most researched disease caused by a eukaryotic organism. RNA-sequencing has become the go-to method to study gene expression and is heavily implemented in malaria studies to understand disease mechanisms during an infection. While publishing scientific work that involves sequencing of samples, scientists are required to submit the sample sequences onto NCBI in the Sequence Read Archive (SRA) or on one of the mirror databases - ENA and DDBJ. For understanding any disease, and in our case, malaria, these databases are an excellent resource for a vast amount of sequence information and metadata.

We used sequence database, SRA, to obtain relevant malaria RNA-seq studies to perform a meta-analysis and infer host-parasite interactions. 

Steps 1 - 6 constitute the upstream processing of data and as results, we get the bipartite or host-parasite interactions from each study that is analysed. In this repo, I have uploaded the blood overall and core networks in .rds and .csv formats where possible with orthogroup names and gene IDs. If you want to run the analysis downstream of retrieving these networks, please feel free to download the networks and run any script atfer Step 6. The scripts were written to run with the networks in orthogroup names.

Step 1:  Find single copy orthologous genes for hosts and for parasites. Here I used OrthoFinder1.0.6.  I put all protein sequences of a species in one file and all the proteome files of the species of interest in one folder to run OrthoFinder. The sequences for _Plasmodium sp._ were obtained from PlasmoDB and for the host from Ensembl. Host orthogroups could also be obtained from ensembl biomart. Single copy orthologous genes are given an orthogroup ID, such as "h\_OG000456" for host and "p\_OG003456" for parasite.  The way I had arranged these orthogoups were: column 1 with orthogroup ID followed by the gene names from individual species in the same row. To run OrthoFinder, I ran: ./orthofinder -f [path to folder with all parasite protein sequences] -t  6

The host and parasite orthogroups and their corresponsing genes in the different species are in files host\_orthogroups.txt and parasite_orthogroups.txt.

Step 2:  To get host-parasite interactions from eligible studies, first list out which studies are eligible, the host organism and the parasite organism, separated by comma into positive_experiments.txt. This file is required later in several steps.
Second, for each study, retreive the metadata from its ENA (European Nucleotide Archive) website and store it as <studyID>_accession.txt. The 7th column in this file has the ftp links to download the runs in fastq format. In another file, runs\_<studyID>.txt, store the runIDs and studyIDs of a study in two columns separated by comma. Later, I will write scripts to make these two files automatically. 
In the command line, run sh TranscriptomeAnalysisENA.sh <txt file with all studies for analysis in a column>
This will first make a separate file only with the ftp links called download.txt inside studyID/ and then it will call script DownloadENA.sh to download all the runs in a study. Next, it invokes doallENA.sh which in turn calls other scripts to index genome files, map the reads in a run and then quantify gene expression.

Step 3: The studyIDs can be submitted as user input in the command line to post\_process\_countfile.R to bring all the runs of a study together, extract their protein coding genes because that's what I want to characterise here and to keep those that belong to an orthogroup created in Step 1. If thresholds are to be applied on the runs being included for analysis, I distinguished them as <studyID>.ortho.data.<threshold>.RData, eg., DRP000987.ortho.data.str.RData. The script for applying thresholds are not included in this repo at the moment and I will add an example later. For the explanation on thresholds, please see "Selection of runs for analysis" https://www.biorxiv.org/content/10.1101/576116v6.

Step 4: After post-processing, Correlation_dataset.R can be run with the name of the dataset, eg DRP000987.ortho.data.str, as an input in the command line. This script performs 10000 permutation tests to get uncorrected p-values as defined in "Identification of co-expressed genes via correlation techniques" section of https://www.biorxiv.org/content/10.1101/576116v6. The numerical scores obtained at the end of this script can be converted to "p-values" by doing (score + 1)/10000. For my needs, I have kept the original scores, meaning that a score of 0 refers to the most robustly inferred association/interaction.

Step 5: From these results, by running Correlation_analyse.R, the bipartite part of the network can be extracted. I have chosen only those associations that had a score of 0. This takes the form of an edge table, with the first column as host orthogroups, second column for parasite orthogroups and the third column for the correlation strength. With this script, I also made the heatmaps showing the common (or overlapping/intersecting) edges between studies, the significance of these overlaps and the jaccard indices for each intersection set between a pair of studies (fig 2) and plots showing the the reduction in the number of edges with score 0 with increase in the number of permutations (fig 4C) and the distribution of permutation scores across the range of correlation coefficients (fig 4D). 
Right now, the script uses the blood overall permutation test results as an example. Once I add all the bipartite edge files from all the studies to this repository, I will make the downstream part of the script more reproducible.

Step 6: Next, using the overall network and the individual dataset networks in Core_network.R, I consructed the blood core network as described in section "A 'core' network of evolutionarily conserved interactions" in the paper. This also includes the construction of fig 4A.

The core network and the overall networks are uploaded to the repo. The following steps can be performed with these two networks. 

Step 7: Run topGO.R to perform GO analysis of the host and parasite genes in the overall or core networks. From the command line where the networks were downloaded, run `Rscript topGO.R <network>` (network = `overall` for an analysis of the overall network or `core` to analyse the core network). 

Step 8: Run GO_association.R to find GO biological processes associated or interacting between the host and the parasite. From the command line where the networks were downloaded, run `Rscript  GO_association.R <network>` (network = `overall` for an analysis of the overall network or `core` to analyse the core network). This can be run only after topGO.R has been for the desired network.

Step 9: Run Modularity_Immunemarker.R to find modules in the overall or core network, the immune genes present and their enrichment in the network based on Vallejo et al., 2018. From the command line where the networks were downloaded, run `Rscript Modularity_Immunemarker.R <network>` (network = `overall` for an analysis of the overall network or `core` to analyse the core network). Link to the supplementary file in Vallejo et. al, 2018 that was used to find the immune genes in the network: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6203889/bin/mmc5.xlsx

Step 10: Run `Rscript beta_models.R <network>` (network = `overall` for an analysis of the overall network or `core` to analyse the core network) to make the beta regression models with network properties. For both overall and core network, this script, right now, only makes the models for the parasite-parasite edges. 

The following packages were required to run all the scripts in this repository. (The packages required for each script is written on top of each script.)
```
library(tidyverse)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)
library(stringr)
library(stargazer)
library(topGO)
library(biomaRt)
library(org.Mm.eg.db)
library(WGCNA)
library(reshape2)
library(UpSetR)
library(grid)
library(limma)
library(EnsDb.Hsapiens.v79)
library(SRAdb)
library(GEOquery)
library(pkgconfig)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(foreach)
library(doParallel)
```
