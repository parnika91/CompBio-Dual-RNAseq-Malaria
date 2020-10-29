# CompBio-Dual-RNAseq-Malaria
Computational analysis of host parasite interactions in malaria from existing data


Malaria is one of the most researched disease caused by a eukaryotic organism. RNA-sequencing has become the go-to method to study gene expression and is heavily implemented in malaria studies to understand disease mechanisms during an infection. While publishing scientific work that involves sequencing of samples, scientists are required to submit the sample sequences onto NCBI in the Sequence Read Archive (SRA) or on one of the mirror databases - ENA and DDBJ. For understanding any disease, and in our case, malaria, these databases are an excellent resource for a vast amount of sequence information and metadata.

We used sequence database, SRA, to obtain relevant malaria RNA-seq studies to perform a meta-analysis and infer host-parasite interactions. 

Step 1:  Find single copy orthologous genes for hosts and for parasites. Here I used OrthoFinder1.0.6.  I put all protein sequences of a species in one file and all the proteome files of the species of interest in one folder to run OrthoFinder. The sequences for _Plasmodium sp._ were obtained from PlasmoDB and for the host from Ensembl. Host orthogroups could also be obtained from ensembl biomart. Single copy orthologous genes are given an orthogroup ID, such as "h\_OG000456" for host and "p\_OG003456" for parasite.  The way I had arranged these orthogoups were: column 1 with orthogroup ID followed by the gene names from individual species in the same row. To run OrthoFinder, I did: ./orthofinder -f [path to folder with all parasite protein sequences] -t  6

Step 2:  To get host-parasite interactions from eligible studies, first list out which studies are eligible, the host organism and the parasite organism, separated by comma into positive_experiments.txt. This file is required later in several steps.
Second, for each study, retreive the metadata from its ENA (European Nucleotide Archive) website and store it as <studyID>_accession.txt. The 7th column in this file has the ftp links to download the runs in fastq format. In another file, runs\_<studyID>.txt, store the runIDs and studyIDs of a study in two columns separated by comma. Later, I will write scripts to make these two files automatically. 
In the command line, run sh TranscriptomeAnalysisENA.sh <txt file with all studies for analysis in a column>
This will first make a separate file only with the ftp links called download.txt inside studyID/ and then it will call script DownloadENA.sh to download all the runs in a study. Next, it invokes doallENA.sh which in turn calls other scripts to index genome files, map the reads in a run and then quantify gene expression.

Step 3: The studyIDs can be submitted as user input in the command line to post\_process\_countfile.R to bring all the runs of a study together, extract their protein coding genes because that's what I want to characterise here and to keep those that belong to an orthogroup created in Step 1. If thresholds are to be applied on the runs being included for analysis, I dinstinguished them as <studyID>.ortho.data.<threshold>.RData, eg., DRP000987.ortho.data.str.RData. The script for applying thresholds are not included in this repo at the moment and I will add an example later.

Step 4: After post-processing, Correlation_dataset.R can be run with the name of the dataset, eg DRP000987.ortho.data.str, as an input in the command line. This script performs 10000 permutation tests to get "p-values" as defined in "Identification of co-expressed genes via correlation techniques" section of https://www.biorxiv.org/content/10.1101/576116v5. The numerical scores obtained at the end of this script can be converted to "p-values" by doing (score + 1)/10000. For my needs, I have kept the original scores, meaning that a score of 0 refers to the most robustly inferred association/interaction.

Step 5: From these results, the bipartite part of the network can be extracted. I have chosen only those associations that had a score of 0. This takes the form of an edge table, with the first column as host orthogroups, second column for parasite orthogroups and the third column for the correlation strength. I will add an example. I call these files <studyID>\_<threshold>\_bipartite.RData, eg, DRP000987\_str\_bipartite.RData.

Step 6: Now, this bipartite network file can be used with topGO.R to perform GO analysis of the associated host and parasite genes.

More details of how each script works is mentioned within each script.
