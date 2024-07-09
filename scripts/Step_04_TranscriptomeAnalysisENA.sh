#!/bin/bash

<<'////' #starting multiline comment
#####################################################################
#
#   This is the main script that needs to be run in order to 
#   download sequence data, map them onto reference genomes and quantify gene expression
#   
#   To run this script, do:
#   sh Step_04_TranscriptomeAnalysisENA.sh <txt filename with all studies to analyse in a column>
#
#   A study ID looks like this: SRP123456 or ERP123465 or DRP124356
# 	To use this pipeline, the mapping tool STAR 2.7.5a was installed
#   Aspera and SRAtoolkit was installed for rapid downloading of the fastq files
#	In the last step, R will be required for counting, along with these packages: 
#	From Bioconductor: GenomicRanges, GenomicFeatures, GenomicAlignments and rtracklayer
#	Others: dplyr
#   
#   As of now, the files that need to be prepared in advance are:
#   1. <studyID>_accession.txt: metadata file from ENA website consisting of ftp links to the fastq files
#   2. runs_<studyID>.txt: text file with two columns - one for runs IDs and the other for studyIDs, separated by commas
#   3. positive_experiments: file with 3 columns - studyIDs, host sepcies (human, mouse or monkey) and parasite sepcies
#   
#	The genome fasta files and annotation files in GTF formatting are downloded from Ensembl
#	For a particular host-parasite system, for example, human infection with P. falciparum,
#	1. concatenate the fasta sequences, eg, cat human.fa Pfalciparum.fa > humanPfalciparim.fa
#	2. concatenate GTF files of the two organisms after removing the commented part from each. 
#	
#	The names of the chromosome / contigs in these concatenated files should be the same, 
#	otherwise the reads cannot be counted as falling in a particular lcation on the chromosome.
#	Make sure to download the same strains for each organism. Also check that they are named in the same way in the fa and gtf files by ensembl
# 
#   The indices of the the reference genomes can be created with STAR in advance and stored in ../Genomes/indices/$host$parasite/ relative to the scripts/ folder, e.g., ../Genomes/indices/humanPberghei/
#   Or, they can be dynamically created if they are not found 
#	Genomes/ folder has 3 subfolders: fasta/, annotation/ (for gtf files) and indices/
#
#	In parallel, OrthoFinder1.0.6 was used to get the single copy orthologous genes among the host and among the parasite species
#	To run OrthoFinder: ./orthofinder -f [path to folder with all parasite protein sequences] -t  6
#	(see post_process_countfile.R for more details)
#
#####################################################################
////

# Accept a list of studyIDs to be analysed as a .txt file arranged in a column

studies=$1
cat $studies >> studies_submitted_for_analysis.txt
cat $studies > current_studies_for_analysis.txt

st=$(cat $studies)

printf "Current studies for analysis are \n$st\n"

# to store all output files from the mapping step into one temp folder that can be deleted after the analysis and does not need to be backed up
mkdir ../tmp

for studyID in $st; do

	mkdir ../$studyID

	# As of now the runs_$studyID file has to be made manually, but I will try to put a python script here to extract this data from the ENA or SRA website
	mv runs\_$studyID ../$studyID/

	# The metadata is saved from ENA website in a file called $studyID_accession.txt which includes the ftp link for downloading the sequences in the 7th column
	# The next line selects the 7th column and puts in a file called download.txt within the study folder 
	# For paired end reads, the column provides two URLs per row separated by a space. For single end reads, there's just one URL in each row.

	awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' ../$studyID\_accession.txt | \
	cut -f7 | awk -F ";" 'OFS=" " {print $1, $2}' | \
	awk NF > ../$studyID/download.txt

 	# The links in the download.txt file can be used qith aspera and/or sratoolkit to download the fastq files
 	# downloading several files in parallel:

	parallel --eta -j 16 --link bash --verbose Step_04a_DownloadENAfastq.sh ::: $(cat ../$studyID/download.txt)

	# All the fastq files are downloaded first - a bottleneck(!). 
	# I will try to make it so that we don't have to wait all the files to be downloaded before the mapping can begin.
	# But might as well - the downloading takes a bit of time unless well-parallelised.
	# The script doallENA.sh maps the reads in a fastq file and sends the BAM files to an R script for gene expression quantification

  	parallel --eta -j 3 --link bash --verbose Step_04b_0_doallENA.sh ::: $(cat ../$studyID/runs\_$studyID.txt)
done

exit
