#!/bin/bash

# With this script, the indices for the reference genomes can be made with STAR

hp=$1

echo "I am looking for $hp indices"

# Check if indices already exist
if [ -f ../Genomes/indices/$hp/genomeParameters.txt ]; then
  echo "Indexed files exist"
else
  echo "Indexing $hp.fa now..."
    # with gff3 files:
    #index="STAR --runThreadN 8 \
    # --runMode genomeGenerate \
    # --genomeDir /SAN/Plasmo_compare/Genomes/indices/$hp \
    # --genomeFastaFiles /SAN/Plasmo_compare/Genomes/fasta/{1} \
    # --sjdbGTFfile /SAN/Plasmo_compare/Genomes/annotation/{2} \
    # --sjdbGTFtagExonParentTranscript Parent \
    # --sjdbGTFtagExonParentGene Parent \
    # --limitGenomeGenerateRAM 210000000000" # 195GB=210000000000bytes

    # index
    # Got gtf files, you dont need to specify gene parent and transcript parent
    index="STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ../Genomes/indices/$hp \
    --genomeFastaFiles ../Genomes/fasta/{1} \
    --sjdbGTFfile ../Genomes/annotation/{2} \
    --limitGenomeGenerateRAM 210000000000"

    parallel --xapply $index {1} {2} ::: $hp.fa ::: $hp.gtf
fi

echo "Finsihed indexing $hp.fa"
