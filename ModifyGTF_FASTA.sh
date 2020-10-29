#! /bin/bash/

# This is only an example:
# The objective is to 
# 1. Remove all commented lines from the GTF files
# 2. Renaming the chromosomes/contigs to be the same in the fasta and gtf files.
# As an example, the P.vivax files could contain the prefix Pvi_chr for all the sequences/locations


# remove all commented lines from GTF file of Plasmodium vivax
grep -v '^#.' Pvivax.gtf > temp && mv temp Pvivax.gtf

# attach prefix to gtf
sed -e 's/^/Pvi_chr/' Pvivax.gtf > temp && mv temp Pvivax.gtf


# attach prefix to fasta
sed -e 's/>/>Pvi_chr/' Pvivax.fa > temp && mv temp Pvivax.fa

