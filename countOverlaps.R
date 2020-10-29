require(GenomicRanges)
require(GenomicFeatures)
require(GenomicAlignments)
require(rtracklayer)

# take host, parasite, run and study names as arguments from cmd
options(echo=TRUE)
args <- commandArgs(TRUE)
host <- args[1]
para <- args[2]
run <- args[3]
study <- args[4]

# Read GTF annotation file, eg, humanPberghei.gtf
gff3 <- rtracklayer::import(paste0("Genomes/annotation/", host, para, ".gtf",collapse=''), format = "gtf")

# Keep only exonic regions in the GTF files
gff3<- gff3[gff3$type%in%"exon",]

# Make GRanges object
gff3 <- GenomicRanges::makeGRangesFromDataFrame(gff3)

setMethod("seqinfo", "BamFile", function (x) {
 h <- scanBamHeader(x, what = "targets")[["targets"]]
 h <- h[!duplicated(names(h))]
 Seqinfo(names(h), h)
})

# Read the BAM file for the run
bamFile <- readGAlignments(paste0("tmp/",run,"Aligned.sortedByCoord.out.bam", collapse=''))

# count the number of times a read overlaps with an exonic region
counts <- countOverlaps(gff3, bamFile)

# Put genome location and count data into a data frame
counts_df <- data.frame(gff3, counts)

write.table(counts_df, file=paste0(study,"/countWithGFF3_",run,".txt"), sep='\t', row.names = FALSE)