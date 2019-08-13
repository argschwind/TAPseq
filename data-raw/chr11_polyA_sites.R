library(TASCseq)
library(GenomicRanges)

## infer polyA sites from example drop-seq data and add to package for examples

# protein-coding exons of genes within chr11 region
data("chr11_genes")
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)

# bam file containing aligned Drop-seq reads
dropseq_bam <- system.file("extdata", "chr11_k562_dropseq.bam", package = "TASCseq")

# infer polyA sites
chr11_polyA_sites <- inferPolyASites(target_genes, bam = dropseq_bam, polyA_downstream = 50,
                                     wdsize = 100, min_cvrg = 1)

# save data as RData files in data directory
usethis::use_data(chr11_polyA_sites, overwrite = TRUE)
