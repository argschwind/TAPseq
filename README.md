# TASCseq
R package to design PCR primers for TASC-seq.

Let's create sequence templates for primer design for a target gene panel:
```
library(TASCseq)
library(GenomicRanges)
library(BiocParallel)

# gene annotations for chromosome 11 genomic region
data("chr11_genes")
target_genes <- chr11_genes

# convert to GRangesList containing annotated exons per gene
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)

# K562 Drop-seq read data (this is just a small example file within the R package)
dropseq_bam <- system.file("extdata", "chr11_k562_dropseq.bam", package = "TASCseq")

# register backend for parallelization
register(MulticoreParam(workers = 5))

# infer polyA sites from Drop-seq data
polyA_sites <- inferPolyASites(target_genes, bam = dropseq_bam, polyA_downstream = 50,
                               wdsize = 100, min_cvrg = 1)

# truncate transcripts at inferred polyA sites
truncated_txs <- truncateTxsPolyA(target_genes, polyA_sites = polyA_sites, parallel = TRUE)

# next, we extract the DNA sequences for the truncated transcripts

# the last thing we need to do to create sequence templates for primer desing is to add the Drop-seq
# primer at the 3' end
```
