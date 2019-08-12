# TASCseq
R package to design PCR primers for TASC-seq.

Let's prepare a target gene panel for primer design:
```
library(TASCseq)
library(GenomicRanges)

# gene annotations for chromosome 11 genomic region
data("chr11_genes")
target_genes <- chr11_genes

# convert to GRangesList containing annotated exons per gene
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)

# K562 Drop-seq read data (this is just a small example file within the R package)
dropseq_bam <- system.file("extdata", "chr11_k562_dropseq.bam", package = "TASCseq")

# infer polyA sites from Drop-seq data
polyA_sites <- infer_polyA_sites(target_genes, bam = dropseq_bam)

```
