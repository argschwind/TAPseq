# TASCseq
R package to design PCR primers for TASC-seq.

Let's create sequence templates for primer design for a target gene panel. First we try to infer
likely polyadenylation (polyA) sites for the selected target genes. We then truncate target gene
transcripts at the polyA sites, as we assume that this is their 3' end.
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
                               wdsize = 100, min_cvrg = 1, parallel = TRUE)

# truncate transcripts at inferred polyA sites
truncated_txs <- truncateTxsPolyA(target_genes, polyA_sites = polyA_sites, parallel = TRUE)
```
To create sequence templates for primer design, we extract the genome sequence for the truncated
transcripts. Here we use Bioconductor's BSgenome, but the genome sequence could also be imported
from a fasta file into a DNAStringSet object.
```
library(BSgenome)

# human genome (hg38) BSgenome object (needs to be istalled separately from Bioconductor)
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# change chromosome names to ENSEMBL style...
seqnames(hg38) <- sub("chr", "", seqnames(hg38))

# get sequences for all target transcripts on chr11
tx_seqs <- getTxsSeq(truncated_txs, genome = hg38)

# the last thing we need to do to create sequence templates for primer desing is to add the Drop-seq
# primer at the 3' end

# create reverse complement of Drop-seq primer to add to 3' end of transcripts 
ds_primer <- "TTTTTTTAAGCAGTGGTATCAACGCAGAGTACJJJJJJJJJJJJNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
ds_primer <- gsub(ds_primer, pattern = "J", replacement = "N")  # replace "J" by "N"
ds_primer <- reverseComplement(DNAString(ds_primer))

# add bead_seq to transcript sequences to create sequence templates
seq_templates <- xscat(tx_seqs, ds_primer)
names(seq_templates) <- names(tx_seqs)
```

Primer3 uses boulder-IO records for input and output (see: http://primer3.org/manual.html). TASCseq
implements TsIO and TsIOList objects, which store Primer3 input and output for TASC-seq primer
design in R's S4 class system. They serve as the users interface to Primer3 during primer design.

First we create Primer3 input for outer forward primer design based on the generated sequence
tamplates. We provide the reverse primer used in all PCR reactions, so that Primer3 can pick primers
suitable to work with it.
```
# reverse primer used in all PCR reactions.
reverse_primer <- "AAGCAGTGGTATCAACGCAGAGT"

# create TASCseq IO for outer forward primers with similar melting temperatures like the rev. primer
outer_primers <- TASCseqInput(seq_templates, reverse_primer = reverse_primer,
                              product_size_range = c(350, 500), primer_num_return = 5,
                              primer_opt_tm = 62, primer_min_tm = 57, primer_max_tm = 65)
                              
# design 5 outer primers for each target gene
outer_primers <- designPrimers(outer_primers)

# the TsIO objects in outer_primers now contain the designed primers and expected amplicons
tascseq_primers(outer_primers$HBE1)
pcr_products(outer_primers$HBE1)

# these can also be accessed for all genes
tascseq_primers(outer_primers)
pcr_products(outer_primers)
```

We have now successfully designed 5 outer primers per target gene.
