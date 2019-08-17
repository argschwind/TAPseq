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

# change seqlevels from UCSC to Ensembl style...
seqnames <- seqnames(hg38)
seqnames <- sub("chr", "", seqnames)
seqnames[seqnames == "M"] <- "MT"
seqnames(hg38) <- seqnames

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

We have now successfully designed 5 outer primers per target gene. To select the best primer per
gene we could just pick the primer with the lowest penalty score. However we will now use BLAST to
try to estimate potential off-target priming for every primer. We will use this to then select the
best primer per gene, i.e. the primer with the fewest off-targets.

First we need to prepare a BLAST database containg the genome sequence as well as the sequences of
annotated transcripts. We then used this database to run BLAST to estimate off-target priming.
```
# download and import gencode hg38 annotations
url <- "ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.chr.gtf.gz"
annot <- import(url, format = "gtf")

# extract exon annotations for protein-coding genes to build transcripts
exons <- annot[annot$type == "exon" & annot$gene_biotype == "protein_coding"]
        
# create fasta file with transcripts and genome sequences. this will take a while and requires
# around 1GB free storage! here we put it into a temporary directory, but to avoid having to
# recreate this file everytime we design primers, we could save it somewhere else.
tmp_dir <- tempdir()
fasta_gt <- file.path(tmp_dir, "hg38_genome_and_transcripts.fasta.gz")
createFasta(genome = hg38, annot = exons, output_fasta = fasta_gt, compress = TRUE)

# create blast database
blastdb <- file.path(tmp_dir, "blastdb")
createBLASTDb(blastdb = blastdb, fasta = fasta_gt, compression = "gzip")

# now we can blast our outer primers against the create database, which contains all chromosome
# sequences and all annotated protein-coding transcripts.
outer_primers <- blastPrimers(outer_primers, blastdb = blastdb, annot = exons, max_mismatch = 0,
                              min_aligned = 0.75, tmpdir = tmp_dir)
                              
# the primers now contain the number of estimated off-targets
tascseq_primers(outer_primers$HBE1)
```


To finalize our set of outer primers, we want to choose the best primer per target gene. We use the
BLAST results to select primers with the fewest exonic, intronic and intergeni off-target hits (in
that order).
```
# select best primer per target gene
best_outer_primers <- pickPrimers(outer_primers, n = 1, by = "off_targets")
```


Now we need to repeat the same procedure to design a set of inner nested primers. We extract the PCR
amplicons of the outer primers and use those as sequence templates.
```
# get amplicons of outer primers and strip primer is from names
inner_templates <- pcr_products(best_outer_primers)
names(inner_templates) <- sub(".primer_left_\\d+", "", names(inner_templates))

# create new TsIO objects for inner primers, note the different product size
inner_primers <- TASCseqInput(inner_templates, reverse_primer = reverse_primer,
                              product_size_range = c(150, 300), primer_num_return = 5,
                              primer_opt_tm = 62, primer_min_tm = 57, primer_max_tm = 65)
                              
# design inner primers
inner_primers <- designPrimers(inner_primers)

# blast inner primers
inner_primers <- blastPrimers(inner_primers, blastdb = blastdb, annot = exons, max_mismatch = 0,
                              min_aligned = 0.75, tmpdir = tmp_dir)
                              
# pick best primer per target gene
best_inner_primers <- pickPrimers(inner_primers, n = 1, by = "off_targets")
```


Done! We succesfully designed TASC-seq outer and inner primer panels for our target genes.
```
tascseq_primers(best_outer_primers)
tascseq_primers(best_inner_primers)
```
