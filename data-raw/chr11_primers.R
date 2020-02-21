library(TAPseq)
library(BSgenome)
library(rtracklayer)

## create TsIOList example data containing designed primers and blast hits

# design primers for chr11 target genes
data("chr11_truncated_txs_seq")
data("chr11_truncated_txs")
chr11_primers <- TAPseqInput(chr11_truncated_txs_seq, product_size_range = c(350, 500),
                             target_annot = chr11_truncated_txs)
chr11_primers <- designPrimers(chr11_primers)

# human genome (hg38) BSgenome object
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# get annotations for BLAST
annot_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"
annot <- import(annot_url, format = "gtf")
blast_exons <- annot[annot$type == "exon" & annot$gene_type == "protein_coding"]

# build BLAST database
blastdb <- file.path(tempdir(), "blastdb")
createBLASTDb(genome = hg38, annot = blast_exons, blastdb = blastdb)

# blast primers
chr11_primers <- blastPrimers(chr11_primers, blastdb = blastdb, max_mismatch = 0,
                              min_aligned = 0.75, primer_targets = "gene_name")

# save data as RData files in data directory
usethis::use_data(chr11_primers, overwrite = TRUE)
