library(TAPseq)
library(Biostrings)

## create sequence templates for truncated transcripts of target genes within chr11 region

# truncated transcripts for chr11
data("chr11_truncated_txs")

# chr11 sequence loaded from fasta file
chr11_seq_fasta <- system.file("extdata", "chr11_sequence.fasta.gz", package = "TAPseq")
chr11_seq <- readDNAStringSet(chr11_seq_fasta)

# create sequence templates for truncated transcripts
chr11_sequence_templates <- TAPseqSeqTemplates(chr11_truncated_txs, genome = chr11_seq)

# save data as RData files in data directory
usethis::use_data(chr11_sequence_templates, overwrite = TRUE)
