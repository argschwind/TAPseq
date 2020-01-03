library(TAPseq)
library(Biostrings)
library(GenomicRanges)

## get sequences of truncated chr11 trabscripts and add to package for examples

# truncated transcripts for chr11
data("chr11_truncated_txs")

# chr11 sequence loaded from fasta file
chr11_seq_fasta <- system.file("extdata", "chr11_sequence.fasta.gz", package = "TAPseq")
chr11_seq <- readDNAStringSet(chr11_seq_fasta)

# get sequences for all target transcripts in chr11 region
chr11_truncated_txs_seq <- getTxsSeq(chr11_truncated_txs, genome = chr11_seq)

# save data as RData files in data directory
usethis::use_data(chr11_truncated_txs_seq, overwrite = TRUE)
