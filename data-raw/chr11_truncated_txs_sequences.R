library(TAPseq)
library(BSgenome)

## get sequences of truncated chr11 trabscripts and add to package for examples

# truncated transcripts for chr11
data("chr11_truncated_txs")

# human genome (hg38) BSgenome object (needs to be installed separately from Bioconductor)
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# get sequences for all target transcripts in chr11 region
chr11_truncated_txs_seq <- getTxsSeq(chr11_truncated_txs, genome = hg38)

# save data as RData files in data directory
usethis::use_data(chr11_truncated_txs_seq, overwrite = TRUE)
