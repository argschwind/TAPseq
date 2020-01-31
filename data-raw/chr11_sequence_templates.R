library(TAPseq)
library(Biostrings)
library(BSgenome)

## create sequence templates for truncated transcripts of target genes within chr11 region

# truncated transcripts for chr11
data("chr11_truncated_txs")

# hg38 BSgenome object
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# create sequence templates for truncated transcripts
chr11_sequence_templates <- TAPseqSeqTemplates(chr11_truncated_txs, genome = hg38)

# save data as RData files in data directory
usethis::use_data(chr11_sequence_templates, overwrite = TRUE)
