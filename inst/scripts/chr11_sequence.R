library(BSgenome)
library(Biostrings)

## extract chr11 genome sequence for sequence template generation examples

# get hg38 human genome sequence
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# extract sequence for the relevant portion of chr11
chr11_region_end <- 25000000
chr11_seq <- hg38$chr11
chr11_seq <- subseq(chr11_seq, 1, chr11_region_end)

# convert to DNAStringSet
chr11_seq <- DNAStringSet(list("chr11" = chr11_seq))

# save to .fasta file
writeXStringSet(chr11_seq, filepath = "inst/extdata/chr11_sequence.fasta.gz", compress = TRUE)
