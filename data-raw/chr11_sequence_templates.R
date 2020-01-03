library(TAPseq)
library(Biostrings)

## create sequence templates for truncated transcripts of target genes within chr11 region

# truncated transcript sequences
data("chr11_truncated_txs_seq")

# create reverse complement of Drop-seq primer to add to 3' end of transcripts
ds_primer <- "TTTTTTTAAGCAGTGGTATCAACGCAGAGTACJJJJJJJJJJJJNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
ds_primer <- gsub(ds_primer, pattern = "J", replacement = "N")  # replace "J" by "N"
ds_primer <- reverseComplement(DNAString(ds_primer))

# add bead_seq to transcript sequences to create sequence templates
chr11_sequence_templates <- xscat(chr11_truncated_txs_seq, ds_primer)
names(chr11_sequence_templates) <- names(chr11_truncated_txs_seq)

# save data as RData files in data directory
usethis::use_data(chr11_sequence_templates, overwrite = TRUE)
