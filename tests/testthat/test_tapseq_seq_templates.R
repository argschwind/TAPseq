context("Create TAP-seq sequence templates from transcript annotations")

library(Biostrings)
library(BSgenome)

# truncated transcripts for chr11
data("chr11_truncated_txs")

# hg38 BSgenome object
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# create DNAStringSet of containing the sequence of chromosome 11 region
chr11_seq <- DNAStringSet(list("chr11" = getSeq(hg38, "chr11", end = 25000000)))

# chr11 sequence templates
data("chr11_sequence_templates")

# run tests ----------------------------------------------------------------------------------------

test_that("getTxsSeq extracts correct sequence from DNAStringSet object", {
  outGR <- TAPseqSeqTemplates(chr11_truncated_txs[[1]], genome = chr11_seq)
  expect_equal(as.character(outGR), as.character(chr11_sequence_templates[[1]]))
  outGRL <- TAPseqSeqTemplates(chr11_truncated_txs[1:2], genome = chr11_seq)
  expect_equal(as.character(outGRL), as.character(chr11_sequence_templates[1:2]))
})

test_that("getTxsSeq correctly handles beads_oligo argument", {
  out <- TAPseqSeqTemplates(chr11_truncated_txs[[1]], genome = chr11_seq, beads_oligo = "CCCCCCCC")
  out_3p <- as.character(subseq(out, 875, 932))
  expect_identical(out_3p, "GAATAAAGCCACTGTATGATTCTCTTAATAGCTATACATTAATCCTGTTTGGGGGGGG")
  expect_error(TAPseqSeqTemplates(chr11_truncated_txs[[1]], chr11_seq, beads_oligo = "NONSENSE"),
               "key 79 \\(char 'O'\\) not in lookup table")
  expect_error(TAPseqSeqTemplates(chr11_truncated_txs[[1]], chr11_seq, beads_oligo = NULL),
               "beads_oligo cannot be NULL")
})
