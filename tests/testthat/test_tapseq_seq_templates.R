context("Create TAP-seq sequence templates from transcript annotations")

library(Biostrings)

# truncated transcripts for chr11
data("chr11_truncated_txs")

# chr11 sequence loaded from fasta file
chr11_seq_fasta <- system.file("extdata", "chr11_sequence.fasta.gz", package = "TAPseq")
chr11_seq <- readDNAStringSet(chr11_seq_fasta)

# chr11 sequence templates
data("chr11_sequence_templates")

# run tests ----------------------------------------------------------------------------------------

test_that("getTxsSeq extracts correct sequence from DNAStringSet object", {
  outGR <- TAPseqSeqTemplates(chr11_truncated_txs[[1]], genome = chr11_seq)
  expect_equivalent(outGR, chr11_sequence_templates[[1]])
  outGRL <- TAPseqSeqTemplates(chr11_truncated_txs[1:2], genome = chr11_seq)
  expect_equal(outGRL, chr11_sequence_templates[1:2])
})

test_that("getTxsSeq correctly handles beads_oligo argument", {
  out <- TAPseqSeqTemplates(chr11_truncated_txs[[1]], genome = chr11_seq, beads_oligo = "CCCCCCCC")
  out_3p <- as.character(subseq(out, 875, 932))
  expect_identical(out_3p, "GAATAAAGCCACTGTATGATTCTCTTAATAGCTATACATTAATCCTGTTTGGGGGGGG")
  expect_error(TAPseqSeqTemplates(chr11_truncated_txs[[1]], chr11_seq, beads_oligo = "NONSENSE"),
               "key 79 \\(char 'O'\\) not in lookup table")
})
