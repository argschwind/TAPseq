context("Get transcript sequences")

library(Biostrings)
library(GenomicRanges)

# truncated transcripts for chr11
data("chr11_truncated_txs")

# chr11 sequence loaded from fasta file
chr11_seq_fasta <- system.file("extdata", "chr11_sequence.fasta.gz", package = "TAPseq")
chr11_seq <- readDNAStringSet(chr11_seq_fasta)

# expected output
data("chr11_truncated_txs_seq")
expect_out <- chr11_truncated_txs_seq

# run tests
test_that("getTxsSeq extracts correct sequence from DNAStringSet object", {
  outGR <- getTxsSeq(chr11_truncated_txs[[1]], genome = chr11_seq)
  expect_equivalent(outGR, expect_out[[1]])
  outGRL <- getTxsSeq(chr11_truncated_txs, genome = chr11_seq)
  expect_equal(outGRL, expect_out)
})

test_that("getTxsSeq aborts if genome is incorrect class", {
  expect_error(getTxsSeq(chr11_truncated_txs, genome = chr11_seq[[1]]),
               "genome must be of class BSgenome or DNAStringSet!")
  expect_error(getTxsSeq(chr11_truncated_txs, genome = data.frame(1)),
               "genome must be of class BSgenome or DNAStringSet!")
})

test_that("getTxsSeq aborts if transcript strand contains '*'", {

  # create input with '*' in strand
  input <- chr11_truncated_txs[1:2]
  strand <- strand(input[[1]])
  strand[2] <- "*"
  strand(input[[1]]) <- strand

  # test error
  expect_error(getTxsSeq(input, genome = chr11_seq), "Incorrect strand information in transcripts!")
  expect_error(getTxsSeq(input[[1]], genome = chr11_seq),
               "Incorrect strand information in transcripts!")
})

test_that("getTxsSeq aborts if transcripts have inconsistent strand information", {

  # create input with conflicting strand
  input <- chr11_truncated_txs[1:2]
  strand <- strand(input[[1]])
  strand[2] <- "-"
  strand(input[[1]]) <- strand

  # test error
  expect_error(getTxsSeq(input, genome = chr11_seq), "Incorrect strand information in transcripts!")
  expect_error(getTxsSeq(input[[1]], genome = chr11_seq),
               "Incorrect strand information in transcripts!")
})
