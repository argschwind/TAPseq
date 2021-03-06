context("Get transcript sequences")

library(Biostrings)
library(GenomicFeatures)
library(BSgenome)

# truncated transcripts for chr11
data("chr11_truncated_txs")

# hg38 BSgenome object
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# create DNAStringSet of containing the sequence of chromosome 11 region
chr11_seq <- DNAStringSet(list("chr11" = getSeq(hg38, "chr11", end = 25000000)))

# expected output for the first 2 truncated transcripts --------------------------------------------

txs <- chr11_truncated_txs[1:2]

# get indices of transcripts on positive and negative strand
txs_pos <- which(all(strand(txs) == "+"))
txs_neg <- which(all(strand(txs) == "-"))

# order exons of each transcript according to order in transcript (5' -> 3')
txs[txs_pos] <- sort(txs[txs_pos], decreasing = FALSE)
txs[txs_neg] <- sort(txs[txs_neg], decreasing = TRUE)

# get sequences of transcripts
expect_out <- extractTranscriptSeqs(chr11_seq, transcripts = txs)

# run tests ----------------------------------------------------------------------------------------

test_that("getTxsSeq extracts correct sequence from DNAStringSet object", {
  outGR <- getTxsSeq(chr11_truncated_txs[[1]], genome = chr11_seq)
  expect_equal(as.character(outGR), as.character(expect_out[[1]]))
  outGRL <- getTxsSeq(chr11_truncated_txs[1:2], genome = chr11_seq)
  expect_equal(as.character(outGRL), as.character(expect_out))
})

test_that("getTxsSeq extracts correct sequence from BSgenome object", {
  outGR <- getTxsSeq(chr11_truncated_txs[[1]], genome = hg38)
  expect_equal(as.character(outGR), as.character(expect_out[[1]]))
  outGRL <- getTxsSeq(chr11_truncated_txs[1:2], genome = hg38)
  expect_equal(as.character(outGRL), as.character(expect_out))
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

test_that("getTxsSeq aborts if chromosomes of transcripts are not found in genome", {
  faulty_genome <- chr11_seq
  names(faulty_genome) <- "11"
  expect_error(getTxsSeq(chr11_truncated_txs, genome = faulty_genome),
               "Not all chromosomes in transcripts found in genome object!")
})
