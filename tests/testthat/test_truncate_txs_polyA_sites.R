context("Truncate transcripts at polyA sites")

library(GenomicRanges)

# load input annotations
data("chr11_genes")
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)

# example polyA sites for these genes
data("chr11_polyA_sites")
chr11_polyA_sites <- chr11_polyA_sites[names(chr11_polyA_sites) != "ARFIP2"]

# get expected output
data("chr11_truncated_txs")
expect_out <- chr11_truncated_txs

# adjust expected output of ARFIP2 to no overlapping polyA site
expect_out[["ARFIP2"]] <- c(expect_out[["ARFIP2"]][1], expect_out[["ARFIP2"]])
mcols(expect_out[["ARFIP2"]])[1, "exon_number"] <- "9"
mcols(expect_out[["ARFIP2"]])[, "tag"] <- as.character(NA)
ranges(expect_out[["ARFIP2"]]) <- IRanges(
  start = c(6474683, 6477718, 6478041, 6478579, 6478738, 6479140, 6479972, 6480323, 6481231),
  end = c(6477268, 6477892, 6478198, 6478652, 6478959, 6479357, 6480068, 6480463, 6481479)
)


## test that correct output is returned ------------------------------------------------------------

# run tests
test_that("truncateTxsPolyA() returns correct output and format for GRanges", {
  output <- truncateTxsPolyA(target_genes[[1]], polyA_sites = chr11_polyA_sites)
  expect_true(is(output, "GRanges"))
  expect_equal(output, expect_out[[1]])
})

test_that("truncateTxsPolyA() returns correct output and format for GRangesList", {
  output <- truncateTxsPolyA(target_genes[2:4], polyA_sites = chr11_polyA_sites)
  expect_true(is(output, "GRangesList"))
  expect_equal(output, expect_out[2:4], check.attributes = FALSE)
})

test_that("truncateTxsPolyA() returns correct output if no overlap is found (1 tx as input)", {
  # input for one transcript that does not overlap with any polyA site
  input <- target_genes[[2]][mcols(target_genes[[2]])[, "transcript_id"] == "ENST00000614314.4"]
  output <- truncateTxsPolyA(input, polyA_sites = chr11_polyA_sites)
  expect_equal(input, output)
})

test_that("transcript_id = NULL is handeled correctly", {
  # create input for specific transcripts
  tx_ids <- c("ENST00000299576.9", "ENST00000614314.4")
  input <- endoapply(target_genes[1:2], FUN = function(x) {
    x[mcols(x)[, "transcript_id"] %in% tx_ids]
  })

  # when truncating with transcript_id set to NULL, a message is raised
  # GRanges input
  expect_message(
    outGR <- truncateTxsPolyA(input[[1]], polyA_sites = chr11_polyA_sites, transcript_id = NULL),
    "transcript_id = NULL, assuming all exons to be from same transcript."
  )

  # GRangesList input
  expect_message(
    outGRL <- truncateTxsPolyA(input, polyA_sites = chr11_polyA_sites, transcript_id = NULL),
    "transcript_id = NULL, assuming all exons to be from same transcript."
  )

  # all exons per gene will be treated as coming from the same transcript
  # GRanges input
  expect_output <- truncateTxsPolyA(input[[1]], polyA_sites = chr11_polyA_sites)
  expect_equal(outGR, expect_output)

  # GRangesList input
  expect_output <- truncateTxsPolyA(input, polyA_sites = chr11_polyA_sites)
  expect_equal(outGRL, expect_output)

  # if genes contain several transcripts that overlap, but transcript_id is set to NULL, this
  # aborts with an error message
  expect_error(
    truncateTxsPolyA(target_genes[[1]], polyA_sites = chr11_polyA_sites, transcript_id = NULL),
    "Overlapping exons found for transcripts:"
  )
  expect_error(
    truncateTxsPolyA(target_genes[1:2], polyA_sites = chr11_polyA_sites, transcript_id = NULL),
    "Overlapping exons found for transcripts:"
  )

})

## TO DO: tests for specific parameters (different extend_3prime_end, polyA_select, ignore_strand)

## test warnings and errors ------------------------------------------------------------------------

test_that("truncateTxPolyA() aborts if polyA_sites has wrong format", {
  expect_error(truncateTxsPolyA(target_genes[[1]], polyA_sites = GRangesList(chr11_polyA_sites)),
               "polyA_sites needs to be of class GRanges!")
  expect_error(truncateTxsPolyA(target_genes, polyA_sites = GRangesList(chr11_polyA_sites)),
               "polyA_sites needs to be of class GRanges!")
})

test_that("truncateTxPolyA() aborts if transcript_id column is not found", {
  expect_error(
    truncateTxsPolyA(target_genes[[1]], polyA_sites = chr11_polyA_sites, transcript_id = "tx_id"),
    "transcript_id column 'tx_id' not found!"
    )
  expect_error(
    truncateTxsPolyA(target_genes, polyA_sites = chr11_polyA_sites, transcript_id = "tx_id"),
    "transcript_id column 'tx_id' not found!"
    )
})

test_that("NA transcript id values are detected and warning is raised", {
  input <- target_genes[[1]]
  mcols(input)[2, "transcript_id"] <- as.character(NA)
  expect_warning(truncateTxsPolyA(input, polyA_sites = chr11_polyA_sites),
                 "NA transcript_id values found! These exons are dropped.")
  expect_warning(truncateTxsPolyA(GRangesList("AKIP1" = input), polyA_sites = chr11_polyA_sites),
                 "NA transcript_id values found! These exons will be dropped for:")
})

test_that("truncateTxsPolyA() aborts with overlapping exons", {
  input <- target_genes[[1]]
  start(input[2]) <- 8911123
  expect_error(truncateTxsPolyA(input, polyA_sites = chr11_polyA_sites),
               "Overlapping exons found for transcripts:")
  expect_error(truncateTxsPolyA(GRangesList("AKIP1" = input), polyA_sites = chr11_polyA_sites),
               "Overlapping exons found for transcripts:")
})

test_that("truncateTxsPolyA() aborts with conflicting strand information", {
  input <- target_genes[[1]]
  strand(input[5]) <- "-"
  expect_error(truncateTxsPolyA(input, polyA_sites = chr11_polyA_sites),
               "Conflicting strand information in provided transcripts!")
  expect_error(truncateTxsPolyA(GRangesList("AKIP1" = input), polyA_sites = chr11_polyA_sites),
               "Conflicting strand information in provided transcripts for:")
})

test_that("truncateTxsPolyA() aborts with chromosome information", {
  input <- target_genes[[1]]
  seqnames <- seqnames(input)
  seqnames[2] <- "chr12"
  seqnames(input) <- seqnames
  expect_error(truncateTxsPolyA(input, polyA_sites = chr11_polyA_sites),
               "Conflicting chromosome \\(seqname\\) information in provided transcripts!")
  expect_error(truncateTxsPolyA(GRangesList("AKIP1" = input), polyA_sites = chr11_polyA_sites),
               "Conflicting chromosome \\(seqname\\) information in provided transcripts for:")
})
