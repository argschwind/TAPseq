context("Infer polyA sites from Drop-seq/10X data")

library(GenomicRanges)

# load input annotations and get drop-seq bam file path
data("chr11_genes")
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)
dropseq_bam <- system.file("extdata", "chr11_k562_dropseq.bam", package = "TAPseq")

# get expected output
data("chr11_polyA_sites")
expect_out <- chr11_polyA_sites

# run tests
test_that("inferPolyASites() returns correct output and format", {
  output <- inferPolyASites(target_genes, bam = dropseq_bam, polyA_downstream = 50,
                           wdsize = 100, min_cvrg = 1)
  expect_true(is(output, "GRanges"))
  expect_length(output, length(expect_out))
  expect_equal(output, expect_out, check.attributes = FALSE)

})

test_that("inferPolyASites() aborts if input has wrong format", {
  expect_error(inferPolyASites(1, bam = dropseq_bam),
               "genes needs to be a GRangesList containing annotations per gene!")
  expect_error(inferPolyASites(chr11_genes, bam = dropseq_bam),
               "genes needs to be a GRangesList containing annotations per gene!")

})

test_that("inferPolyASites() aborts if strand information is incorrect", {
  strand(target_genes[[13]]) <- "*"
  expect_error(inferPolyASites(target_genes, bam = dropseq_bam),
               "Strand not '\\+' or '\\-' for at least 1 gene!")
})
