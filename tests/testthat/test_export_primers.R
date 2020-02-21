context("Export primers to data.frames or BED tracks")

# example primers
data("chr11_primers")
input <- pickPrimers(chr11_primers[1:2], n = 1, by = "off_targets")

# test primerDataFrame() ---------------------------------------------------------------------------

exp_out1 <- data.frame(seq_id = "AKIP1", seq_len = 924, start = 557, end = 576, primer_len = 20,
                       primer_id = "AKIP1.primer_left_4", penalty = 0.473386,
                       sequence = "CAGAGGCGAGTCGAAGCTGC", tm = 63.473, gc_percent = 65,
                       self_any_th = 18.61, self_end_th = 0, hairpin_th = 0, end_stability = 5.25,
                       pcr_product_size = 448, intergenic_off_targets = 0, intronic_off_targets = 0,
                       exonic_off_targets = 0, stringsAsFactors = FALSE)

exp_out2 <- data.frame(seq_id = "ARFIP2", seq_len = 1888, start = 1538, end = 1557, primer_len = 20,
                       primer_id = "ARFIP2.primer_left_2", penalty = 0.642407,
                       sequence = "CTGGGGCCTGACACCAGTTT", tm = 62.358, gc_percent = 60,
                       self_any_th = 14.98, self_end_th = 0, hairpin_th = 45.21,
                       end_stability = 2.66, pcr_product_size = 431, intergenic_off_targets = 4,
                       intronic_off_targets = 3, exonic_off_targets = 0, stringsAsFactors = FALSE)
exp_out2 <- rbind(exp_out1, exp_out2)

test_that("primerDataFrame() creates correct output", {
  out1 <- primerDataFrame(input[[1]])
  out2 <- primerDataFrame(input)
  expect_equal(out1, exp_out1)
  expect_equal(out2, exp_out2)
})

# test createPrimerTrack() -------------------------------------------------------------------------

exp_out_track1 <- data.frame(chrom = "chr11", start = 8914881, end = 8914901,
                             name = "AKIP1.primer_left_4", score = 0.473386, strand = "+",
                             thickStart = 8914881, thickEnd = 8914901, itemRgb = "0,0,0",
                             blockCount = 1, blockSizes = "20", blockStarts = "0",
                             stringsAsFactors = FALSE)

exp_out_track2 <- data.frame(chrom = "chr11", start = 6480404, end = 6480424,
                             name = "ARFIP2.primer_left_2", score = 0.642407, strand = "-",
                             thickStart = 6480404, thickEnd = 6480424, itemRgb = "0,0,0",
                             blockCount = 1, blockSizes = "20", blockStarts = "0",
                             stringsAsFactors = FALSE)
exp_out_track2 <- rbind(exp_out_track1, exp_out_track2)

test_that("createPrimerTrack() creates correct output", {
  out1 <- createPrimerTrack(input[[1]])
  out2 <- createPrimerTrack(input)
  expect_equal(out1, exp_out_track1)
  expect_equal(out2, exp_out_track2)
})

test_that("createPrimerTrack() aborts if target annotations are missing", {
  input_no_annot <- input
  target_annot(input_no_annot[[1]]) <- GRanges()
  expect_error(createPrimerTrack(input_no_annot[[1]]),
               "Input object does not contain valid target annotations")
  expect_error(createPrimerTrack(input_no_annot),
               "Not all TsIO objects in input valid contain target annotations")
})
