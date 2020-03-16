context("TsIO and TsIOList class creation")

library(Biostrings)
library(GenomicRanges)

# chr11 truncated transcript sequences
data("chr11_truncated_txs_seq")
txs_seqs <- chr11_truncated_txs_seq[1:2]
txs_ids  <- names(txs_seqs)

# chr11 truncated transcript annotations
data("chr11_truncated_txs")
txs_annot <- chr11_truncated_txs[1:2]

# 10x beads-oligo-dt sequence
beads_oligo <- "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

# reverse primer used in all PCR reactions
reverse_primer <- "CTACACGACGCTCTTCCGATCT"

# create TsIO objects with sequence inputs both as character strings and DNAString objects and
# optional target annotations
tsio_obj1 <- TsIO(target_sequence = as.character(txs_seqs[[1]]), beads_oligo = beads_oligo,
                  reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                  product_size_range = c(350, 500))
tsio_obj2 <- TsIO(target_sequence = txs_seqs[[2]], beads_oligo = DNAString(beads_oligo),
                  reverse_primer = DNAString(reverse_primer), sequence_id = txs_ids[2],
                  target_annot = chr11_truncated_txs[[2]], product_size_range = c(150, 300),
                  primer_num_return = 10, min_primer_region = 150, primer_opt_tm = 62,
                  primer_min_tm = 57, primer_max_tm = 65)

# create TsIOList object
tsiolist_obj <- TsIOList(tsio_obj1, tsio_obj2)

## test that correct output is returned ------------------------------------------------------------

test_that("TsIO class slots have correct format", {
  expect_true(is(tsio_obj1, "TsIO"))
  expect_true(is(sequence_id(tsio_obj1), "character"))
  expect_true(is(target_sequence(tsio_obj1), "DNAString"))
  expect_true(is(target_sequence(tsio_obj2), "DNAString"))
  expect_true(is(beads_oligo(tsio_obj1), "DNAString"))
  expect_true(is(beads_oligo(tsio_obj2), "DNAString"))
  expect_true(is(reverse_primer(tsio_obj1), "DNAString"))
  expect_true(is(reverse_primer(tsio_obj2), "DNAString"))
  expect_true(is.integer(product_size_range(tsio_obj1)))
  expect_length(product_size_range(tsio_obj1), 2)
  expect_true(is(target_annot(tsio_obj1), "GRanges"))
  expect_true(is(target_annot(tsio_obj2), "GRanges"))
  expect_true(is(primer_num_return(tsio_obj1), "integer"))
  expect_true(is(min_primer_region(tsio_obj1), "integer"))
  expect_true(is(primer_opt_tm(tsio_obj1), "integer"))
  expect_true(is(primer_min_tm(tsio_obj1), "integer"))
  expect_true(is(primer_max_tm(tsio_obj1), "integer"))
  expect_true(is(primer_num_return(tsio_obj2), "integer"))
  expect_true(is(min_primer_region(tsio_obj2), "integer"))
  expect_true(is(primer_opt_tm(tsio_obj2), "integer"))
  expect_true(is(primer_min_tm(tsio_obj2), "integer"))
  expect_true(is(primer_max_tm(tsio_obj2), "integer"))
})

test_that("TsIO class slots have correct values", {
  expect_equal(sequence_id(tsio_obj1), "AKIP1")
  expect_equal(as.character(target_sequence(tsio_obj1)), as.character(txs_seqs[[1]]))
  expect_equal(as.character(target_sequence(tsio_obj2)), as.character(txs_seqs[[2]]))
  expect_equal(as.character(beads_oligo(tsio_obj1)), beads_oligo)
  expect_equal(as.character(beads_oligo(tsio_obj2)), beads_oligo)
  expect_equal(as.character(reverse_primer(tsio_obj1)), reverse_primer)
  expect_equal(as.character(reverse_primer(tsio_obj2)), reverse_primer)
  expect_equal(product_size_range(tsio_obj1), c(350, 500))
  expect_equal(target_annot(tsio_obj1), GRanges())
  expect_equal(target_annot(tsio_obj2), chr11_truncated_txs[[2]])
  expect_equal(primer_num_return(tsio_obj1), 5)
  expect_equal(primer_num_return(tsio_obj2), 10)
  expect_equal(min_primer_region(tsio_obj1), 100)
  expect_equal(min_primer_region(tsio_obj2), 150)
  expect_equal(primer_opt_tm(tsio_obj1), NA_integer_)
  expect_equal(primer_min_tm(tsio_obj1), NA_integer_)
  expect_equal(primer_max_tm(tsio_obj1), NA_integer_)
  expect_equal(primer_opt_tm(tsio_obj2), 62)
  expect_equal(primer_min_tm(tsio_obj2), 57)
  expect_equal(primer_max_tm(tsio_obj2), 65)
})

test_that("TsIOList class data has correct format", {
  expect_length(tsiolist_obj, 2)
  expect_true(is(tsiolist_obj, "TsIOList"))
  expect_true(is(tsiolist_obj[[1]], "TsIO"))
  expect_true(is(tsiolist_obj[[2]], "TsIO"))
})

test_that("TsIOList class data has correct values", {
  expect_equal(tsiolist_obj[[1]], tsio_obj1)
  expect_equal(tsiolist_obj[[2]], tsio_obj2)
})

test_that("sequence_template returns correct output", {

  # from TsIO object
  beads_oligo_revcomp <- reverseComplement(DNAString(beads_oligo))
  expect_out1 <- as.character(xscat(txs_seqs[[1]], beads_oligo_revcomp))
  expect_equal(as.character(sequence_template(tsio_obj1)), expect_out1)

  # from TsIOList object
  expect_out2 <- as.character(xscat(txs_seqs, beads_oligo_revcomp))
  expect_equal(as.character(sequence_template(tsiolist_obj)), expect_out2)

})

test_that("other TsIOList accessor functions return correct data", {
  expect_equal(sequence_id(tsiolist_obj), txs_ids)
  target_seq <- target_sequence(tsiolist_obj)
  expect_true(is(target_seq, "DNAStringSet"))
  expect_length(target_seq, 2)
  expect_equal(as.character(target_seq[[1]]), as.character(txs_seqs[[1]]))
  expect_equal(as.character(target_seq[[2]]), as.character(txs_seqs[[2]]))
  annot <- target_annot(tsiolist_obj)
  expect_true(is(annot, "GRangesList"))
  expect_length(annot, 2)
  expect_equal(annot, GRangesList(GRanges(), chr11_truncated_txs[[2]]))
})

## test warnings and errors ------------------------------------------------------------------------

test_that("TsIO() raises errors/warnings with wrong input formats", {
  expect_error(TsIO(target_sequence = 1, beads_oligo = beads_oligo, reverse_primer = reverse_primer,
                    sequence_id = txs_ids[1], product_size_range = c(350, 500)),
               "no method or default for coercing \"numeric\" to \"DNAString\"")
  expect_error(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = 1,
                    reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                    product_size_range = c(350, 500)),
               "no method or default for coercing \"numeric\" to \"DNAString\"")
  expect_error(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo, reverse_primer = 1,
                    sequence_id = txs_ids[1], product_size_range = c(350, 500)),
               "no method or default for coercing \"numeric\" to \"DNAString\"")
  expect_error(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                    reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                    target_annot = txs_annot[[2]], product_size_range = c(350, 500)),
               "exons in target_annot are incompatible with target_sequence")
  expect_error(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                    reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                    target_annot = c(txs_annot[[2]], txs_annot[[2]][8]),
                    product_size_range = c(350, 500)),
               "overlapping exons found in target_annot")
  expect_error(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                    reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                    min_primer_region = 500, product_size_range = c(350, 500)),
               "product_size_range too narrow to allow min_primer_range")
  expect_error(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                    reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                    product_size_range = 500),
               "product_size_range needs to be an integer vector of length 2")
  expect_warning(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                      reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                      product_size_range = c("x", "y")),
                 "NAs introduced by coercion")
  expect_warning(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                      reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                      product_size_range = c(350, 500), primer_num_return = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                      reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                      product_size_range = c(350, 500), min_primer_region = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                      reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                      product_size_range = c(350, 500), primer_opt_tm = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                      reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                      product_size_range = c(350, 500), primer_min_tm = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(target_sequence = txs_seqs[[1]], beads_oligo = beads_oligo,
                      reverse_primer = reverse_primer, sequence_id = txs_ids[1],
                      product_size_range = c(350, 500), primer_max_tm = "a"),
                 "NAs introduced by coercion")
})

test_that("TsIOList aborts with wrong input formats", {
  expect_error(TsIOList(tsio_obj1, "a"), "must be a list containing TsIO objects")
  expect_error(TsIOList(tsio_obj1, 1), "must be a list containing TsIO objects")
})
