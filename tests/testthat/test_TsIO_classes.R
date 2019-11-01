context("TsIO and TsIOList classes")

library(Biostrings)

# chr11 sequence templates
data("chr11_sequence_templates")
txs_seqs <- chr11_sequence_templates[1:2]
txs_ids  <- names(txs_seqs)

# reverse primer used in all PCR reactions
reverse_primer <- "AAGCAGTGGTATCAACGCAGAGT"

# create TsIO objects with sequence input both as character strings and DNAString objects
tsio_obj1 <- TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                  reverse_primer = reverse_primer, product_size_range = c(350, 500))
tsio_obj2 <- TsIO(sequence_template = txs_seqs[[2]], sequence_id = txs_ids[2],
                  reverse_primer = DNAString(reverse_primer), product_size_range = c(150, 300),
                  primer_num_return = 10, min_primer_region = 200, primer_opt_tm = 62,
                  primer_min_tm = 57, primer_max_tm = 65)

# create TsIOList object
tsiolist_obj <- TsIOList(tsio_obj1, tsio_obj2)

## test that correct output is returned ------------------------------------------------------------

test_that("TsIO class data has correct format", {
  expect_equivalent(class(tsio_obj1), "TsIO")
  expect_equivalent(class(sequence_id(tsio_obj1)), "character")
  expect_equivalent(class(sequence_template(tsio_obj1)), "DNAString")
  expect_equivalent(class(sequence_template(tsio_obj2)), "DNAString")
  expect_equivalent(class(reverse_primer(tsio_obj1)), "DNAString")
  expect_equivalent(class(reverse_primer(tsio_obj2)), "DNAString")
  expect_equivalent(class(product_size_range(tsio_obj1)), "integer")
  expect_length(product_size_range(tsio_obj1), 2)
  expect_equivalent(class(primer_num_return(tsio_obj1)), "integer")
  expect_equivalent(class(min_primer_region(tsio_obj1)), "integer")
  expect_equivalent(class(primer_opt_tm(tsio_obj1)), "integer")
  expect_equivalent(class(primer_min_tm(tsio_obj1)), "integer")
  expect_equivalent(class(primer_max_tm(tsio_obj1)), "integer")
  expect_equivalent(class(primer_num_return(tsio_obj2)), "integer")
  expect_equivalent(class(min_primer_region(tsio_obj2)), "integer")
  expect_equivalent(class(primer_opt_tm(tsio_obj2)), "integer")
  expect_equivalent(class(primer_min_tm(tsio_obj2)), "integer")
  expect_equivalent(class(primer_max_tm(tsio_obj2)), "integer")
})

test_that("TsIO class data has correct values", {
  expect_equal(sequence_id(tsio_obj1), "AKIP1")
  expect_equal(as.character(sequence_template(tsio_obj1)), as.character(txs_seqs[[1]]))
  expect_equal(as.character(sequence_template(tsio_obj2)), as.character(txs_seqs[[2]]))
  expect_equal(as.character(reverse_primer(tsio_obj1)), reverse_primer)
  expect_equal(as.character(reverse_primer(tsio_obj2)), reverse_primer)
  expect_equal(product_size_range(tsio_obj1), c(350, 500))
  expect_equal(primer_num_return(tsio_obj1), 5)
  expect_equal(primer_num_return(tsio_obj2), 10)
  expect_equal(min_primer_region(tsio_obj1), 100)
  expect_equal(min_primer_region(tsio_obj2), 200)
  expect_equal(primer_opt_tm(tsio_obj1), NA_integer_)
  expect_equal(primer_min_tm(tsio_obj1), NA_integer_)
  expect_equal(primer_max_tm(tsio_obj1), NA_integer_)
  expect_equal(primer_opt_tm(tsio_obj2), 62)
  expect_equal(primer_min_tm(tsio_obj2), 57)
  expect_equal(primer_max_tm(tsio_obj2), 65)
})

test_that("TsIOList class data has correct format", {
  expect_length(tsiolist_obj, 2)
  expect_equivalent(class(tsiolist_obj), "TsIOList")
  expect_equivalent(class(tsiolist_obj[[1]]), "TsIO")
  expect_equivalent(class(tsiolist_obj[[2]]), "TsIO")
})

test_that("TsIOList class data has correct values", {
  expect_equal(tsiolist_obj[[1]], tsio_obj1)
  expect_equal(tsiolist_obj[[2]], tsio_obj2)
})

test_that("TsIOList accessor function return correct data", {
  seq_templ <- sequence_template(tsiolist_obj)
  expect_equivalent(class(seq_templ), "DNAStringSet")
  expect_length(seq_templ, 2)
  expect_equal(as.character(seq_templ[[1]]), as.character(txs_seqs[[1]]))
  expect_equal(as.character(seq_templ[[2]]), as.character(txs_seqs[[2]]))
})

## test warnings and errors ------------------------------------------------------------------------

test_that("TsIO() raises errors/warnings with wrong input formats", {
  expect_error(TsIO(sequence_template = 1, sequence_id = txs_ids[1],
                    reverse_primer = reverse_primer, product_size_range = c(350, 500)),
               "no method or default for coercing \"numeric\" to \"DNAString\"")
  expect_error(TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                    reverse_primer = 1, product_size_range = c(350, 500)),
               "no method or default for coercing \"numeric\" to \"DNAString\"")
  expect_error(TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                    reverse_primer = reverse_primer, product_size_range = c("x", "y")),
               "product_size_range must be a numeric vector of length 2!")
  expect_warning(TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                      reverse_primer = reverse_primer, product_size_range = c(350, 500),
                      primer_num_return = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                      reverse_primer = reverse_primer, product_size_range = c(350, 500),
                      min_primer_region = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                      reverse_primer = reverse_primer, product_size_range = c(350, 500),
                      primer_opt_tm = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                      reverse_primer = reverse_primer, product_size_range = c(350, 500),
                      primer_min_tm = "a"),
                 "NAs introduced by coercion")
  expect_warning(TsIO(sequence_template = as.character(txs_seqs[[1]]), sequence_id = txs_ids[1],
                      reverse_primer = reverse_primer, product_size_range = c(350, 500),
                      primer_max_tm = "a"),
                 "NAs introduced by coercion")
})

test_that("TsIOList aborts with wrong input formats", {
  expect_error(TsIOList(tsio_obj1, "a"), "all elements in 'x' must be TsIO objects")
  expect_error(TsIOList(tsio_obj1, 1), "all elements in 'x' must be TsIO objects")
})
