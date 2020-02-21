context("TAPseq input from sequence templates and create boulder IO records")

library(Biostrings)

# chr11 truncated transcript sequences
data("chr11_truncated_txs_seq")
txs_seqs <- chr11_truncated_txs_seq[1:2]
txs_ids  <- names(txs_seqs)

# chr11 truncated transcript annotations
data("chr11_truncated_txs")
txs_annot <- chr11_truncated_txs[1:2]

# Drop-seq beads-oligo-dt and right primer sequence
ds_oligo <- "TTTTTTTAAGCAGTGGTATCAACGCAGAGTACNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
ds_rev_primer <- "AAGCAGTGGTATCAACGCAGAGT"

# create TsIOList from target sequences with default arguments
obj1 <- TAPseqInput(txs_seqs, product_size_range = c(350, 500))

# create TsIOList from target sequences with altered arguments
obj2 <- TAPseqInput(target_sequences = txs_seqs, beads_oligo = ds_oligo,
                    reverse_primer = ds_rev_primer, product_size_range = c(150, 300),
                    target_annot = txs_annot, primer_num_return = 10, min_primer_region = 150,
                    primer_opt_tm = 62, primer_min_tm = 57, primer_max_tm = 65)

# test TAPseqInput() -------------------------------------------------------------------------------

test_that("TAPseqInput() output has correct format", {
  expect_true(is(obj1, "TsIOList"))
  expect_true(is(obj1[[1]], "TsIO"))
})

test_that("TAPseqInput() adds correct sequences to TsIOList", {
  expect_length(obj1, length(txs_seqs))
  expect_equal(as.character(target_sequence(obj1)[[1]]), as.character(txs_seqs[[1]]))
  expect_equal(as.character(target_sequence(obj1)[[2]]), as.character(txs_seqs[[2]]))
  expect_equal(names(obj1), names(txs_seqs))
})

test_that("TAPseqInput() adds correct default arguments", {
  expect_equal(product_size_range(obj1[[1]]), c(350L, 500L))
  expect_equal(primer_num_return(obj1[[1]]), 5L)
  expect_equal(min_primer_region(obj1[[1]]), 100L)
  expect_equal(primer_opt_tm(obj1[[1]]), 63L)
  expect_equal(primer_min_tm(obj1[[1]]), 59L)
  expect_equal(primer_max_tm(obj1[[1]]), 66L)
})

test_that("TAPseqInput() adds correct optional arguments", {
  expect_equal(ranges(target_annot(obj2)), ranges(txs_annot))
  expect_equal(as.character(beads_oligo(obj2[[2]])), ds_oligo)
  expect_equal(as.character(reverse_primer(obj2[[2]])), ds_rev_primer)
  expect_equal(product_size_range(obj2[[2]]), c(150L, 300L))
  expect_equal(primer_num_return(obj2[[2]]), 10L)
  expect_equal(min_primer_region(obj2[[2]]), 150L)
  expect_equal(primer_opt_tm(obj2[[2]]), 62L)
  expect_equal(primer_min_tm(obj2[[2]]), 57L)
  expect_equal(primer_max_tm(obj2[[2]]), 65L)
})

test_that("TAPseqInput() aborts with incorrect input format", {

  # wrong input format
  expect_error(TAPseqInput(target_sequences = txs_seqs[[1]], product_size_range = c(350, 500)),
    "sequence_templates must be a named DNAStringSet object")
  expect_error(TAPseqInput(target_sequences = txs_seqs, product_size_range = c(350, 500),
                           target_annot = txs_annot[1]),
    "target_annot must be a named GRangesList object of the same length as target_sequences")
  expect_error(TAPseqInput(target_sequences = txs_seqs, product_size_range = c(350, 500),
                           target_annot = c("a", "b")),
    "target_annot must be a named GRangesList object of the same length as target_sequences")

  # malformed input names
  txs_seqs_nn <- txs_seqs
  names(txs_seqs_nn) <- NULL
  expect_error(TAPseqInput(target_sequences = txs_seqs_nn, product_size_range = c(350, 500),
                           target_annot = txs_annot),
    "sequence_templates must be a named DNAStringSet object")

  txs_annot_mn <- txs_annot
  names(txs_annot_mn) <- c("AKIP1", "ARFIP1")
  expect_error(TAPseqInput(target_sequences = txs_seqs, product_size_range = c(350, 500),
                           target_annot = txs_annot_mn),
    "Names of target_annot and target_sequences are not the same")

  names(txs_annot_mn) <- NULL
  expect_error(TAPseqInput(target_sequences = txs_seqs, product_size_range = c(350, 500),
                           target_annot = txs_annot_mn),
               "Names of target_annot and target_sequences are not the same")

})

# test createIORecord() ----------------------------------------------------------------------------

# expected boulder IO records for TsIO objects
expect_out1 <- c(
  "SEQUENCE_ID=AKIP1",
  paste0("SEQUENCE_TEMPLATE=", as.character(sequence_template(obj1[[1]]))),
  "SEQUENCE_PRIMER_REVCOMP=CTACACGACGCTCTTCCGATCT",
  "PRIMER_NUM_RETURN=5",
  "PRIMER_OPT_TM=63",
  "PRIMER_MIN_TM=59",
  "PRIMER_MAX_TM=66",
  "PRIMER_PICK_LEFT_PRIMER=1",
  "PRIMER_PICK_RIGHT_PRIMER=0",
  "SEQUENCE_EXCLUDED_REGION=0,504 655,349",
  "="
)

# expected output with added thermo_params_path argument
expect_out2 <- c(expect_out1[1:10],
                 "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/path/to/primer3_config",
                 "=")

# expected boulder IO record for a TsIOList object
expect_out3 <- c(
  expect_out1,
  "SEQUENCE_ID=ARFIP2",
  paste0("SEQUENCE_TEMPLATE=", as.character(sequence_template(obj1[[2]]))),
  expect_out1[3:9],
  "SEQUENCE_EXCLUDED_REGION=0,1468 1619,349",
  "="
)

test_that("createIORecord() creates correct output from TsIO object", {
  output <- createIORecord(obj1[[1]])
  expect_true(is(output, "character"))
  expect_identical(output, expect_out1)
})

test_that("createIORecord() handles optional thermo_params_path parameter correctly", {
  output <- createIORecord(obj1[[1]], thermo_params_path = "/path/to/primer3_config")
  expect_true(is(output, "character"))
  expect_identical(output, expect_out2)
})

test_that("createIORecord() creates correct output from TsIOList object", {
  output <- createIORecord(obj1)
  expect_true(is(output, "character"))
  expect_identical(output, expect_out3)
})
