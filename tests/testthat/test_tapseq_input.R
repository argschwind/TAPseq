context("TAPseq input from sequence templates and create boulder IO records")

library(Biostrings)

# sequence templates
data("chr11_sequence_templates")
seq_templ <- chr11_sequence_templates[1:2]

# reverse primer used in all PCR reactions
reverse_primer <- "AAGCAGTGGTATCAACGCAGAGT"

# create TsIOList from sequence templates
obj <- TAPseqInput(sequence_templates = seq_templ, reverse_primer = reverse_primer,
                   product_size_range = c(350, 500), primer_num_return = 10,
                   min_primer_region = 150, primer_opt_tm = 62, primer_min_tm = 57,
                   primer_max_tm = 65)

# test TAPseqInput() -------------------------------------------------------------------------------

test_that("TAPseqInput() output has correct format", {
  expect_equivalent(class(obj), "TsIOList")
  expect_equivalent(class(obj[[1]]), "TsIO")
})

test_that("TAPseqInput() adds correct sequence templates to output", {
  expect_length(obj, length(seq_templ))
  expect_equivalent(class(sequence_template(obj)), "DNAStringSet")
  expect_length(sequence_template(obj), 2)
  expect_equal(as.character(sequence_template(obj)[[1]]), as.character(seq_templ[[1]]))
  expect_equal(as.character(sequence_template(obj)[[2]]), as.character(seq_templ[[2]]))
  expect_equal(names(obj), names(seq_templ))
})

test_that("TAPseqInput() adds common arguments to output", {
  expect_equal(product_size_range(obj[[1]]), c(350L, 500L))
  expect_equal(primer_num_return(obj[[1]]), 10L)
  expect_equal(min_primer_region(obj[[1]]), 150L)
  expect_equal(primer_opt_tm(obj[[1]]), 62L)
  expect_equal(primer_min_tm(obj[[1]]), 57L)
  expect_equal(primer_max_tm(obj[[1]]), 65L)
})

test_that("TAPseqInput() aborts with incorrect input format", {
  expect_error(TAPseqInput(sequence_templates = seq_templ[[1]], reverse_primer = reverse_primer,
                           product_size_range = c(350, 500)),
               "sequence_templates must be a DNAStringSet object")
})

# test createIORecord() ----------------------------------------------------------------------------

# expected boulder IO records for TsIO objects
expect_out1 <- c(
  "SEQUENCE_ID=AKIP1",
  paste0("SEQUENCE_TEMPLATE=", as.character(seq_templ[[1]])),
  "SEQUENCE_PRIMER_REVCOMP=AAGCAGTGGTATCAACGCAGAGT",
  "PRIMER_NUM_RETURN=10",
  "PRIMER_OPT_TM=62",
  "PRIMER_MIN_TM=57",
  "PRIMER_MAX_TM=65",
  "PRIMER_PICK_LEFT_PRIMER=1",
  "PRIMER_PICK_RIGHT_PRIMER=0",
  "SEQUENCE_EXCLUDED_REGION=0,499 650,356",
  "="
)

expect_out2 <- c(expect_out1[1:7],
                 "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/path/to/primer3_config",
                 expect_out1[8:11])

# expected boulder IO record for a TsIOList object
expect_out3 <- c(
  expect_out1,
  "SEQUENCE_ID=ARFIP2",
  paste0("SEQUENCE_TEMPLATE=", as.character(seq_templ[[2]])),
  "SEQUENCE_PRIMER_REVCOMP=AAGCAGTGGTATCAACGCAGAGT",
  "PRIMER_NUM_RETURN=10",
  "PRIMER_OPT_TM=62",
  "PRIMER_MIN_TM=57",
  "PRIMER_MAX_TM=65",
  "PRIMER_PICK_LEFT_PRIMER=1",
  "PRIMER_PICK_RIGHT_PRIMER=0",
  "SEQUENCE_EXCLUDED_REGION=0,1463 1614,356",
  "="
)

test_that("createIORecord() creates correct output from TsIO object", {
  output <- createIORecord(obj[[1]])
  expect_equal(class(output), "character")
  expect_identical(output, expect_out1)
})

test_that("createIORecord() handles optional thermo_params_path parameter correctly", {
  output <- createIORecord(obj[[1]], thermo_params_path = "/path/to/primer3_config")
  expect_equal(class(output), "character")
  expect_identical(output, expect_out2)
})

test_that("createIORecord() creates correct output from TsIOList object", {
  output <- createIORecord(obj)
  expect_equal(class(output), "character")
  expect_identical(output, expect_out3)
})
