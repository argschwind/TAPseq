context("Design primers using Primer3 and parse output")

library(Biostrings)
library(GenomicRanges)

# chr11 truncated transcript sequences
data("chr11_truncated_txs_seq")
txs_seqs <- chr11_truncated_txs_seq[1:2]
txs_ids  <- names(txs_seqs)

# create TsIOList from sequence templates
tapseq_io <- TAPseqInput(target_sequences = txs_seqs, product_size_range = c(350, 500))

# test parsePrimer3Output() ------------------------------------------------------------------------

# raw Primer3 output
primer3_output <- c(
  "SEQUENCE_ID=AKIP1",
  paste0("SEQUENCE_TEMPLATE=", sequence_template(tapseq_io[[1]])),
  "SEQUENCE_PRIMER_REVCOMP=CTACACGACGCTCTTCCGATCT",
  "PRIMER_NUM_RETURN=2",
  "PRIMER_PICK_LEFT_PRIMER=1",
  "PRIMER_PICK_RIGHT_PRIMER=0",
  "SEQUENCE_EXCLUDED_REGION=0,499 650,356",
  "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/path/to/primer3_config/",
  "PRIMER_LEFT_NUM_RETURNED=2",
  "PRIMER_RIGHT_NUM_RETURNED=0",
  "PRIMER_INTERNAL_NUM_RETURNED=0",
  "PRIMER_PAIR_NUM_RETURNED=0",
  "PRIMER_LEFT_0_PENALTY=0.119155",
  "PRIMER_LEFT_0_SEQUENCE=AGACATCCCTTGGTCCTGGA",
  "PRIMER_LEFT_0=621,20",
  "PRIMER_LEFT_0_TM=59.881",
  "PRIMER_LEFT_0_GC_PERCENT=55.000",
  "PRIMER_LEFT_0_SELF_ANY_TH=0.00",
  "PRIMER_LEFT_0_SELF_END_TH=0.00",
  "PRIMER_LEFT_0_HAIRPIN_TH=34.58",
  "PRIMER_LEFT_0_END_STABILITY=3.8600",
  "PRIMER_LEFT_1_PENALTY=0.441533",
  "PRIMER_LEFT_1_SEQUENCE=GAGTCGAAGCTGCACATGTG",
  "PRIMER_LEFT_1=563,20",
  "PRIMER_LEFT_1_TM=59.558",
  "PRIMER_LEFT_1_GC_PERCENT=55.000",
  "PRIMER_LEFT_1_SELF_ANY_TH=21.83",
  "PRIMER_LEFT_1_SELF_END_TH=21.83",
  "PRIMER_LEFT_1_HAIRPIN_TH=0.00",
  "PRIMER_LEFT_1_END_STABILITY=3.2100",
  "=",
  "SEQUENCE_ID=ARFIP2",
  paste0("SEQUENCE_TEMPLATE=", sequence_template(tapseq_io[[2]])),
  "SEQUENCE_PRIMER_REVCOMP=CTACACGACGCTCTTCCGATCT",
  "PRIMER_NUM_RETURN=2",
  "PRIMER_PICK_LEFT_PRIMER=1",
  "PRIMER_PICK_RIGHT_PRIMER=0",
  "SEQUENCE_EXCLUDED_REGION=0,1463 1614,356",
  "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/path/to/primer3_config/",
  "PRIMER_LEFT_NUM_RETURNED=2",
  "PRIMER_RIGHT_NUM_RETURNED=0",
  "PRIMER_INTERNAL_NUM_RETURNED=0",
  "PRIMER_PAIR_NUM_RETURNED=0",
  "PRIMER_LEFT_0_PENALTY=0.035056",
  "PRIMER_LEFT_0_SEQUENCE=CAGGGTGTGGGAGATTGGAC",
  "PRIMER_LEFT_0=1468,20",
  "PRIMER_LEFT_0_TM=60.035",
  "PRIMER_LEFT_0_GC_PERCENT=60.000",
  "PRIMER_LEFT_0_SELF_ANY_TH=0.00",
  "PRIMER_LEFT_0_SELF_END_TH=0.00",
  "PRIMER_LEFT_0_HAIRPIN_TH=0.00",
  "PRIMER_LEFT_0_END_STABILITY=4.0200",
  "PRIMER_LEFT_1_PENALTY=0.106521",
  "PRIMER_LEFT_1_SEQUENCE=ACCAGTTTTGCCCACATTGC",
  "PRIMER_LEFT_1=1549,20",
  "PRIMER_LEFT_1_TM=59.893",
  "PRIMER_LEFT_1_GC_PERCENT=50.000",
  "PRIMER_LEFT_1_SELF_ANY_TH=0.00",
  "PRIMER_LEFT_1_SELF_END_TH=0.00",
  "PRIMER_LEFT_1_HAIRPIN_TH=0.00",
  "PRIMER_LEFT_1_END_STABILITY=3.5600",
  "="
)

test_that("parsePrimer3Output() parses Primer3 output correctly", {

  # parse raw Primer3 output and add to TsIOList object
  output <- parsePrimer3Output(tapseq_io, primer3_output)

  # expected parsed primers for the AKIP1
  expect_primers <- IRanges(start = c(622, 564), end = c(641, 583))
  names(expect_primers) <- paste0("AKIP1.primer_left_", 0:1)
  mcols(expect_primers) <- DataFrame(penalty = c(0.119155, 0.441533),
                                     sequence = c("AGACATCCCTTGGTCCTGGA", "GAGTCGAAGCTGCACATGTG"),
                                     tm = c(59.881, 59.558),
                                     gc_percent = c(55.000, 55.000),
                                     self_any_th = c(0.00, 21.83),
                                     self_end_th = c(0.00, 21.83),
                                     hairpin_th = c(34.58, 0.00),
                                     end_stability = c(3.8600, 3.2100))

  # test that parsed primers for the AKIP1 are correct
  primers <- tapseq_primers(output[[1]])
  expect_equivalent(class(primers), "IRanges")
  expect_length(primers, 2)
  expect_equal(primers, expect_primers)

  # expected pcr products
  seq_templ <- sequence_template(tapseq_io)
  expect_pcr_prod <- c(DNAStringSet(subseq(seq_templ[[1]], 622, length(seq_templ[[1]]))),
                       DNAStringSet(subseq(seq_templ[[1]], 564, length(seq_templ[[1]]))))
  names(expect_pcr_prod) <- paste0("AKIP1.primer_left_", 0:1)

  # test that inferred PCR products are correct
  pcr_prod <- pcr_products(output[[1]])

  # test that parsePrimer3Output() infers correct PCR products
  expect_equivalent(class(pcr_prod), "DNAStringSet")
  expect_length(pcr_prod, 2)
  expect_equal(as.character(pcr_prod[[1]]), as.character(expect_pcr_prod[[1]]))
  expect_equal(as.character(pcr_prod[[2]]), as.character(expect_pcr_prod[[2]]))
  expect_equal(names(pcr_prod), names(expect_pcr_prod))

})

test_that("parsePrimer3Output() handles errors correctly", {

  # change primers in Primer3 output so that one doesn't match the template and one matches 3 times
  primer3_output[14] <- "PRIMER_LEFT_0_SEQUENCE=CGACATCCCTTGGTCCTGGA"
  primer3_output[54] <- "PRIMER_LEFT_1_SEQUENCE=ACCAG"

  # parse modified Primer3 output
  msg <- capture_messages(output <- parsePrimer3Output(tapseq_io, primer3_output))

  # check that correct error messages are returned
  expect_match(msg, "Error in parse_primer\\(\\) for: AKIP1", all = FALSE)
  expect_match(msg, "Error: No matching sequence found for primer 'primer_left_0' in template!",
               all = FALSE)
  expect_match(msg, "Error in parse_primer\\(\\) for: ARFIP2", all = FALSE)
  expect_match(msg, "Error: More than 1 matching sequence found for primer 'primer_left_1'",
               all = FALSE)

  # check that correct primers and pcr products are returned
  primers <- tapseq_primers(output)
  expect_equivalent(class(primers), "IRanges")
  expect_length(primers, 2)
  expect_equal(names(primers), c("AKIP1.primer_left_1", "ARFIP2.primer_left_0"))

  pcr_prod <- pcr_products(output)
  expect_equivalent(class(pcr_prod), "DNAStringSet")
  expect_length(pcr_prod, 2)
  expect_equal(names(pcr_prod), c("AKIP1.primer_left_1", "ARFIP2.primer_left_0"))

})

# test designPrimers() -----------------------------------------------------------------------------

test_that("designPrimers() returns output in correct format", {

  # design primers
  output <- designPrimers(tapseq_io)
  primers <- tapseq_primers(output)
  expect_equivalent(class(primers), "IRanges")
  expect_length(primers, 10)
  expect_equal(names(primers),
               c(paste0("AKIP1.primer_left_", 0:4), paste0("ARFIP2.primer_left_", 0:4)))
  expect_equal(names(mcols(primers)),
               c("penalty", "sequence", "tm", "gc_percent", "self_any_th", "self_end_th",
                 "hairpin_th", "end_stability"))
  expect_equal(as.character(sapply(mcols(primers), class)),
               c("numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric",
                 "numeric"))

  # test that parsed primers are correct (equivalent data structure, not exact values)
  pcr_prod <- pcr_products(output)
  expect_equivalent(class(pcr_prod), "DNAStringSet")
  expect_length(pcr_prod, 10)
  expect_equal(names(pcr_prod),
               c(paste0("AKIP1.primer_left_", 0:4), paste0("ARFIP2.primer_left_", 0:4)))
})
