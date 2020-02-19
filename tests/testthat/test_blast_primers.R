context("Estimate primer off-targets using BLAST")

library(Biostrings)
library(BSgenome)
library(GenomicRanges)

# test get_gt_sequences() --------------------------------------------------------------------------

# create DNAStringSet of containing 2 chromosomes of the hg38 genome
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
genome <- getSeq(hg38, c("chr11", "chrUn_KI270583v1"))

# get example exon annotations
data("chr11_truncated_txs")
exons <- unlist(chr11_truncated_txs[1:2])

# expected transcript sequences for first 2 truncated transcripts
truncated_txs_seq <- getTxsSeq(chr11_truncated_txs[1:2], genome = genome)

# expected names of output DNAStringSet object
expect_names <- c("lcl|chromosome;chr11",
                  "lcl|transcript;ENSG00000132254.12-MERGED;NA;ENSG00000132254.12;ARFIP2",
                  "lcl|transcript;ENSG00000166452.12-MERGED;NA;ENSG00000166452.12;AKIP1")

test_that("get_gt_sequences() obtains correct genome and transcript sequences", {
  gt_seqs <- get_gt_sequences(genome = genome, annot = exons)
  expect_true(is(gt_seqs, "DNAStringSet"))
  expect_length(gt_seqs, 3)
  expect_identical(names(gt_seqs), expect_names)
  expect_identical(as.character(gt_seqs[[1]]), as.character(getSeq(hg38, "chr11")))
  expect_identical(as.character(gt_seqs[[2]]), as.character(truncated_txs_seq[["ARFIP2"]]))
  expect_identical(as.character(gt_seqs[[3]]), as.character(truncated_txs_seq[["AKIP1"]]))
})

test_that("get_gt_sequences() with option standard_chromosomes = FALSE", {
  gt_seqs <- get_gt_sequences(genome = genome, annot = exons, standard_chromosomes = FALSE)
  expect_true(is(gt_seqs, "DNAStringSet"))
  expect_length(gt_seqs, 4)
  expect_identical(names(gt_seqs),
                   c(expect_names[1], "lcl|chromosome;chrUn_KI270583v1", expect_names[2:3]))
})

# test blastPrimers helper functions ---------------------------------------------------------------

# input BLAST hits
input <- data.frame(qaccver = "1-AKIP1.primer_left_0", qlen = 20,
           saccver = c("lcl|transcript;ENST00000534506.5;AKIP1-211;ENSG00000166452.12;AKIP1",
                       "lcl|chromosome;chr11", "lcl|chromosome;chr6", "lcl|chromosome;chr6",
                       "lcl|chromosome;chr16", "lcl|chromosome;chr16"),
           length = c(20, 20, 20, 15, 17, 15), mismatch = c(0, 1, 0, 0, 0, 0),
           qend = c(20, 20, 20, 18, 19, 20),
           sstart = c(375, 8917303, 90116893, 79172447, 63872962, 650328),
           send = c(394, 8917322, 90116874, 79172461, 63872978, 650342),
           stringsAsFactors = FALSE)

# input annot
input_annot <- GRanges(seqnames = c(rep("chr6", 9), rep("chr11", 5), rep("chr16", 41)),
                 ranges = IRanges(
                   start = c(90296480, 90271849, 90252513, 90206569, 90088961, 90008602, 89950270,
                             89938144, 89926528, 8911139, 8911444, 8914826, 8917287, 8919337,
                             649365, 649763, 649991, 650254, 650539, 650995, 651199, 651644, 651827,
                             652467, 653341, 653525, 653746, 655029, 655307, 655573, 655773, 656302,
                             656732, 657091, 657762, 658183, 658525, 658896, 659086, 659245, 660058,
                             660612, 661051, 661342, 661597, 661891, 662220, 662679, 665679, 665950,
                             666220, 666455, 666673, 666905, 667432),
                   end = c(90296843, 90271941, 90252590, 90206681, 90089109, 90008856, 89951862,
                           89938350, 89932890, 8911223, 8911671, 8914930, 8917367, 8920084, 649426,
                           649854, 650167, 650362, 650709, 651103, 651266, 651747, 652039, 652535,
                           653451, 653670, 653803, 655147, 655468, 655703, 655889, 656537, 656871,
                           657221, 657892, 658344, 658653, 659011, 659126, 659376, 660161, 660714,
                           661172, 661501, 661787, 662059, 662331, 662844, 665801, 666124, 666350,
                           666598, 666792, 666989, 667830)),
                 strand = c(rep("-", 9), rep("+", 46)),
                 transcript_id = c(rep("ENST00000257749.9", 9), rep("ENST00000299576.9", 5),
                                   rep("ENST00000293879.9", 41)),
                 transcript_name = c(rep("BACH2-201", 9), rep("AKIP1-201", 5),
                                     rep("WDR90-201", 41)),
                 gene_id = c(rep("ENSG00000112182.15", 9), rep("ENSG00000166452.12", 5),
                             rep("ENSG00000161996.19", 41)),
                 gene_name = c(rep("BACH2", 9), rep("AKIP1", 5), rep("WDR90", 41)))

# expected output
expect_out <- data.frame(qaccver = "1-AKIP1.primer_left_0", qlen = 20,
                subject_type = c(rep("chromosome", 5), "transcript"),
                id = c("chr11", "chr6", "chr6", "chr16", "chr16", "ENST00000534506.5"),
                name = c(rep(as.character(NA), 5), "AKIP1-211"), length = c(20, 20, 15, 17, 15, 20),
                mismatch = c(1, 0, 0, 0, 0, 0), qend = c(20, 20, 18, 19, 20, 20),
                sstart = c(8917303, 90116893, 79172447, 63872962, 650328, 375),
                send = c(8917322, 90116874, 79172461, 63872978, 650342, 394),
                hit_id = c("hit_2", "hit_3", "hit_4", "hit_5", "hit_6", "hit_1"),
                gene_id = c("ENSG00000166452.12", "ENSG00000112182.15", NA, NA,
                            "ENSG00000161996.19", "ENSG00000166452.12"),
                type = c("exonic", "intronic", "intergenic", "intergenic", "exonic", "exonic"),
                gene_name = c("AKIP1", "BACH2", NA, NA, "WDR90", "AKIP1"),
                stringsAsFactors = FALSE)

test_that("annotate_blast_hits() returns correct output", {
  out <- annotate_blast_hits(input, annot = input_annot)
  expect_true(is(out, "data.frame"))
  expect_identical(out, expect_out)
})

test_that("filter_3p_hits() filters BLAST hits correctly", {
  out1 <- filter_3p_hits(expect_out, max_mismatch = 0, min_aligned = 0.75)
  out2 <- filter_3p_hits(expect_out, max_mismatch = 1, min_aligned = 0.75)
  out3 <- filter_3p_hits(expect_out, max_mismatch = 0, min_aligned = 1)
  expect_identical(out1, expect_out[c(2, 5, 6), ])
  expect_identical(out2, expect_out[c(1, 2, 5, 6), ])
  expect_identical(out3, expect_out[c(2, 6), ])
})

test_that("get_primer_target_genes() correctly identifies primer targets", {
  out <- get_primer_target_genes(expect_out, annot = input_annot, primer_targets = "gene_name")
  expect_identical(out$target_gene_id, rep("ENSG00000166452.12", 6))
})

test_that("get_primer_target_genes() warns if primer targets can't be inferred from primer_id", {
  expect_out2 <- expect_out
  expect_out2$qaccver <- sub("AKIP1.", "", expect_out2$qaccver)
  expect_warning(out <- get_primer_target_genes(expect_out2, annot = input_annot,
                                                primer_targets = "gene_name"),
                 "Can't infer targets from primer IDs for \\(some\\) primers.")
  expect_identical(out$target_gene_id, rep(as.character(NA), 6))
})

test_that("get_primer_target_genes() warns if primer targets are not found in primer_targets", {
  expect_warning(out <- get_primer_target_genes(expect_out, annot = input_annot,
                                                primer_targets = "transcript_id"),
                 "Can't find primer targets for all primers with primer_targets set to")
  expect_identical(out$target_gene_id, rep(as.character(NA), 6))
})

test_that("count_primer_hits() counts primer (off-) targets correctly", {
 out <- count_primer_hits(expect_out)
 expect_out2 <- data.frame(intergenic_off_targets = 2, intronic_off_targets = 1,
                          exonic_off_targets = 2, row.names = "1-AKIP1.primer_left_0")
 expect_equal(out, expect_out2)
})
