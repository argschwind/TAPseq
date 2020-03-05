library(TAPseq)
library(BSgenome)
library(rtracklayer)

## create internal data for vignette, where Primer3 and BLAST are not run when building the package

# design outer and inner primers -------------------------------------------------------------------

# selected target genes
data("chr11_genes")
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)
target_genes <- target_genes[18:27]

# chr11 truncated transcripts sequences and annotations for selected target genes
data("chr11_truncated_txs_seq")
data("chr11_truncated_txs")
data("chr11_genes")
chr11_truncated_txs_seq <- chr11_truncated_txs_seq[names(target_genes)]
chr11_truncated_txs <- chr11_truncated_txs[names(target_genes)]

# design outer primers for chr11 target genes
outer_primers <- TAPseqInput(chr11_truncated_txs_seq, product_size_range = c(350, 500),
                             target_annot = chr11_truncated_txs)
outer_primers <- designPrimers(outer_primers)

# design inner primers for chr11 target genes
inner_primers <- TAPseqInput(chr11_truncated_txs_seq, product_size_range = c(150, 300),
                             target_annot = chr11_truncated_txs)
inner_primers <- designPrimers(inner_primers)

# blast primers ------------------------------------------------------------------------------------

# human genome BSgenome object (needs to be istalled from Bioconductor)
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# chromosome 11 sequence
chr11_genome <- DNAStringSet(getSeq(hg38, "chr11"))
names(chr11_genome) <- "chr11"

# create blast database
blastdb <- file.path(tempdir(), "blastdb")
createBLASTDb(genome = chr11_genome, annot = unlist(target_genes), blastdb = blastdb)

# blast primers
outer_primers <- blastPrimers(outer_primers, blastdb = blastdb, max_mismatch = 0,
                              min_aligned = 0.75, primer_targets = "gene_name")
inner_primers <- blastPrimers(inner_primers, blastdb = blastdb, max_mismatch = 0,
                              min_aligned = 0.75, primer_targets = "gene_name")

# primer complementarity ---------------------------------------------------------------------------

# select best outer and inner primers based on number of off-targets
best_outer_primers <- pickPrimers(outer_primers, n = 1, by = "off_targets")
best_inner_primers <- pickPrimers(inner_primers, n = 1, by = "off_targets")

# calculate primer complementarity
outer_comp <- checkPrimers(best_outer_primers)
inner_comp <- checkPrimers(best_inner_primers)

# save data ----------------------------------------------------------------------------------------

# create list wit data used in vignette
vignette_data <- list(outer_primers = outer_primers, best_inner_primers = best_inner_primers,
                      outer_comp = outer_comp, inner_comp = inner_comp)

# save data as RData files in data directory
usethis::use_data(vignette_data, overwrite = TRUE, internal = TRUE)
