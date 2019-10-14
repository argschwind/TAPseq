library(TAPseq)
library(GenomicRanges)

## truncate transcripts at polyA sites and add to package for examples

# protein-coding exons of genes within chr11 region
data("chr11_genes")
target_genes <- split(chr11_genes, f = chr11_genes$gene_name)

# polyA sites for chr11 region
data("chr11_polyA_sites")

# truncate transcripts at inferred polyA sites
chr11_truncated_txs <- truncateTxsPolyA(target_genes, polyA_sites = chr11_polyA_sites,
                                        extend_3prime_end = 0, polyA_select = "downstream",
                                        ignore_strand = FALSE)

# save data as RData files in data directory
usethis::use_data(chr11_truncated_txs, overwrite = TRUE)
