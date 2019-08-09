library(rtracklayer)

## download GENCODE hg38 annotations, filter for exons of protein-coding genes in chr11 region

# genomic region coordinates on chromosome 11
chr11_region <- GRanges(seqnames = "11", ranges = IRanges(start = 155000, end = 23000000))

# import gencode hg38 annotations
annot_url <- "ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.chr.gtf.gz"
annot <- import(annot_url, format = "gtf")

# extract annotations within chr11 region
annot_region <- subsetByOverlaps(annot, ranges = chr11_region)

# only retain exons of protein-coding genes on chr21p
chr11_genes <- annot_region[annot_region$type == "exon" &
                              annot_region$gene_biotype == "protein_coding" &
                              annot_region$transcript_biotype == "protein_coding"
                            ]

# expressed genes (to be removed once gene selection functions are available)
expr_genes <- c("AKIP1", "ARFIP2", "BET1L", "C11orf58", "CARS", "CD151", "CD81", "CHID1", "COPB1",
  "CSNK2A3", "CTR9", "CTSD", "DEAF1", "EIF3F", "EIF4G2", "FANCF", "FAR1", "GTF2H1", "HRAS",
  "HTATIP2", "IFITM1", "IFITM2", "ILK", "IPO7", "LDHA", "MOB2", "MRPL17", "MRPL23", "NAP1L4",
  "NCR3LG1", "NUCB2", "NUP98", "PDE3B", "PGAP2", "PHRF1", "PKP3", "PNPLA2", "POLR2L", "PRMT3",
  "PSMA1", "PSMD13", "RASSF7", "RHOG", "RIC8A", "RNH1", "RPL27A", "RPLP2", "RPS13", "RRAS2", "RRM1",
  "RRP8", "SAAL1", "SIGIRR", "SIRT3", "SLC25A22", "SPTY2D1", "STIM1", "SVIP", "TAF10", "TALDO1",
  "TEAD1", "TIMM10B", "TMEM41B", "TMEM9B", "TOLLIP", "TPP1", "TRIM5", "TSG101", "TSPAN4", "TSSC4",
  "USP47", "WEE1", "ZBED5", "ZNF143", "ZNF195", "ZNF215", "HBB", "HBD", "HBE1", "HBG1", "HBG2")

# only select exons of expressed genes
chr11_genes <- chr11_genes[chr11_genes$gene_name %in% expr_genes]

# save data as RData files in data directory
usethis::use_data(chr11_genes, overwrite = TRUE)
