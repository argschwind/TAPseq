library(Seurat)
library(Matrix)

## create Seurat object containing cell population example data

# download 10x data
url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41556-019-0439-6/MediaObjects/41556_2019_439_MOESM4_ESM.zip"
file <- tempfile(fileext = ".zip")
download.file(url, destfile = file)

# load full data set
load(unz(file, "NicheData10x.rda"))
NicheData10x <- UpdateSeuratObject(NicheData10x)

# filter for cell types to use
use_cells <- c("Ery/Mk prog.", "Erythroblasts", "Gran/Mono prog.", "LMPPs", "Mono prog.",
               "Neutro prog.", "Neutrophils", "pro-B")
NicheData10x_filt <- subset(NicheData10x, idents = use_cells)

# create new object without any meta data
counts <- GetAssayData(NicheData10x_filt, slot = "counts", assay = "RNA")
cell_idents <- Idents(NicheData10x_filt)
object <- CreateSeuratObject(counts = counts)
Idents(object) <- cell_idents

# get top 5% cells per population (~180cells)
n_txs <- colSums(object)
cell_idents <- cell_idents[names(sort(n_txs, decreasing = TRUE))]
idents_split <- split(cell_idents, f = cell_idents)
idents_top <- lapply(idents_split, FUN = function(x) {
  head(x, n = length(x) * 0.05)
})

# create vector with cell ids for these cells
names(idents_top) <- NULL
top_cells <- names(unlist(idents_top))

# subset object to these cells
bone_marrow_genex <- subset(object, cells = top_cells)

# remove any genes with less than 10 total transcripts
txs <- rowSums(GetAssayData(bone_marrow_genex))
bone_marrow_genex <- subset(bone_marrow_genex, features = names(txs[txs > 10]))

# save data as RData files in data directory
usethis::use_data(bone_marrow_genex, overwrite = TRUE)
