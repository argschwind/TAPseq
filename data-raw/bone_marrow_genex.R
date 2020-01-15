library(Seurat)

## create Seurat object containing cell population example data

# load full data set.
## THIS NEEDS TO BE REPLACED BY DOWNLOAD FROM URL!!
load("/g/steinmetz/gschwind/software/DEV/TAPseq_add/UseCells.Robj")
UseCells <- UpdateSeuratObject(UseCells)  # update from Seurat v2 to v3

# create new object
counts <- GetAssayData(UseCells, slot = "counts", assay = "RNA")
cell_idents <- Idents(UseCells)
object <- CreateSeuratObject(counts = counts)
Idents(object) <- cell_idents

# subsample cells to about 20% of cells (~1000 cells)
set.seed("20200115")
idents_split <- split(cell_idents, f = cell_idents)
idents_sampled <- lapply(idents_split, FUN = function(x) {
  sample(x, size = length(x) * 0.2)
})

# create vector with cell ids for these cells
names(idents_sampled) <- NULL
sampled_cells <- names(unlist(idents_sampled))

# subset object to these cells
bone_marrow_genex <- subset(object, cells = sampled_cells)

# save data as RData files in data directory
usethis::use_data(bone_marrow_genex, overwrite = TRUE)
