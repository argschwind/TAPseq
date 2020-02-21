#' Select target genes
#'
#' Select target genes that serve as markers for cell populations using a linear model with lasso
#' regularization. How well a selected set of target genes discriminates between cell populations
#' can be assessed in an intuitive way using UMAP visualization.
#'
#' @param object Seurat object containing single-cell RNA-seq data from which best marker genes for
#'   different cell populations should be learned. Needs to contain population identities for all
#'   cell.
#' @param expr_percentile Expression percentiles that candidate target genes need to fall into.
#'   Default is 60\% to 99\%, which excludes bottom 60\% and top 1\% expressed genes from markers.
#' @param targets Desired number of target genes. Approximately this many target genes will be
#'   returned. If set to NULL, the optimal number of target genes will be estimated using a
#'   cross-valdation approach. Warning: The number of target genes might end up being very large!
#' @param target_genes (character) Target gene names.
#' @param npcs (integer) Number of principal components to use for UMAP.
#' @return A character vector containing selected target gene identifiers.
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # example of mouse bone marrow 10x gene expression data
#' data("bone_marrow_genex")
#'
#' # identify approximately 100 target genes that can be used to identify cell populations
#' target_genes <- selectTargetGenes(bone_marrow_genex, targets = 100)
#'
#' # automatically identify the number of target genes to best identify cell populations using
#' # cross-validation. caution: this can lead to very large target gene panels!
#' target_genes_cv <- selectTargetGenes(bone_marrow_genex)
#'
#' # create UMAP plots to compare cell type identification based on full dataset and selected 100
#' # target genes
#' plotTargetGenes(bone_marrow_genex, target_genes = target_genes)
#' }
#' @name selectTargetGenes
NULL

#' @rdname selectTargetGenes
#' @export
selectTargetGenes <- function(object, targets = NULL, expr_percentile = c(0.6, 0.99)) {

  # check that Seurat and glmnet are installed
  if (requireNamespace("Seurat", quietly = TRUE) == FALSE) {
   stop("Seurat package must be installed for this optional method!")
  }

  if (requireNamespace("glmnet", quietly = TRUE) == FALSE) {
    stop("glmnet package must be installed for this optional method!")
  }

  # normalize raw data
  object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", assay = "RNA")

  # identify marker genes for different cell populations
  message("Finding all marker genes")
  all_markers <- Seurat::FindAllMarkers(object = object, only.pos = TRUE, assay = "RNA")

  # get genes that are within specified expression percentiles
  counts <- Seurat::GetAssayData(object, slot = "counts", assay = "RNA")
  avg_expr <- Matrix::rowMeans(counts)
  expr_quantiles <- stats::quantile(avg_expr, probs = expr_percentile)
  filt_genes  <- names(avg_expr[avg_expr >= expr_quantiles[1] & avg_expr <= expr_quantiles[2]])

  # extract data only on cell population markers within specified expression percentiles
  filt_markers <- all_markers[all_markers$gene %in% filt_genes, ]
  filt_object  <- subset(object, features = unique(filt_markers$gene))

  # extract gene expression and cell identities
  gene_expr <- Seurat::GetAssayData(filt_object, slot = "data", assay = "RNA")
  cell_pops <- Seurat::Idents(filt_object)

  # fit lasso glm
  if (is.null(targets)) {

    message("Fitting cross-validation model")
    cv_model <- glmnet::cv.glmnet(x = Matrix::t(gene_expr), y = cell_pops, family = "multinomial",
                                  alpha = 1)

    # get optimal lambda estimated through cv procedure and model with full data
    target_lambda <- cv_model$lambda.1se
    model <- cv_model$glmnet.fit

  } else {

    message("Fitting model")
    model <- glmnet::glmnet(x = Matrix::t(gene_expr), y = cell_pops, family = "multinomial",
                            alpha = 1)

    # define lambda that results in at least the requested number of targets
    target_df <- which.min(abs(targets - model$df))
    target_lambda <- model$lambda[target_df]
  }

  # get gene names with non-zero coefficients for the chosen lambda
  message("Extract target genes")
  coefs <- stats::coef(model, s = target_lambda)
  genes <- lapply(coefs, FUN = function(x) rownames(x)[x[, 1] != 0] )
  genes <- unique(unlist(genes))
  genes[genes != "(Intercept)"]

}

#' @rdname selectTargetGenes
#' @export
plotTargetGenes <- function(object, target_genes, npcs = 15) {

  # check that Seurat and cowplot are installed
  if (requireNamespace("Seurat", quietly = TRUE) == FALSE) {
    stop("Seurat package must be installed for this optional method!")
  }

  if (requireNamespace("cowplot", quietly = TRUE) == FALSE) {
    stop("cowplot package must be installed for this optional method!")
  }

  # extract data on target genes only
  targets_object <- subset(object, features = target_genes)

  # normalize raw data, find variable genes and scale for full dataset
  message("Prepare data")
  object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", assay = "RNA",
                                  verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, assay = "RNA", verbose = FALSE)
  object <- Seurat::ScaleData(object, assay = "RNA", verbose = FALSE)

  # normalize raw data, find variable genes and scale for target genes dataset
  targets_object <- Seurat::NormalizeData(targets_object, normalization.method = "LogNormalize",
                                  assay = "RNA", verbose = FALSE)
  targets_object <- Seurat::FindVariableFeatures(targets_object, assay = "RNA", verbose = FALSE)
  targets_object <- Seurat::ScaleData(targets_object, assay = "RNA", verbose = FALSE)

  # compute PCA
  message("Compute PCA")
  object <- Seurat::RunPCA(object, npcs = npcs * 3, assay = "RNA", verbose = FALSE)
  targets_object <- Seurat::RunPCA(targets_object, npcs = npcs * 3, assay = "RNA", approx = FALSE,
                                   verbose = FALSE)

  # compute UMAPs
  message("Compute UMAP")
  object <- Seurat::RunUMAP(object, reduction = "pca", dims = seq_len(npcs))
  targets_object <- Seurat::RunUMAP(targets_object, reduction = "pca", seq_len(npcs))

  # create UMAP plots
  p1 <- Seurat::DimPlot(object, reduction = "umap") +
    ggplot2::ggtitle("Full dataset")
  p2 <- Seurat::DimPlot(targets_object, reduction = "umap") +
    ggplot2::ggtitle(paste(length(target_genes), "target genes"))

  # create plot containing both umap plots
  legend <- cowplot::get_legend(p1)
  cowplot::plot_grid(p1 + ggplot2::theme(legend.position = "none"),
                     p2 + ggplot2::theme(legend.position = "none"),
                     legend,
                     axis = "b",
                     nrow = 1,
                     rel_widths = c(1, 1, 0.75))

}
