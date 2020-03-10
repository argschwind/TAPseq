#' Chromosome 11 genes
#'
#' GENCODE exon annotations of all protein-coding genes within a genomic region on human chromosome
#' 11.
#' @format object of \code{\link[GenomicRanges:GRanges-class]{GRanges}} class.
"chr11_genes"

#' Chromosome 11 polyA sites
#'
#' Example polyadenylation sites for expressed protein-coding genes within human chromosome 11
#' genomic region. This dataset was created using \code{\link[TAPseq]{inferPolyASites}} on available
#' K562 Drop-seq data. In a real world appliaction these sites would have to be pruned manually
#' before further use.
#'
#' @format object of \code{\link[GenomicRanges:GRanges-class]{GRanges}} class.
"chr11_polyA_sites"

#' Chromosome 11 truncated transcripts
#'
#' Annotations of target gene transcripts within human chromosome 11 region that were truncated at
#' inferred polyA sites using \code{\link[TAPseq]{truncateTxsPolyA}}.
#'
#' @format object of \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} class.
"chr11_truncated_txs"

#' Chromosome 11 truncated transcript sequences
#'
#' Sequences of truncated transcripts within human chromosome 11 region that were extracted using
#' \code{\link[TAPseq]{getTxsSeq}}.
#'
#' @format object of \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} class.
"chr11_truncated_txs_seq"

#' Chromosome 11 primers
#'
#' Example of a \code{\link[TAPseq:TsIOList-class]{TsIOList}} object containing input and output for
#' chromosome 11 genes primer design.
#'
#' @format object of \code{\link[TAPseq:TsIOList-class]{TsIOList}} class.
"chr11_primers"

#' Mouse bone marrow 10x data
#'
#' Subset of a 10x mouse bone marrow dataset taken from Baccin et al., 2019
#' (\url{https://www.nature.com/articles/s41556-019-0439-6}). Contains gene expression and cell type
#' data for 362 cells.
#'
#' @format object of \code{\link[Seurat:Seurat-class]{Seurat}} class.
"bone_marrow_genex"
