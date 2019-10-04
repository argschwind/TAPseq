#' Chromosome 11 genes
#'
#' GENCODE exon annotations of all protein-coding genes within a genomic region on chromosome 11.
#' @format object of \code{GRangesList} class.
"chr11_genes"

#' Chromosome 11 polyA sites
#'
#' Example polyadenylation sites for expressed protein-coding genes in chromosome 11 genomic region.
#' This dataset was created using \code{\link[TAPseq]{inferPolyASites}} on available K562 Drop-seq
#' data. In a real world appliaction these sites would have to be pruned manually before further
#' use.
#'
#' @format object of \code{GRanges} class.
"chr11_polyA_sites"

#' Chromosome 11 sequence
#'
#' The human hg38 sequence of chromosome 11 as \code{\link[Biostrings]{DNAStringSet}} object. This
#' data is intended for examples for \code{\link[TAPseq]{getTxsSeq}}. In a real world application,
#' such an object would probably contain all chromosome sequences loaded from a .fasta file.
#'
#' @format object of \code{DNAStringSet} class.
"chr11_seq"
