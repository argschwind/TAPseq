#' Chromosome 11 genes
#'
#' GENCODE exon annotations of all protein-coding genes within a genomic region on human chromosome
#' 11.
#' @format object of \code{GRanges} class.
"chr11_genes"

#' Chromosome 11 polyA sites
#'
#' Example polyadenylation sites for expressed protein-coding genes within human chromosome 11
#' genomic region. This dataset was created using \code{\link[TAPseq]{inferPolyASites}} on available
#' K562 Drop-seq data. In a real world appliaction these sites would have to be pruned manually
#' before further use.
#'
#' @format object of \code{GRanges} class.
"chr11_polyA_sites"

#' Chromosome 11 truncated transcripts
#'
#' Annotations of target gene transcripts within human chromosome 11 region that were truncated at
#' inferred polyA sites using \code{\link[TAPseq]{truncateTxsPolyA}}.
#'
#' @format object of \code{GRangesList} class.
"chr11_truncated_txs"

#' Chromosome 11 truncated transcript sequences
#'
#' Sequences of truncated transcripts within human chromosome 11 region that were extracted using
#' \code{\link[TAPseq]{getTxsSeq}}.
#'
#' @format object of \code{DNAStringSet} class.
"chr11_truncated_txs_seq"
