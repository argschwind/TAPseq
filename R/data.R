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

#' Chromosome 11 sequence templates
#'
#' Sequence templates for primer design examples. The templates consist of chromosome 11 truncated
#' transcript sequences with the Drop-seq primer (polyA capture) added to the 3' ends.
#'
#' @format object of \code{DNAStringSet} class.
"chr11_sequence_templates"

#' Chromosome 11 primers
#'
#' Example of a \code{\link[TAPseq]{TsIOList}} object containing input and output for chromosome 11
#' genes primer design.
#'
#' @format object of \code{TsIOList} class.
"chr11_primers"

#' Mouse bone marrow 10x data
#'
#' Subset of a 10x mouse bone marrow dataset taken from Baccin et al., 2019.
#' (\url{https://www.nature.com/articles/s41556-019-0439-6}).
#'
#' @format object of \code{Seurat} class.
"bone_marrow_genex"
