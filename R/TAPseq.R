#' TAPseq: R-package to design primers for TAP-seq
#'
#' This package provides functions to select transcript isoforms and design PCR primers for
#' TAP-seq.
#'
#' @section Installation:
#' In order to use the full functionality, Primer3 and BLAST need to be installed and added to PATH.
#' Furthermore, the \code{primer3_config} directory containing important files for Primer3 should be
#' located in the same directory as the \code{primer3_core} executable. If this is not practical,
#' all functions interacting with Primer3 have arguments to specify the paths to these files.
#'
#' For more information on installation see: \url{https://github.com/argschwind/TAPseq}.
#'
#' @docType package
#' @name TAPseq
#'
#' @import methods
#'
#' @import GenomicRanges
#' @import BiocGenerics
#' @import BSgenome
#'
#' @importFrom dplyr bind_rows left_join distinct
#' @importFrom tidyr separate pivot_wider
#'
#' @importFrom S4Vectors endoapply mendoapply queryHits subjectHits mcols DataFrame
#' @importFrom IRanges IRanges IRangesList Views
#' @importFrom Biostrings DNAString DNAStringSet start end subseq reverseComplement matchPattern
#' @importClassesFrom Biostrings DNAString DNAStringSet
#' @importClassesFrom IRanges IRanges
NULL
