#' TASCseq: R-package to design primers for TASC-seq
#'
#' This package provides functions to select transcript isoforms and design PCR primers for
#' TASC-seq.
#'
#' @section Installation:
#' In order to use the full functionality, Primer3 needs to be installed and added to PATH.
#' Furthermore, the \code{primer3_config} directory containing important config files should be
#' located in the same directory as the \code{primer3_core} executable. If this is not practical,
#' all functions interacting woth Primer3 have arguments to specify the path to these files.
#'
#' For more information on installation see: \url{https://github.com/argschwind/TASCseq}.
#'
#' @docType package
#' @name TASCseq
#'
#' @import methods
#'
#' @import GenomicRanges
#' @import BiocGenerics
#'
#' @importFrom IRanges IRanges Views
#'
NULL
