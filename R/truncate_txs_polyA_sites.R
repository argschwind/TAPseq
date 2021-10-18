#' Truncate transcripts at polyA sites
#'
#' Truncate transcripts at overlapping polyadenylation (polyA) sites to infer  likely 3' ends of
#' transcripts. This is crucial to correctly design TAP-seq primers that amplify fragments of
#' specific lengths. Typically the exons of all annotated transcripts per target gene are provided
#' as input. If a polyA site overlaps a single transcript of a given gene, this transcript is
#' truncated and returned. In case a polyA site overlaps multiple transcripts of the same gene, a
#' "metatranscript" consisting of all annotated exons of the overlapping transcripts is generated
#' and truncated. No statements about expressed transcripts can be made if no overlapping polyA
#' sites are found for any transcripts of a gene. In that case a "meta transcript" consisting of
#' the merged exons of that gene is generated and returned.
#'
#' @param transcripts A \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#'   \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object containing  exons of the
#'   transcripts to be truncated. Transcripts for multiple genes can be provided as \code{GRanges}
#'   objects within a \code{GRangesList}.
#' @param polyA_sites A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object containing the
#'   polyA sites. This needs to contain a metadata entry names "score" if the option
#'   \code{polyA_select = "score"} is used. PolyA sites can be either obtained via running
#'   \code{\link[TAPseq]{inferPolyASites}} or imported from an existing .bed file
#'   (\code{\link[rtracklayer]{BEDFile}}).
#' @param extend_3prime_end Specifies how far (bp) 3' ends of transcripts should be extended when
#'   looking for overlapping polyA sites (default = 0). This enables capturing of polyA sites that
#'   occur downstream of annotated 3' ends.
#' @param polyA_select Specifies which eurisic should be used to select the polyA site used to
#'   truncate the transcripts if multiple overlapping polyA sites are found. By default
#'   \code{"downstream"} is used which choses the most downstream polyA site. \code{"score"} selects
#'   the polyA site with the highest score, which correspons to the read coverage when using
#'   \code{\link[TAPseq]{inferPolyASites}} to estimate polyA sites.
#' @param transcript_id (character) Name of the column in the metadata of \code{transcripts}
#'   providing transcript id for each exon (default: \code{"transcript_id"}). Set to \code{NULL} to
#'   ignore transcript ids and assume that all exons per gene belong to the same transcript.
#' @param gene_id,exon_number (character) Optional names of columns in metadata of
#'   \code{transcripts} containing gene id and exon number. These are only used to create new
#'   metadata when merging multiple transcripts into a meta transcript.
#' @param ignore_strand (logical) Specifies whether the strand of polyA sites should be ignored when
#'   looking for overlapping polyA sites. Default is \code{FALSE} and therefore only polyA sites on
#'   the same strand as the transcripts are considered. PolyA sites with strand \code{*} has the
#'   same effect as \code{ignore_strand = TRUE}.
#' @param parallel (logical) Triggers parallel computing using the \code{BiocParallel} package.
#'   This requires that a parallel back-end was registered prior to executing the function.
#'   (default: \code{FALSE}).
#' @return Either a \code{\link[GenomicRanges:GRanges-class]{GRanges}} or
#'   \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object containing the truncated
#'   transcripts.
#' @examples
#' library(GenomicRanges)
#'
#' # protein-coding exons of genes within chr11 region
#' data("chr11_genes")
#' target_genes <- split(chr11_genes, f = chr11_genes$gene_name)
#'
#' # only retain first 2 target genes, because truncating transcripts is currently computationally
#' # quite costly. try using BiocParallel for parallelization (see ?truncateTxsPolyA).
#' target_genes <- target_genes[1:2]
#'
#' # example polyA sites for these genes
#' data("chr11_polyA_sites")
#'
#' # truncate target genes at most downstream polyA site (default)
#' truncated_txs <- truncateTxsPolyA(target_genes, polyA_sites = chr11_polyA_sites)
#'
#' # change polyA selection to "score" (read coverage of polyA sites) and extend 3' end of target
#' # genes by 50 bp (see ?truncateTxsPolyA).
#' truncated_txs <- truncateTxsPolyA(target_genes, polyA_sites = chr11_polyA_sites,
#'                                   polyA_select = "score", extend_3prime_end = 50)
#' @export
setGeneric("truncateTxsPolyA",
           function(transcripts, polyA_sites, extend_3prime_end = 0,
                    polyA_select = c("downstream", "upstream", "score"),
                    transcript_id = "transcript_id", gene_id = "gene_id",
                    exon_number = "exon_number", ignore_strand = FALSE, parallel = FALSE)
             standardGeneric("truncateTxsPolyA")
)

#' @describeIn truncateTxsPolyA Truncate transcripts of one gene provided as \code{GRanges} object
#' @export
setMethod("truncateTxsPolyA", "GRanges",
          function(transcripts, polyA_sites, extend_3prime_end = 0,
                   polyA_select = c("downstream", "upstream", "score"),
                   transcript_id = "transcript_id", gene_id = "gene_id",
                   exon_number = "exon_number", ignore_strand = FALSE, parallel = FALSE) {

    message("Verifying input...")

    # abort if polyA_sites has wrong format
    if (!is(polyA_sites, "GRanges")) {
      stop("polyA_sites needs to be of class GRanges!", call. = FALSE)
    }

    # check for transcript_id column and add dummy transcript_id if set to NULL
    if (is.null(transcript_id)) {
      message("transcript_id = NULL, assuming all exons to be from same transcript.")
      tx_id_col <- paste0("dummy_tx_id_", floor(stats::runif(1, 1000000, 9999999)))
      mcols(transcripts)[[tx_id_col]] <- "transcript1"
    } else {
      missing_tx_ids <- check_missing_tx_id(transcripts, transcript_id = transcript_id)
      if (missing_tx_ids) {
        stop("transcript_id column '", transcript_id, "' not found!", call. = FALSE)
      } else {
        tx_id_col <- transcript_id
      }
    }

    # check for NA in transcript ids and raise warning if any are found
    na_tx_ids <- check_na_tx_id(transcripts, transcript_id = tx_id_col)
    if (na_tx_ids) {
      warning("NA transcript_id values found! These exons are dropped.", call. = FALSE)
    }

    # abort if transcripts have conflicting strand information
    if (check_conflict_strand(transcripts)) {
      stop("Conflicting strand information in provided transcripts!", call. = FALSE)
    }

    # abort if transcripts map to more than one chromosome
    if (check_conflict_chr(transcripts)) {
      stop("Conflicting chromosome (seqname) information in provided transcripts!", call. = FALSE)
    }

    # check for transcripts with overlapping exons and abort if any are found
    overlapping_exons <- check_overlap_exons(transcripts, transcript_id = tx_id_col)
    if (length(overlapping_exons) > 0) {
      stop("Overlapping exons found for transcripts:\n", print_max(overlapping_exons),
           call. = FALSE)
    }

    # if all input is fine, truncate transcripts at overlapping polyA sites
    message("Truncating transcripts...")
    output <- truncate_tx_polyA(transcripts, polyA_sites = polyA_sites,
                                extend_3prime_end = extend_3prime_end,
                                polyA_select = polyA_select, transcript_id = tx_id_col,
                                gene_id = gene_id,
                                exon_number = exon_number,
                                ignore_strand = ignore_strand)

    # remove dummy transcript id if added
    if (is.null(transcript_id)) {
      mcols(output)[tx_id_col] <- NULL
    }

    message("Done!")
    return(output)

  }
)

#' @describeIn truncateTxsPolyA Truncate transcripts of multiple genes provided as
#'   \code{GRangesList}
#' @export
setMethod("truncateTxsPolyA", "GRangesList",
          function(transcripts, polyA_sites, extend_3prime_end = 0,
                   polyA_select = c("downstream", "upstream", "score"),
                   transcript_id = "transcript_id", gene_id = "gene_id",
                   exon_number = "exon_number", ignore_strand = FALSE, parallel = FALSE) {

    message("Verifying input...")

    # abort if polyA_sites has wrong format
    if (!is(polyA_sites, "GRanges")) {
      stop("polyA_sites needs to be of class GRanges!", call. = FALSE)
    }

    # check for transcript_id column and add dummy transcript_id if set to NULL
    if (is.null(transcript_id)) {
      message("transcript_id = NULL, assuming all exons to be from same transcript.")
      tx_id_col <- paste0("dummy_tx_id_", floor(stats::runif(1, 1000000, 9999999)))
      transcripts <- endoapply(X = transcripts, FUN = function(x) {
        mcols(x)[[tx_id_col]] <- "transcript1"
        return(x)
      })
    } else {
      missing_tx_ids <- vapply(transcripts, FUN = check_missing_tx_id, FUN.VALUE = logical(1),
                               transcript_id = transcript_id)
      if (any(missing_tx_ids)) {
        stop("transcript_id column '", transcript_id, "' not found!", call. = FALSE)
      } else {
        tx_id_col <- transcript_id
      }
    }

    # check for NA in transcript ids
    na_tx_ids <- vapply(transcripts, FUN = check_na_tx_id, FUN.VALUE = logical(1),
                        transcript_id = tx_id_col)

    # raise warning if any NA transcript ids are found
    if (any(na_tx_ids)) {
      warning("NA transcript_id values found! These exons will be dropped for:\n",
              print_max(names(na_tx_ids)[na_tx_ids]), call. = FALSE)
    }

    # abort if transcripts have conflicting strand information
    conflicting_strand <- vapply(transcripts, FUN = check_conflict_strand, FUN.VALUE = logical(1))
    if (any(conflicting_strand)) {
      stop("Conflicting strand information in provided transcripts for:\n",
           names(conflicting_strand)[conflicting_strand], call. = FALSE)
    }

    # abort if transcripts map to more than one chromosome
    conflicting_chr <- vapply(transcripts, FUN = check_conflict_chr, FUN.VALUE = logical(1))
    if (any(conflicting_chr)) {
      stop("Conflicting chromosome (seqname) information in provided transcripts for:\n",
           names(conflicting_chr)[conflicting_chr], call. = FALSE)
    }

    # check for transcripts with overlapping exons
    if (parallel == TRUE) {
      overlapping_exons <- BiocParallel::bplapply(transcripts, FUN = check_overlap_exons,
                                                  transcript_id = tx_id_col)
    } else {
      overlapping_exons <- lapply(transcripts, FUN = check_overlap_exons,
                                  transcript_id = tx_id_col)
    }

    # convert to gene - transcript pairs for printing and abort if any were found
    overlapping_txs <- list_to_print(overlapping_exons)
    if (length(overlapping_txs) > 0) {
      stop("Overlapping exons found for transcripts:\n", print_max(overlapping_txs), call. = FALSE)
    }

    # if all input is fine, truncate transcripts at overlapping polyA sites
    message("Truncating transcripts...")
    if (parallel == TRUE) {
      output <- BiocParallel::bplapply(X = transcripts, FUN = truncate_tx_polyA,
                                       polyA_sites = polyA_sites, polyA_select = polyA_select,
                                       extend_3prime_end = extend_3prime_end,
                                       transcript_id = tx_id_col, gene_id = gene_id,
                                       exon_number = exon_number, ignore_strand = ignore_strand)

    } else {
      output <- lapply(X = transcripts, FUN = truncate_tx_polyA, polyA_sites = polyA_sites,
                       polyA_select = polyA_select, extend_3prime_end = extend_3prime_end,
                       transcript_id = tx_id_col, gene_id = gene_id, exon_number = exon_number,
                       ignore_strand = ignore_strand)
    }

    # convert ouptut to GRangesList, remove dummy transcript id if added
    output <- GenomicRanges::GRangesList(output)
    if (is.null(transcript_id)) {
      output <- endoapply(output, FUN = function(x) {
        mcols(x)[[tx_id_col]] <- NULL
        x
      })
    }

    message("Done!")
    return(output)

  }
)

## HELPER FUNCTIONS ================================================================================

# check for missing transcript id column
check_missing_tx_id <- function(txs, transcript_id) {
  if (is.null(mcols(txs)[[transcript_id]])) {
    TRUE
  } else {
    FALSE
  }
}

# check if a gene contains NA transcript ids (these transcripts will be dropped!)
check_na_tx_id <- function(txs, transcript_id) {
  if (any(is.na(mcols(txs)[[transcript_id]]))) {
    TRUE
  } else {
    FALSE
  }
}

# check if transcripts for a given gene have conflicting strand information (txs as GRanges object)
check_conflict_strand <- function(txs) {
  tx_strand <- unique(strand(txs))
  if (length(tx_strand) > 1) {
    TRUE
  } else {
    FALSE
  }
}

# check if transcripts map to more than one chromosome (seqnames)
check_conflict_chr <- function(txs) {
  tx_chrs <- unique(seqnames(txs))
  if (length(tx_chrs) > 1) {
    TRUE
  } else {
    FALSE
  }
}

# check for transcripts with overlapping exons and return their transcript ids
check_overlap_exons <- function(txs, transcript_id) {
  txs <- split(txs, f = mcols(txs)[[transcript_id]])
  overlap_exons <- vapply(txs, FUN = function(x) any(countOverlaps(x) != 1), FUN.VALUE = logical(1))
  names(overlap_exons)[overlap_exons]
}

# convert list of vectors to "name-value" pairs for printing
list_to_print <- function(x) {
  name_reps <- vapply(x, FUN = length, FUN.VALUE = integer(1))
  names <- rep(names(x), times = name_reps)
  values <- unlist(x)
  paste(names, values, sep = " - ")
}

# create a string of a maximumn of n elements from x, separated by new lines
print_max <- function(x, n = 19) {
  output <- paste0(utils::head(x, n), collapse = "\n")
  if (length(x) > n) {
    output <- paste(output, paste(length(x) - n, "more..."), sep = "\n")
  }
  return(output)
}
