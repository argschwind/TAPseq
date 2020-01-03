#' Check primers for complementarity
#'
#' Check a TAP-seq primer set, i.e. outer or inner primers for a target gene panel, for potential
#' complementarity issues when multiplexing. Uses Primer3's \code{check_primers} functionality.
#'
#' @param object A \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object containing
#'   designed primers.
#' @param primer_opt_tm,primer_min_tm,primer_max_tm Optimal, minumum and maximum primer melting
#'   temperature.
#' @param thermo_params_path Optional path (character) to the \code{primer3_config} directory. Only
#'   required when using Primer3 < 2.5.0.
#' @param primer3_core Path (character) to the \code{primer3_core} executable. Usually this is
#'   inferred when loading/attaching the package.
#' @return A \code{\link[base]{data.frame}} with \code{check_primers} results.
#' @seealso \url{http://primer3.org/manual.html} for Primer3 manual.
#' @examples
#' \dontrun{
#' library(TAPseq)
#' library(ggplot2)
#'
#' # chr11 primers example data
#' data("chr11_primers")
#'
#' # pick best primers based on predicted off-targets
#' best_primers <- pickPrimers(chr11_primers, n = 1, by = "off_targets")
#'
#' # check for complementarity
#' comp <- checkPrimers(best_primers)
#'
#' # plot complementarity scores for every pair. the lines indicate complementarity scores of 47,
#' # the default value applied by Primer3 to identify high complementarity primer pairs
#' ggplot(comp, aes(x = primer_pair_compl_any_th, y = primer_pair_compl_end_th)) +
#'   geom_hline(aes(yintercept = 47), colour = "darkgray", linetype = "dashed") +
#'   geom_vline(aes(xintercept = 47), colour = "darkgray", linetype = "dashed") +
#'   geom_point(alpha = 0.25) +
#'   theme_bw()
#' }
#' @export
setGeneric("checkPrimers",
           function(object, primer_opt_tm = NA, primer_min_tm = NA, primer_max_tm = NA,
                    thermo_params_path = NA, primer3_core = getOption("TAPseq.primer3_core"))
             standardGeneric("checkPrimers")
)

#' @describeIn checkPrimers Check primers from \code{TsIO} objects.
#' @export
setMethod("checkPrimers", "TsIO", function(object, primer_opt_tm, primer_min_tm, primer_max_tm,
                                           thermo_params_path, primer3_core) {

  # extract primers
  primers <- tapseq_primers(object)

  # abort if no desinged primers are found
  if (length(primers) < 2) stop("At least 2 TAP-seq primers needed!", call. = FALSE)

  # check primers using Primer3
  check_primers(primers, primer_opt_tm = primer_opt_tm, primer_min_tm = primer_min_tm,
                primer_max_tm = primer_max_tm, thermo_params_path = thermo_params_path,
                primer3_core = primer3_core)

})

#' @describeIn checkPrimers Check primers from \code{TsIOList} objects.
#' @export
setMethod("checkPrimers", "TsIOList", function(object, primer_opt_tm, primer_min_tm, primer_max_tm,
                                               thermo_params_path, primer3_core) {

  # extract primers
  primers <- tapseq_primers(object)

  # abort if no desinged primers are found
  if (length(primers) < 2) stop("At least 2 TAP-seq primers needed!", call. = FALSE)

  # check primers using Primer3
  check_primers(primers, primer_opt_tm = primer_opt_tm, primer_min_tm = primer_min_tm,
                primer_max_tm = primer_max_tm, thermo_params_path = thermo_params_path,
                primer3_core = primer3_core)

})


# HELPER FUNCTIONS =================================================================================

# check primers for a set stored in an IRanges object (output of tapseq_primers())
check_primers <- function(primers, primer_opt_tm, primer_min_tm, primer_max_tm,
                          thermo_params_path, primer3_core) {

  # get all primer sequences
  primer_seqs <- structure(mcols(primers)$sequence, names = names(primers))

  # create all possible primer pairs
  pairs <- utils::combn(primer_seqs, m = 2, simplify = FALSE)

  # create IO record for all primer pairs and catch any errors or warnings
  message("Creating input for Primer3...")
  io <- primer_pairs_io(pairs, primer_opt_tm = primer_opt_tm, primer_min_tm = primer_min_tm,
                        primer_max_tm = primer_max_tm,
                        primer_thermodynamic_parameters_path = thermo_params_path)

  # design primers
  message("Running Primer3...")
  primer3_output <- system2(command = primer3_core, input = io, stdout = TRUE)

  # parse output into list
  message("Processing output...")
  primer3_output <- parse_primer3_output(primer3_output)

  # process output for every pair
  output <- lapply(primer3_output, FUN = process_output_record)
  output <- do.call(rbind, output)
  row.names(output) <- NULL

  # reurn output
  message("Done!")
  return(output)

}

# create boulder-io records for a list of primer pairs
primer_pairs_io <- function(primer_pairs, ...) {

  # create boulder io records for every pair and catch errors and warnings
  io <- lapply(primer_pairs, FUN = function(pair) {
    tryCatch({
      check_primers_io(pair, ...)
    }, error = function(e) {
      message("Error in check_primers_io() for primer pair: ")
      message(e, "")
      return(NULL)
    }, warning = function(w) {
      message("Warning in check_primers_io() for primer pair: ")
      message(w, "")
      return(NULL)
    })
  })

  # convert from list to vector
  unlist(io)

}

# create a boulder-io record for a primer pair
check_primers_io <- function(primer_pair, ...) {

  # id for primer pair
  pair_id <- paste(names(primer_pair), collapse = "-")

  # create reverse complement of first primer
  primer2_revcomp <- as.character(Biostrings::reverseComplement(DNAString(primer_pair[1])))

  # abort if primers are reverse complements
  if (primer_pair[2] == primer2_revcomp) {
    stop("Primers are reverse complements in pair: ", pair_id, call. = FALSE)
  }

  # create basic input list
  io <- list("sequence_id" = pair_id,
             "sequence_primer" = primer_pair[[1]],
             "sequence_primer_revcomp"= primer_pair[[2]],
             "primer_task" = "check_primers",
             "primer_explain_flag" = 1)

  # add any additional arguments
  io <- append(io, list(...))

  # remove any NA values
  io <- io[!is.na(io)]

  # transform to vector where each element contains one input line ("tag=value" format)
  io <- paste0(toupper(names(io)), "=", io)

  # add "=" record separator at the end of the record
  c(io, "=")

}

# process Primer3 check_pair output for one output record (1 pair)
process_output_record <- function(record) {

  # get pair id and split into individual primer ids
  pair_id <- record[["sequence_id"]]
  primer1 <- sub("(.+primer_left_\\d+)-(.+primer_left_\\d+)", "\\1", pair_id)
  primer2 <- sub("(.+primer_left_\\d+)-(.+primer_left_\\d+)", "\\2", pair_id)

  # get pair id and primer sequences
  primer1_seq <- record[["sequence_primer"]]
  primer2_seq <- record[["sequence_primer_revcomp"]]

  # extract results
  extract <- c("primer_left_0_penalty", "primer_right_0_penalty", "primer_pair_0_penalty",
               "primer_pair_0_compl_any_th", "primer_pair_0_compl_end_th")
  results <- structure(as.numeric(record[extract]), names = extract)

  # create output data.frame
  data.frame(primer1, primer2, primer1_seq, primer2_seq,
             primer1_penalty = results[["primer_left_0_penalty"]],
             primer2_penalty = results[["primer_right_0_penalty"]],
             primer_pair_penalty = results[["primer_pair_0_penalty"]],
             primer_pair_compl_any_th = results[["primer_pair_0_compl_any_th"]],
             primer_pair_compl_end_th = results[["primer_pair_0_compl_end_th"]],
             stringsAsFactors = FALSE)

}
