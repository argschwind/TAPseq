#' Create TAPseq input from target sequences
#'
#' This function creates input for TAP-seq primer design from a DNAStringSet containing the
#' target sequences (typically transcript sequences).
#'
#' @param target_sequences A named \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#'   containing all target sequences.
#' @param product_size_range Numerical vector of length 2 specifying the desired length of the
#'   resulting amplicons.
#' @param beads_oligo Beads-oligo-dT sequence for the used droplet sequencing protocol (10x,
#'   Drop-seq). If nothing is specified (\code{beads_oligo = NA}), the 10x V3 Beads-oligo-dT
#'   sequence is used. Can be changed if primers are for instance designed for Drop-seq. Any barcode
#'   bases need to be replaced by \code{N}.
#' @param reverse_primer Reverse primer sequence used for all PCR reactions. Default is the 10x
#'   primer sequence: \code{CTACACGACGCTCTTCCGATCT}.
#' @param target_annot (optional) A named \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
#'   object with transcript annotations in case the targets are transcripts of gene loci.
#'   If provided, each \code{\link[GenomicRanges:GRanges-class]{GRanges}} within the list
#'   should contain all exons of one targeted transcripts. Names need to be the same as for
#'   \code{target_sequences}.
#' @param primer_num_return How many forward primers should be designed? (default: 5)
#' @param min_primer_region Minimum sequence length required for primer design. Mostly relevant in
#'   case a sequence template is too short to allow the specified \code{product_size_range}.
#' @param primer_opt_tm,primer_min_tm,primer_max_tm Optimal, minumum and maximum primer melting
#'   temperature. Set to NA to use Primer3s default values.
#' @return \code{\link[TAPseq:TsIOList-class]{TsIOList}} object.
#' @examples
#' # chromosome 11 truncated transcript sequences and annotations
#' data("chr11_truncated_txs_seq")
#'
#' # create TsIOList object for primer design from target sequences
#' obj <- TAPseqInput(chr11_truncated_txs_seq, product_size_range = c(350, 500))
#' obj
#'
#' # transcript annotations can be added for optional genome browser tracks of designed primers
#' data("chr11_truncated_txs")
#' obj <- TAPseqInput(chr11_truncated_txs_seq, product_size_range = c(350, 500),
#'                    target_annot = chr11_truncated_txs)
#'
#' # create input for primer design with Drop-seq instead of default 10x
#' ds_oligo <- "TTTTTTTAAGCAGTGGTATCAACGCAGAGTACNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
#' ds_rev_primer <- "AAGCAGTGGTATCAACGCAGAGT"
#' ds_obj <- TAPseqInput(chr11_truncated_txs_seq, beads_oligo = ds_oligo,
#'                       reverse_primer = ds_rev_primer, product_size_range = c(350, 500),
#'                       primer_opt_tm = 62, primer_min_tm = 57, primer_max_tm = 65)
#' @export
TAPseqInput <- function(target_sequences, product_size_range,
                        beads_oligo = NA, reverse_primer = "CTACACGACGCTCTTCCGATCT",
                        target_annot = NULL, primer_num_return = 5, min_primer_region = 100,
                        primer_opt_tm = 63, primer_min_tm = 59, primer_max_tm = 66) {

  # make sure that sequence template has the right format
  if (!is(target_sequences, "DNAStringSet") | is.null(names(target_sequences))) {
    stop("sequence_templates must be a named DNAStringSet object", call. = FALSE)
  }

  # make sure that target_annot has the right format (if provided)
  if (!is.null(target_annot)) {

    # check that it is a names GRangesList
    if (!is(target_annot, "GRangesList") | length(target_annot) != length(target_sequences)) {
      stop("target_annot must be a named GRangesList object of the same length as target_sequences",
           call. = FALSE)
    }

    # check names are the same as target_sequences
    targets <- sort(names(target_sequences))
    if (!identical(sort(names(target_annot)), targets)) {
      stop("Names of target_annot and target_sequences are not the same", call. = FALSE)
    }

    # sort target_sequences and target_annots by name
    target_sequences <- target_sequences[targets]
    target_annot <- target_annot[targets]

  } else {

    # create empty list in case target_annot was not provided
    target_annot <- vector(mode = "list", length = length(target_sequences))

  }

  # parse beads_oligo argument
  if (is.na(beads_oligo)) beads_oligo <- get_beads_oligo()

  # create TsIO objects for every sequence template
  io <- mapply(FUN = TsIO,
               target_sequence = target_sequences,
               target_annot = target_annot,
               sequence_id = names(target_sequences),
               MoreArgs = list(
                 beads_oligo = beads_oligo, reverse_primer = reverse_primer,
                 product_size_range = product_size_range, primer_num_return = primer_num_return,
                 min_primer_region = min_primer_region, primer_opt_tm = primer_opt_tm,
                 primer_min_tm = primer_min_tm, primer_max_tm = primer_max_tm),
               SIMPLIFY = FALSE)

  # convert to TsIOList object
  TsIOList(io)

}

#' Create boulder IO record
#'
#' Takes a \code{\link[TAPseq:TsIO-class]{TsIO}} or \code{\link[TAPseq:TsIOList-class]{TsIOList}}
#' object and converts it into a boulder IO record for Primer3. Essentially it converts it into a
#' list of character vectors that each contain the tag and the value in the form: "TAG=VALUE". More
#' on this format can be found in the \href{http://primer3.org/manual.html}{Primer3 manual}.
#'
#' This function is usually not needed by the user, because functions such as
#' \code{\link[TAPseq]{designPrimers}} handle IO record generation. However, this function can for
#' instance be useful to generate IO records, write them to a file and pass them to Primer3 in the
#' conventional way.
#'
#' @param object TsIO of TsIOList object for which a Primer3 boulder IO record should be created.
#' @param thermo_params_path Optional path (character) to the \code{primer3_config} directory. Only
#'   required when using Primer3 < 2.5.0.
#' @return A character vector containing the lines of the IO record.
#' @seealso \url{http://primer3.org/manual.html} for Primer3 manual.
#' @examples
#' # chromosome 11 truncated transcript sequences
#' data("chr11_truncated_txs_seq")
#'
#' # create TsIOList object for primer desing from sequence templates
#' obj <- TAPseqInput(chr11_truncated_txs_seq, product_size_range = c(350, 500))
#'
#' # create boulder IO record
#' boulder_io <- createIORecord(obj)
#' head(boulder_io, 11)
#' @export
setGeneric("createIORecord",
           function(object, thermo_params_path = NA) standardGeneric("createIORecord")
)

#' @describeIn createIORecord Create IO record from \code{TsIO} objects.
#' @export
setMethod("createIORecord", "TsIO", function(object, thermo_params_path) {

  # slot names of values for Primer3 that do not need to be processed
  slot_names <- c("sequence_id", "reverse_primer", "primer_num_return", "primer_opt_tm",
                  "primer_min_tm", "primer_max_tm")

  # get values of all slots as one character vector
  io <- vapply(slot_names, FUN = function(s, x) {
    as.character(slot(x, s))
  }, x = object, FUN.VALUE = character(1))

  # add sequence_template
  io <- c(io[1], "sequence_template" = as.character(sequence_template(object)), io[-1])

  # rename reverse_primer to Primer3 argument name
  rev_index <- which(names(io) == "reverse_primer")
  names(io)[rev_index] <- "sequence_primer_revcomp"

  # add Primer3 parameters required to trigger design of only forward primers
  io["primer_pick_left_primer"] <- "1"
  io["primer_pick_right_primer"] <- "0"

  # calculate and add excluded regions
  excluded_regions <- create_excluded_regions(object)
  io["sequence_excluded_region"] <- excluded_regions

  # add thermodynamic parameters path to io
  io <- c(io, "primer_thermodynamic_parameters_path" = thermo_params_path)

  # remove any NA elements, because they should not be part of the primer3 input
  io <- io[!is.na(io)]

  # transform to vector where each element contains one input line ("tag=value" format)
  io <- paste0(toupper(names(io)), "=", io)

  # add "=" record separator at the end of the output vector
  c(io, "=")

})

#' @describeIn createIORecord Create IO record from \code{TsIO} objects.
#' @export
setMethod("createIORecord", "TsIOList", function(object, thermo_params_path) {

  # create IO records for every TsIO object
  io <- lapply(object, FUN = createIORecord, thermo_params_path = thermo_params_path)
  names(io) <- NULL

  # convert into vector
  unlist(io)

})

## HELPER FUNCTIONS ================================================================================

# get beads oligo sequence for standard droplet sequencing protocol. currently this function is only
# used internally and not exported via NAMESPACE
get_beads_oligo <- function(protocol = c("10x_v3", "10x_v2", "dropseq")) {
  protocol <- match.arg(protocol)
  switch(protocol,
    "10x_v3" = "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    "10x_v2" = "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    "dropseq" = "TTTTTTTAAGCAGTGGTATCAACGCAGAGTACNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  )
}

# create excluded regions that respect the specified product size range specified in a TsIO object
create_excluded_regions <- function(object) {

  # ------------------------------------------------------------------------------------------------
  # STEP 1:
  # match reverse primer to sequence and get start and stop positions of right primer binding site
  # ------------------------------------------------------------------------------------------------

  # sequence id and sequence template
  seq_id <- sequence_id(object)
  seq <- sequence_template(object)

  # create reverse complement of right primer
  rc_rev_primer <- Biostrings::reverseComplement(reverse_primer(object))

  # match right primer to seq
  match <- Biostrings::matchPattern(pattern = rc_rev_primer, subject = seq)
  match <- ranges(match)

  # get maximum product length if primer matches once, else abort with informative error message
  if (length(match) == 1) {
    max_product_size <- BiocGenerics::end(match)
    reverse_primer_start <- BiocGenerics::start(match)
    reverse_primer_length <- BiocGenerics::width(match)
  } else if (length(match) > 1) {
    stop("reverse primer matches the sequence multiple times for: ", seq_id, call. = FALSE)
  } else {
    stop("reverse primer doesn't match sequence template for: ", seq_id, call. = FALSE)
  }

  # ------------------------------------------------------------------------------------------------
  # STEP 2:
  # calculate excluded regions, that would generate products outside of the desired product size
  # range if primers would be designed there and specify them as excluded regions for Primer3
  # ------------------------------------------------------------------------------------------------

  # get the specified product size range and minimum primer region
  prod_size_range <- product_size_range(object)
  min_primer_region <- min_primer_region(object)

  # 2 excluded regions are designed, which span the sequence up- and downstream of the region that
  # yields primers with the desired product_size_range
  r1_start <- 0
  r1_length <- max_product_size - prod_size_range[2]

  r2_start <- max_product_size - prod_size_range[1] + 1
  r2_length <- length(seq) - r2_start

  # if the length of the first region is positive, both regions are valid and are created
  if (r1_length > 0) {

    # paste together in primer3 format
    sprintf("%s,%s %s,%s", r1_start, r1_length, r2_start, r2_length)

    # create r2 as only excluded region if the region resulting from r2 is >= the min_primer_region
  } else if (r2_start >= min_primer_region){

    sprintf("%s,%s", r2_start, r2_length)

    # else adapt r2 so that min_primer_size is respected (if possible)
  } else if (min_primer_region < reverse_primer_start) {

    warning("Desired product size range not possible! Product size will be shorter for: ",
            seq_id, call. = FALSE)

    r2_start <- min_primer_region
    r2_length <- length(seq) - r2_start
    sprintf("%s,%s", r2_start, r2_length)

  } else {
    stop("Can't respect min_primer_region, sequence too short for: ", seq_id, call. = FALSE)
  }
}
