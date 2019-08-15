#' Create TASCseq input from sequence templates
#'
#' This function creates input for TASC-seq primer design from a DNAStringSet containing the
#' sequence templates.
#'
#' @param sequence_templates A \code{\link{DNAStringSet}} object containing all sequence templates.
#' @param reverse_primer Reverse primer sequence used for all PCR reactions.
#' @param product_size_range Numerical vector of length 2 specifying the desired length of the
#'   resulting amplicons.
#' @param primer_num_return How many forward primers should be designed? (default: 5)
#' @param min_primer_region Minimum sequence length required for primer design. Relevant in case the
#'   a sequence template is too short to allow the specified \code{product_size_range}.
#' @param primer_opt_tm,primer_min_tm,primer_max_tm Optimal, minumum and maximum primer melting
#'   temperature.
#' @return \code{\link[TASCseq]{TsIOList}} object.
#' @export
TASCseqInput <- function(sequence_templates, reverse_primer, product_size_range,
                         primer_num_return = 5, min_primer_region = 100, primer_opt_tm = NA,
                         primer_min_tm = NA, primer_max_tm = NA) {

  # make sure that sequence template has the right format
  if (!is(sequence_templates, "DNAStringSet")) {
    stop("sequence_templates must be a DNAStringSet object")
  }

  # create TsIO objects for every sequence template
  io <- lapply(sequence_templates, FUN = TsIO, reverse_primer = reverse_primer,
               product_size_range = product_size_range, primer_num_return = primer_num_return,
               min_primer_region = min_primer_region, primer_opt_tm = primer_opt_tm,
               primer_min_tm = primer_min_tm, primer_max_tm = primer_max_tm)

  # set sequence ids, which were lost when creating io list
  io <- mapply(FUN = `sequence_id<-`, io, names(io), SIMPLIFY = FALSE, USE.NAMES = TRUE)

  # convert to TsIOList object
  TsIOList(io)

}

#' Create boulder IO record
#'
#' Takes a \code{\link[TASCseq]{TsIO}} or \code{\link[TASCseq]{TsIOList}} object and converts it
#' into a boulder IO record for Primer3. Essentially it converts it into a list of character vectors
#' that each contain the tag and the value in the form: "TAG=VALUE". More on this format can be
#' found in the \href{http://primer3.org/manual.html}{Primer3 manual}.
#'
#' This function is usually not needed by the user, because functions such as
#' \code{\link[TASCseq]{designPrimers}} handle IO record generation. However, this function can for
#' instance be useful to generate IO records, write them to a file and pass them to Primer3 in the
#' conventional way.
#'
#' @param object TsIO of TsIOList object for which a Primer3 boulder IO record should be created.
#' @param thermo_params_path Path (character) to the \code{primer3_config} directory. Default set to
#'   the same directory where \code{primer3_core} executable is found.
#' @return A character vectors containing the individual lines of the IO record.
#' @seealso \url{http://primer3.org/manual.html} for Primer3 manual.
#' @export
setGeneric("createIORecord",
           function(object, thermo_params_path = getOption("TASCseq.thermodynamic_params_path"))
             standardGeneric("createIORecord")
)

#' @describeIn createIORecord Create IO record from \code{TsIO} objects.
#' @export
setMethod("createIORecord", "TsIO", function(object, thermo_params_path) {

  # get slot names of TsIO class
  slot_names <- slotNames(object)

  # exclude any slots intended for output
  output_slots <- c("tascseq_primers", "pcr_products", "blast_off_targets")
  slot_names <- slot_names[!slot_names %in% output_slots]

  # exclude min_primer_region and product_size_range, since it's are processed separately
  slot_names <- slot_names[!slot_names %in% c("min_primer_region", "product_size_range")]

  # get data of all slots as one character vector
  io <- vapply(slot_names, FUN = function(s, x) {
    as.character(slot(x, s))
  }, x = object, FUN.VALUE = character(1))

  # remove any NA elements, because they should not be part of the primer3 input
  io <- io[!is.na(io)]

  # rename reverse_primer to Primer3 argument name
  rev_index <- which(names(io) == "reverse_primer")
  names(io)[rev_index] <- "sequence_primer_revcomp"

  # add Primer3 parameters required to trigger design of only forward primers
  io["primer_pick_left_primer"] <- 1
  io["primer_pick_right_primer"] <- 0

  # calculate and add excluded regions
  excluded_regions <- create_excluded_regions(object)
  io["sequence_excluded_region"] <- excluded_regions

  # transform to vector where each element contains one input line ("tag=value" format)
  io <- paste0(toupper(names(io)), "=", io)

  # add thermodynamic parameters path to record
  io <- c(io, paste0("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=", thermo_params_path))

  # add "=" record separator at the end of the output list
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
  }else if (length(match) > 1) {
    stop("reverse primer matches the sequence multiple times for: ", seq_id, call. = FALSE)
  }else{
    stop("reverse primer doesn't match sequence template for: ", seq_id)
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
  }else if (r2_start >= min_primer_region){

    sprintf("%s,%s", r2_start, r2_length)

    # else adapt r2 so that min_primer_size is respected (if possible)
  }else if (min_primer_region < reverse_primer_start) {

    warning("Desired product size range not possible! Product size will be shorter for: ",
            seq_id, call. = FALSE)

    r2_start <- min_primer_region
    r2_length <- length(seq) - r2_start
    sprintf("%s,%s", r2_start, r2_length)

  }else{
    stop("Can't respect min_primer_region, sequence too short for: ", seq_id, call. = FALSE)
  }
}
