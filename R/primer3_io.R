#' Design primers
#'
#' Design primers based on \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} objects.
#' Creates boulder-IO records, passes input to Primer3 and parses the output.
#'
#' @param object TsIO of TsIOList object for which primers should be designed.
#' @param thermo_params_path Optional path (character) to the \code{primer3_config} directory. Only
#'   required when using Primer3 < 2.5.0.
#' @param primer3_core Path (character) to the \code{primer3_core} executable. Usually this is
#'   inferred when loading/attaching the package.
#' @return A new \code{TsIO} or \code{TsIOList} object containing all Primer3 output.
#' @seealso \url{http://primer3.org/manual.html} for Primer3 manual and \code{\link[TAPseq]{TsIO}}
#'   for TsIO class objects.
#' @examples
#' library(TAPseq)
#'
#' # chromosome 11 region sequence templates
#' data("chr11_sequence_templates")
#'
#' # reverse primer used in all PCR reactions
#' reverse_primer <- "AAGCAGTGGTATCAACGCAGAGT"
#'
#' # create TsIOList object for the first two sequence templates
#' tapseq_io <- TAPseqInput(chr11_sequence_templates[1:2], reverse_primer = reverse_primer,
#'                          product_size_range = c(350, 500))
#'
#' # design primers
#' tapseq_io <- designPrimers(tapseq_io)
#'
#' # tapseq_io now contains the designed primers in the tapseq_primers slot
#' tapseq_primers(tapseq_io)
#' @export
setGeneric("designPrimers",
           function(object, thermo_params_path = NA,
                    primer3_core = getOption("TAPseq.primer3_core"))
             standardGeneric("designPrimers")
)

#' @describeIn designPrimers Design primers using Primer3 from a \code{TsIO} object
#' @export
setMethod("designPrimers", "TsIO", function(object, thermo_params_path, primer3_core) {

  # create boulder-IO record
  message("Converting to intput for Primer3...")
  io <- createIORecord(object, thermo_params_path = thermo_params_path)

  # design primers
  message("Running Primer3...")
  primer3_output <- system2(command = primer3_core, input = io, stdout = TRUE)

  # parse output and add it to input TsIO object(s)
  message("Processing output...")
  output <- parsePrimer3Output(object, primer3_output)

  # return output
  message("Done!")
  return(output)

})

#' @describeIn designPrimers Design primers using Primer3 from a \code{TsIOList} object
#' @export
setMethod("designPrimers", "TsIOList", function(object, thermo_params_path, primer3_core) {

  # create boulder-IO record
  message("Converting to input for Primer3...")
  io <- createIORecord(object, thermo_params_path = thermo_params_path)

  # design primers
  message("Running Primer3...")
  primer3_output <- system2(command = primer3_core, input = io, stdout = TRUE)

  # parse output and add it to input TsIO object(s)
  message("Processing output...")
  output <- parsePrimer3Output(object, primer3_output)

  # return output
  message("Done!")
  return(output)

})

#' Parse Primer3 Output
#'
#' Parse Primer3 output and add to input \code{\link[TAPseq]{TsIO}} or
#' \code{\link[TAPseq]{TsIOList}} object. This function is usually not used by the user, as
#' \code{\link[TAPseq]{designPrimers}} handels Primer3 output parsing.
#'
#' @param object The \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object used to
#'   design primers. No errors or warnings if this is another \code{TsIO} or \code{TsIOList} object!
#' @param primer3_output Character vector containing raw Primer3 output.
#' @return \code{TsIO} or \code{TsIOList} object with added Primer3 output
#' @examples
#' \dontrun{
#' library(TAPseq)
#'
#' # chromosome 11 region sequence templates
#' data("chr11_sequence_templates")
#'
#' # reverse primer used in all PCR reactions
#' reverse_primer <- "AAGCAGTGGTATCAACGCAGAGT"
#'
#' # create TsIOList object for the first two sequence templates
#' tapseq_io <- TAPseqInput(chr11_sequence_templates[1:2], reverse_primer = reverse_primer,
#'                          product_size_range = c(350, 500))
#'
#' # create boulder IO records
#' io_record <- createIORecord(tapseq_io)
#'
#' # design primers and store raw Primer3 output
#' primer3_core <- getOption("TAPseq.primer3_core")
#' primer3_output <- system2(command = primer3_core, input = io_record, stdout = TRUE)
#'
#' # parse output and add it to input TsIO object(s)
#' tapseq_io <- parsePrimer3Output(tapseq_io, primer3_output)
#' tapseq_primers(tapseq_io)
#' }
#' @keywords internal
parsePrimer3Output <- function(object, primer3_output) {

  # create a TsIOList object if a TsIO object is passed as input
  if (is(object, "TsIO")) {
    TsIO_input <- TRUE
    object <- TsIOList(object)
  }else{
    TsIO_input <- FALSE
  }

  # parse primer3 output into a list
  primer3_output <- parse_primer3_output(primer3_output)

  # parse every record in primer3_output to extract designed primers
  primers <- lapply(primer3_output, FUN = parse_primer3_record)

  # add designed primers to TsIO object(s)
  output <- mendoapply(FUN = `tapseq_primers<-`, object, primers)

  # infer pcr products
  output <- endoapply(output, FUN = infer_pcr_products)

  # return output
  if (TsIO_input == TRUE) {
    output[[1]]
  }else{
    output
  }

}

## HELPER FUNCTIONS ================================================================================

# Parse Primer3 output into a list of lists with named character vectors containing output per entry
parse_primer3_output <- function(x) {

  # get indices of final elements of each record (separators)
  seps <- which(x == "=")

  # drop record separators, as we now have their indices stored in seps
  xnosep <- x[-seps]

  # adjust seps so that they match the last entry of every record. we removed the original separator
  # element for each record, so the separator indices must be adjusted to match the final element of
  # each record: ei = si - i, where ei is the end of the record i, and si is the separator index of
  # record i from seps.
  rend <- seps - seq_along(seps)

  # split x into tags and values and transform into list of values with tags as names
  xsplit <- vapply(X = xnosep, FUN = strsplit, split = "=", FUN.VALUE = list(1),
                   USE.NAMES = FALSE)
  xlist <- vapply(X = xsplit, "[", 2, FUN.VALUE = character(1))
  names(xlist) <- tolower(vapply(X = xsplit, "[", 1, FUN.VALUE = character(1)))

  # assign each element in x_list to a record based on the indices of the last index of each record
  # inferred earlier
  record_indices <- rep.int(x = seq_along(rend), times = diff(c(0, rend)))

  # spit xlist into individual records
  split(xlist, f = record_indices)

}

# Parse one Primer3 output record and extract designed primers and return them as IRanges object
parse_primer3_record <- function(x) {

  # get sequence id and template
  seq_id <- x["sequence_id"]
  seq_template <- x["sequence_template"]

  # get entires in x that contain designed primers
  primers_idx <- grepl(names(x), pattern = "^primer_left_\\d+")
  primers <- x[primers_idx]

  # get primer id  for every primer
  primer_ids <- sub("(primer_left_\\d+).*", "\\1", names(primers))

  # split primers by primer id
  primers_split <- split(primers, f = primer_ids)

  # create IRange objects for every primer
  primer_ranges <- lapply(X = primers_split, FUN = function(x) {
    tryCatch({
      parse_primer(x, seq_template = seq_template, seq_id = seq_id)
    }, warning = function(w){
      message("Warning in parse_primer() for: ", seq_id)
      message(w, "")
      return(IRanges())
    }, error = function(e){
      message("Error in parse_primer() for: ", seq_id)
      message(e, "")
      return(IRanges())
    })
  })

  # convert to one IRanges object
  names(primer_ranges) <- NULL
  unlist(IRangesList(primer_ranges))

}

# parse a primer into a IRanges object containing the binding site in the sequence template and all
# primer3 stats as attached metadata columns
parse_primer <- function(primer, seq_template, seq_id){

  # get primer id
  primer_id <- unique(sub("(primer_left_\\d+).*", "\\1", names(primer)))

  # remove old primer coordinates entry, which has just the primer id as name
  primer <- primer[names(primer) != primer_id]

  # remove primer3 primer ids from names
  names(primer) <- sub(paste0(primer_id, "_"), "", names(primer))

  # prepare primer meta data
  primer_meta <- t(data.frame(primer))
  primer_meta <- data.frame(primer_meta, stringsAsFactors = FALSE)
  primer_meta <- utils::type.convert(primer_meta, as.is = TRUE)

  # extract primer sequence and transform to DNAString
  primer_seq <- DNAString(primer["sequence"])

  # find binding site of primer in sequence
  primer_site <- matchPattern(primer_seq, subject = seq_template)

  # create IRanges object with primer and it's binding site in the sequence template
  if (length(primer_site) == 1){

    # extract primer site coordinates and add meta data
    primer_range <- ranges(primer_site)
    mcols(primer_range) <- primer_meta

    # set name and return parsed primer
    names(primer_range) <- paste(seq_id, primer_id, sep=".")
    return(primer_range)

  }else if (length(primer_site) < 1) {
    stop("No matching sequence found for primer '", primer_id, "' in template!", call. = FALSE)
  }else{
    stop("More than 1 matching sequence found for primer '", primer_id, "' in template!",
         call. = FALSE)
  }

}

# infer pcr products based on designed primers. input must be a TsIO object with designed primers!
infer_pcr_products <- function(object) {

  # get sequence template and designed primers
  seq_template <- sequence_template(object)
  primers <- tapseq_primers(object)

  # find binding site of provided reverse primer
  rev_primer <- reverse_primer(object)
  rc_rev_primer <- reverseComplement(rev_primer)
  rev_primer_site <- matchPattern(rc_rev_primer, subject = seq_template)

  # create pcr products (if any primers are provided)
  pcr_products <- lapply(start(primers), FUN = subseq, x = seq_template, end = end(rev_primer_site))

  # convert to DNAstringSet and set primer ids as names
  pcr_products <- DNAStringSet(pcr_products)
  names(pcr_products) <- names(primers)

  # add pcr products to object and return object
  pcr_products(object) <- pcr_products
  return(object)

}
