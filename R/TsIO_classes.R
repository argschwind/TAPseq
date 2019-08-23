## Define the TsIO class and it's constructor and validation functions =============================

#' TsIO class
#'
#' TsIO objects store TASC-seq input and output for Primer3.
#'
#' The TsIO class is based on the Boulder IO records used by Primer3
#' (\href{http://primer3.sourceforge.net/primer3_manual.htm}{Primer3 manual}). These objects allow
#' storing of required input such as the sequence template and additional Primer3 arguments
#' relevant for TASC-seq primer design. Once Primer3 is used to design TASC-seq primers based on
#' a TsIO object, the output is added to the TsIO object.
#'
#' Use \code{TsIO()} to construct a new TsIO object from scratch.
#'
#' @param sequence_template A \code{\link[Biostrings]{DNAString}} object containing the sequence
#'   template for which primers should be designed.
#' @param sequence_id Name (\code{character}) of the sequence template e.g. gene name.
#' @param reverse_primer Reverse primer sequence used for all PCR reactions.
#' @param product_size_range Numerical vector of length 2 specifying the desired length of the
#'   resulting amplicons.
#' @param primer_num_return How many forward primers should be designed? (default: 5)
#' @param min_primer_region Minimum sequence length required for primer design. Relevant in case the
#'   a sequence template is too short to allow the specified \code{product_size_range}.
#' @param primer_opt_tm,primer_min_tm,primer_max_tm Optimal, minumum and maximum primer melting
#'   temperature.
#' @param tascseq_primers Slot where designed TASC-seq primers are stored. Not set by user.
#' @param pcr_products Slot where PCR products of primers are stored. Not set by user.
#' @param x A \code{TsIO} object.
#' @param value A valid value to assign to the chosen slot.
#' @return A \code{TsIO} object.
#' @seealso \url{http://primer3.org/manual.html} for Primer3 manual.
#'
setClass("TsIO",
  slots = c(
    sequence_id = "character",
    sequence_template = "DNAString",
    reverse_primer = "DNAString",
    product_size_range = "integer",
    primer_num_return = "integer",
    min_primer_region = "integer",
    primer_opt_tm = "integer",
    primer_min_tm = "integer",
    primer_max_tm = "integer",
    tascseq_primers = "IRanges",
    pcr_products = "DNAStringSet"
  ),
  prototype = list(
    sequence_id = NA_character_,
    sequence_template = new("DNAString"),
    reverse_primer = new("DNAString"),
    product_size_range = c(NA_integer_, NA_integer_),
    primer_num_return = NA_integer_,
    min_primer_region = NA_integer_,
    primer_opt_tm = NA_integer_,
    primer_min_tm = NA_integer_,
    primer_max_tm = NA_integer_,
    tascseq_primers = new("IRanges"),
    pcr_products = new("DNAStringSet")
  )
)

#' @rdname TsIO-class
#' @export
TsIO <- function(sequence_template, sequence_id = NA, reverse_primer, product_size_range,
                 primer_num_return = 5, min_primer_region = 100, primer_opt_tm = NA,
                 primer_min_tm = NA, primer_max_tm = NA) {

  # convert sequence_id to character
  sequence_id <- as.character(sequence_id)

  # convert sequence_template and reverse_primer into DNAString objects
  sequence_template <- as(sequence_template, "DNAString")
  reverse_primer <- as(reverse_primer, "DNAString")

  # make sure that product_size_range is a numeric vector of length 2
  if (!is(product_size_range, "numeric") | length(product_size_range) != 2) {
    stop("product_size_range must be a numeric vector of length 2!", call. = FALSE)
  }

  # convert parameters from numeric to integer
  product_size_range <- as(product_size_range, "integer")
  primer_num_return <- as(primer_num_return, "integer")
  min_primer_region <- as(min_primer_region, "integer")
  primer_opt_tm <- as(primer_opt_tm, "integer")
  primer_min_tm <- as(primer_min_tm, "integer")
  primer_max_tm <- as(primer_max_tm, "integer")

  # create new TsIO object
  new("TsIO", sequence_id = sequence_id, sequence_template = sequence_template,
      reverse_primer = reverse_primer,
      product_size_range = sort(product_size_range),
      primer_num_return = primer_num_return,
      min_primer_region = min_primer_region,
      primer_opt_tm = primer_opt_tm,
      primer_min_tm = primer_min_tm,
      primer_max_tm = primer_max_tm)

}

#' Validator function for TsIO objects
#'
#' @noRd
setValidity("TsIO", function(object) {
  if (length(object@sequence_id) != 1) {
    "sequence_id needs to be of length 1"
  }else{
    TRUE
  }
  if (length(object@product_size_range) != 2) {
    "product_size_range needs to be a numeric vector of length 2"
  }else{
    TRUE
  }
  if (length(object@primer_num_return) != 1) {
    "primer_num_return needs to be of length 1"
  }else{
    TRUE
  }
  if (length(object@min_primer_region) != 1) {
    "min_primer_region needs to be of length 1"
  }else{
    TRUE
  }
  if (length(object@primer_opt_tm) != 1) {
    "primer_opt_tm needs to be of length 1"
  }else{
    TRUE
  }
  if (length(object@primer_min_tm) != 1) {
    "primer_min_tm needs to be of length 1"
  }else{
    TRUE
  }
  if (length(object@primer_max_tm) != 1) {
    "primer_max_tm needs to be of length 1"
  }else{
    TRUE
  }
})


#' TsIOList class
#'
#' TsIOList class is a container to store multiple \code{\link[TASCseq]{TsIO}} objects. This enables
#' storing of Primer3 input and output for multiple target genes.
#'
#' @param ... Multiple TsIO objects from which a TsIOList object should be created.
#' @param x A \code{TsIOList} object.
#' @return A \code{TsIOList} object.
#' @seealso \link[TASCseq]{TsIO}
setClass("TsIOList",
         contains = "SimpleList",
         prototype = prototype(elementType = "TsIO"))

#' @rdname TsIOList-class
#' @export
TsIOList <- function(...) {

  # create list containing all passed TsIO objects
  objects <- list(...)

  # if already a list was passed to the function, this created a list of a list and needs to be
  # reverted
  if (length(objects) == 1L && is.list(objects[[1L]])) {
    objects <- objects[[1L]]
  }

  # create and return TsIOList object
  S4Vectors:::new_SimpleList_from_list(Class = "TsIOList", x = objects)

}


## GENERICS FOR TsIO and TsIOList OBJECTS ==========================================================

#' Accessors for TsIO objects
#'
#' A set of functions for getting/setting/modifying the data stored in \code{\link[TASCseq]{TsIO}}
#' or \code{\link[TASCseq]{TsIOList}} class objects.
#'
#' @param x A \code{TsIO} or \code{TsIOList} class object.
#' @param value A valid value to assign to the chosen slot.
#' @name accessors
NULL

#' @rdname accessors
#' @export
setGeneric("sequence_id", function(x) standardGeneric("sequence_id"))

#' @rdname accessors
#' @export
setGeneric("sequence_id<-", function(x, value) standardGeneric("sequence_id<-"))

#' @rdname accessors
#' @export
setGeneric("sequence_template", function(x) standardGeneric("sequence_template"))

#' @rdname accessors
#' @export
setGeneric("sequence_template<-", function(x, value) standardGeneric("sequence_template<-"))

#' @rdname accessors
#' @export
setGeneric("reverse_primer", function(x) standardGeneric("reverse_primer"))

#' @rdname accessors
#' @export
setGeneric("reverse_primer<-", function(x, value) standardGeneric("reverse_primer<-"))

#' @rdname accessors
#' @export
setGeneric("product_size_range", function(x) standardGeneric("product_size_range"))

#' @rdname accessors
#' @export
setGeneric("product_size_range<-", function(x, value) standardGeneric("product_size_range<-"))

#' @rdname accessors
#' @export
setGeneric("primer_num_return", function(x) standardGeneric("primer_num_return"))

#' @rdname accessors
#' @export
setGeneric("primer_num_return<-", function(x, value) standardGeneric("primer_num_return<-"))

#' @rdname accessors
#' @export
setGeneric("min_primer_region", function(x) standardGeneric("min_primer_region"))

#' @rdname accessors
#' @export
setGeneric("min_primer_region<-", function(x, value) standardGeneric("min_primer_region<-"))

#' @rdname accessors
#' @export
setGeneric("primer_opt_tm", function(x) standardGeneric("primer_opt_tm"))

#' @rdname accessors
#' @export
setGeneric("primer_opt_tm<-", function(x, value) standardGeneric("primer_opt_tm<-"))

#' @rdname accessors
#' @export
setGeneric("primer_min_tm", function(x) standardGeneric("primer_min_tm"))

#' @rdname accessors
#' @export
setGeneric("primer_min_tm<-", function(x, value) standardGeneric("primer_min_tm<-"))

#' @rdname accessors
#' @export
setGeneric("primer_max_tm", function(x) standardGeneric("primer_max_tm"))

#' @rdname accessors
#' @export
setGeneric("primer_max_tm<-", function(x, value) standardGeneric("primer_max_tm<-"))

#' @rdname accessors
#' @export
setGeneric("tascseq_primers", function(x) standardGeneric("tascseq_primers"))

#' @keywords internal
setGeneric("tascseq_primers<-", function(x, value) standardGeneric("tascseq_primers<-"))

#' @rdname accessors
#' @export
setGeneric("pcr_products", function(x) standardGeneric("pcr_products"))

#' @keywords internal
setGeneric("pcr_products<-", function(x, value) standardGeneric("pcr_products<-"))


## METHODS =========================================================================================

# TsIO objects -------------------------------------------------------------------------------------

#' @describeIn TsIO Get sequence id
#' @export
setMethod("sequence_id", "TsIO", function(x) x@sequence_id)

#' @describeIn TsIO Set sequence_id
#' @export
setMethod("sequence_id<-", "TsIO", function(x, value) {
  x@sequence_id <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get sequence_template
#' @export
setMethod("sequence_template", "TsIO", function(x) x@sequence_template)

#' @describeIn TsIO Set sequence_template
#' @export
setMethod("sequence_template<-", "TsIO", function(x, value) {
  value <- as(value, "DNAString")
  x@sequence_template <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get reverse_primer
#' @export
setMethod("reverse_primer", "TsIO", function(x) x@reverse_primer)

#' @describeIn TsIO Set reverse_primer
#' @export
setMethod("reverse_primer<-", "TsIO", function(x, value) {
  value <- as(value, "DNAString")
  x@reverse_primer <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get product_size_range
#' @export
setMethod("product_size_range", "TsIO", function(x) x@product_size_range)

#' @describeIn TsIO Set product_size_range
#' @export
setMethod("product_size_range<-", "TsIO", function(x, value) {
  value <- sort(as(value, "integer"))
  x@product_size_range <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_num_return
#' @export
setMethod("primer_num_return", "TsIO", function(x) x@primer_num_return)

#' @describeIn TsIO Set product_size_range
#' @export
setMethod("primer_num_return<-", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_num_return <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get min_primer_region
#' @export
setMethod("min_primer_region", "TsIO", function(x) x@min_primer_region)

#' @describeIn TsIO Set product_size_range
#' @export
setMethod("min_primer_region<-", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@min_primer_region <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_opt_tm
#' @export
setMethod("primer_opt_tm", "TsIO", function(x) x@primer_opt_tm)

#' @describeIn TsIO Set product_size_range
#' @export
setMethod("primer_opt_tm<-", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_opt_tm <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_min_tm
#' @export
setMethod("primer_min_tm", "TsIO", function(x) x@primer_min_tm)

#' @describeIn TsIO Set product_size_range
#' @export
setMethod("primer_min_tm<-", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_min_tm <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_max_tm
#' @export
setMethod("primer_max_tm", "TsIO", function(x) x@primer_max_tm)

#' @describeIn TsIO Set product_size_range
#' @export
setMethod("primer_max_tm<-", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_max_tm <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get tascseq_primers
#' @export
setMethod("tascseq_primers", "TsIO", function(x) x@tascseq_primers)

#' @keywords internal
setMethod("tascseq_primers<-", "TsIO", function(x, value) {
  x@tascseq_primers <- value
  x
})

#' @describeIn TsIO Get pcr_products
#' @export
setMethod("pcr_products", "TsIO", function(x) x@pcr_products)

#' @keywords internal
setMethod("pcr_products<-", "TsIO", function(x, value) {
  x@pcr_products <- value
  x
})


# TsIOList objects ---------------------------------------------------------------------------------

#' @describeIn TsIOList Get sequence_templates
#' @export
setMethod("sequence_template", "TsIOList", function(x) {
  seq_templates <- lapply(x, FUN = sequence_template)
  Biostrings::DNAStringSet(seq_templates)
})

#' @describeIn TsIOList Get tascseq_primers
#' @export
setMethod("tascseq_primers", "TsIOList", function(x) {
  primers <- lapply(x, FUN = tascseq_primers)
  names(primers) <- NULL
  BiocGenerics::unlist(IRanges::IRangesList(primers))
})

#' @describeIn TsIOList Get pcr_products
#' @export
setMethod("pcr_products", "TsIOList", function(x) {
  pcr_prod <- lapply(x, FUN = pcr_products)
  names(pcr_prod) <- NULL
  BiocGenerics::unlist(Biostrings::DNAStringSetList(pcr_prod))
})


## Miscellaneous helper functions ==================================================================

#' Show method for TsIO objects
#'
#' @keywords internal
setMethod("show", "TsIO", function(object) {
  cat(is(object)[[1]], " instance", "\n",
      "  ", length(object@sequence_template), " bp sequence template", "\n",
      "  seqID: ", object@sequence_id, "\n",
      "  right primer: ", as.character(object@reverse_primer), "\n",
      "  specified product size range: ", object@product_size_range[1], "-",
           object@product_size_range[2], " basepairs",
      sep = "")
})
