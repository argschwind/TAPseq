## Define the TsIO class and it's constructor and validation functions =============================

#' TsIO class
#'
#' TsIO objects store TAP-seq Primer3 input and output.
#'
#' The TsIO class is based on the Boulder IO records used by Primer3
#' (\href{http://primer3.sourceforge.net/primer3_manual.htm}{Primer3 manual}). These objects are
#' used to store the sequence templates and parameters needed for TAP-seq primer design. Primers
#' designed  with Primer3 are also stored in the same TsIO objects.
#'
#' Use \code{TsIO()} to construct a new TsIO object from scratch.
#'
#' @param target_sequence A \code{\link[Biostrings:DNAString-class]{DNAString}} or \code{character}
#'   object containing the target sequence for which primers should be designed. Usually a
#'   transcript sequence.
#' @param beads_oligo Beads-oligo-dT sequence for the used droplet sequencing protocol (10x,
#'   Drop-seq).
#' @param reverse_primer Reverse primer sequence used for all PCR reactions.
#' @param product_size_range Numerical vector of length 2 specifying the desired length of the
#'   resulting amplicons.
#' @param sequence_id Name (\code{character}) of the target sequence, e.g. the gene name. It's
#'  adviced to use meaningful sequence ids to savely assign designed primers to their targets.
#' @param target_annot (optional) A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with
#'   transcript annotation in case the target is a transcript of a gene locus. If provided, it
#'   should contain all exons of the targeted transcript.
#' @param primer_num_return How many forward primers should be designed? (default: 5)
#' @param min_primer_region Minimum sequence length required for primer design. Mostly relevant in
#'   case the sequence template is too short to allow the specified \code{product_size_range}.
#' @param primer_opt_tm,primer_min_tm,primer_max_tm Optimal, minumum and maximum primer melting
#'   temperature.
#' @param tapseq_primers Slot where designed TAP-seq primers are stored. Not set by user.
#' @param pcr_products Slot where PCR products of primers are stored. Not set by user.
#' @param x A \code{TsIO} object.
#' @param value A valid value to assign to the chosen slot.
#' @return A \code{TsIO} object.
#' @seealso \url{http://primer3.org/manual.html} for Primer3 manual.
#' @examples
#' # get example transcript sequence
#' data("chr11_truncated_txs_seq")
#' tx_seq <- chr11_truncated_txs_seq[[1]]
#' tx_id <- names(chr11_truncated_txs_seq)[1]
#'
#' # 10x beads-oligo-dt sequence
#' beads_oligo <- "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
#'
#' # reverse primer used in all PCR reactions
#' reverse_primer <- "CTACACGACGCTCTTCCGATCT"
#'
#' # create TsIO object
#' obj <- TsIO(sequence_id = tx_id, target_sequence = tx_seq, beads_oligo = beads_oligo,
#'             reverse_primer = reverse_primer, product_size_range = c(350, 500))
#'
#' # slot values can be accessed using accessor functions
#' sequence_id(obj)
#' sequence_id(obj) <- "Gene1"
#' sequence_id(obj)
#'
#' # the sequence template (target sequence + reverse complement of beads-oligo-dt) for primer
#' # design can be viewed as well
#' sequence_template(obj)
setClass("TsIO",
  slots = c(
    sequence_id = "character",
    target_sequence = "DNAString",
    beads_oligo = "DNAString",
    reverse_primer = "DNAString",
    target_annot = "GRanges",
    product_size_range = "integer",
    primer_num_return = "integer",
    min_primer_region = "integer",
    primer_opt_tm = "integer",
    primer_min_tm = "integer",
    primer_max_tm = "integer",
    tapseq_primers = "IRanges",
    pcr_products = "DNAStringSet"
  ),
  prototype = list(
    sequence_id = NA_character_,
    target_sequence = new("DNAString"),
    beads_oligo = new("DNAString"),
    reverse_primer = new("DNAString"),
    target_annot = new("GRanges"),
    product_size_range = c(NA_integer_, NA_integer_),
    primer_num_return = NA_integer_,
    min_primer_region = NA_integer_,
    primer_opt_tm = NA_integer_,
    primer_min_tm = NA_integer_,
    primer_max_tm = NA_integer_,
    tapseq_primers = new("IRanges"),
    pcr_products = new("DNAStringSet")
  )
)

#' @rdname TsIO-class
#' @export
TsIO <- function(sequence_id, target_sequence, beads_oligo, reverse_primer, product_size_range,
                 target_annot = NULL, primer_num_return = 5, min_primer_region = 100,
                 primer_opt_tm = NA, primer_min_tm = NA, primer_max_tm = NA) {

  # convert sequence_id to character
  sequence_id <- as.character(sequence_id)

  # convert target_sequence, beads_oligo and reverse_primer into DNAString objects
  target_sequence <- as(target_sequence, "DNAString")
  beads_oligo <- as(beads_oligo, "DNAString")
  reverse_primer <- as(reverse_primer, "DNAString")

  # make sure that target_annot is a GRanges object
  target_annot <- GRanges(target_annot)

  # convert other parameters to integers
  product_size_range <- sort(as(product_size_range, "integer"), na.last = TRUE)
  primer_num_return <- as(primer_num_return, "integer")
  min_primer_region <- as(min_primer_region, "integer")
  primer_opt_tm <- as(primer_opt_tm, "integer")
  primer_min_tm <- as(primer_min_tm, "integer")
  primer_max_tm <- as(primer_max_tm, "integer")

  # create new TsIO object
  new("TsIO",
      sequence_id = sequence_id,
      target_sequence = target_sequence,
      beads_oligo = beads_oligo,
      reverse_primer = reverse_primer,
      target_annot = target_annot,
      product_size_range = product_size_range,
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

  # get target annot strand and chromosome, and merge any overlapping exons
  annot <- target_annot(object)
  seqname <- unique(seqnames(annot))
  strand <- unique(strand(annot))
  annot_merge <- reduce(annot)

  # check object slots
  err <- character()
  if (length(seqname) > 1)
    err <- c(err, "multiple seqnames in target_annot")
  if (length(strand) > 1)
    err <- c(err, "conflicting strand information in target_annot")
  if (length(annot) != length(annot_merge))
    err <- c(err, "overlapping exons found in target_annot")
  if (!sum(width(annot_merge)) %in% c(0, length(object@target_sequence)))
    err <- c(err, "exons in target_annot are incompatible with target_sequence")
  if (length(object@sequence_id) != 1)
    err <- c(err, "sequence_id needs to be of length 1")
  if (length(object@product_size_range) != 2)
    err <- c(err, "product_size_range needs to be an integer vector of length 2")
  if (length(object@primer_num_return) != 1)
    err <- c(err, "primer_num_return needs to be of length 1")
  if (length(object@min_primer_region) != 1)
    err <- c(err, "min_primer_region needs to be of length 1")
  if (length(object@primer_opt_tm) != 1)
    err <- c(err, "primer_opt_tm needs to be of length 1")
  if (length(object@primer_min_tm) != 1)
    err <- c(err, "primer_min_tm needs to be of length 1")
  if (length(object@primer_max_tm) != 1)
    err <- c(err, "primer_max_tm needs to be of length 1")

  # check that specified product_size_range and min_primer_range are compatible
  if (length(object@product_size_range) == 2) {
    if (!any(is.na(c(object@product_size_range, object@min_primer_region)))) {
      if (diff(object@product_size_range) < object@min_primer_region) {
        err <- c(err, "product_size_range too narrow to allow min_primer_range")
      }
    }
  }

  # return all slot errors if any are found
  if (length(err) > 0) err else TRUE

})


#' TsIOList class
#'
#' TsIOList class is a container to store multiple \code{\link[TAPseq:TsIO-class]{TsIO}} objects.
#' This enables storing of Primer3 input and output for multiple target genes.
#'
#' @param ... Multiple TsIO objects from which a TsIOList object should be created.
#' @param x A \code{TsIOList} object.
#' @return A \code{TsIOList} object.
#' @seealso \link[TAPseq:TsIO-class]{TsIO}
#' @examples
#' # get example transcript sequences
#' data("chr11_truncated_txs_seq")
#' txs_seqs <- chr11_truncated_txs_seq[1:2]
#' txs_ids <- names(txs_seqs)
#'
#' # 10x beads-oligo-dt sequence
#' beads_oligo <- "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
#'
#' # reverse primer used in all PCR reactions
#' reverse_primer <- "CTACACGACGCTCTTCCGATCT"
#'
#' # create TsIO objects
#' tsio1 <- TsIO(sequence_id = txs_ids[1], target_sequence = txs_seqs[[1]],
#'               beads_oligo = beads_oligo, reverse_primer = reverse_primer,
#'               product_size_range = c(350, 500))
#'
#' tsio2 <- TsIO(sequence_id = txs_ids[2], target_sequence = txs_seqs[[2]],
#'               beads_oligo = beads_oligo, reverse_primer = reverse_primer,
#'               product_size_range = c(350, 500))
#'
#' # create TsIOList object
#' obj <- TsIOList(tsio1 = tsio1, tsio2 = tsio2)
#'
#' # it's noteworthy to mention that when creating a TsIOList from a DNAStringSet of target
#' # sequences, it's easier to use TAPseqInput()
#' ?TAPseqInput
#'
#' # as with TsIO objects, some values can be accessed using accessor functions
#' sequence_template(obj)
setClass("TsIOList", contains = "SimpleList", prototype = "TsIO")

#' @rdname TsIOList-class
#' @export
TsIOList <- function(...) {

  # create list containing all passed TsIO objects
  objects <- list(...)

  # if already a list or TsIOList was passed to the function, this created a list of a list and
  # needs to be reverted
  if (length(objects) == 1L) {
    tmp <- objects[[1L]]
    if (is.list(tmp) | is(tmp, "List")) {
      objects <- tmp
    }
  }

  # create and return TsIOList object
  new("TsIOList", objects, elementType = "TsIO")

}


## GENERICS FOR TsIO and TsIOList OBJECTS ==========================================================

#' Accessors for TsIO objects
#'
#' A set of functions for getting/setting/modifying the data stored in
#' \code{\link[TAPseq:TsIO-class]{TsIO}} or \code{\link[TAPseq:TsIOList-class]{TsIOList}} class
#' objects.
#'
#' @param x A \code{TsIO} or \code{TsIOList} class object.
#' @param value A valid value to assign to the chosen slot.
#' @return Returns the stored value(s) of a slot, or sets a new value
#' @examples
#' # chr11 primers example data
#' data("chr11_primers")
#'
#' # slot values of TsIO objects can be accessed using accessor functions
#' tsio <- chr11_primers[[1]]
#' sequence_id(tsio)
#' sequence_id(tsio) <- "Gene1"
#' sequence_id(tsio)
#'
#' # some slots can only be obtained, but not set as filling these is part of the TAPseq workflow
#' tapseq_primers(tsio)
#' pcr_products(tsio)
#'
#' # sequence templates can be created
#' sequence_template(tsio)
#'
#' # values of TsIOList object slots can be extracted as well, but not set
#' tsio_list <- chr11_primers[1:2]
#' sequence_id(tsio_list)
#' target_sequence(tsio_list)
#' target_annot(tsio_list)
#' tapseq_primers(tsio_list)
#' pcr_products(tsio_list)
#' sequence_template(tsio_list)
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
setGeneric("target_sequence", function(x) standardGeneric("target_sequence"))

#' @rdname accessors
#' @export
setGeneric("target_sequence<-", function(x, value) standardGeneric("target_sequence<-"))

#' @rdname accessors
#' @export
setGeneric("beads_oligo", function(x) standardGeneric("beads_oligo"))

#' @rdname accessors
#' @export
setGeneric("beads_oligo<-", function(x,value) standardGeneric("beads_oligo<-"))

#' @rdname accessors
#' @export
setGeneric("reverse_primer", function(x) standardGeneric("reverse_primer"))

#' @rdname accessors
#' @export
setGeneric("reverse_primer<-", function(x, value) standardGeneric("reverse_primer<-"))

#' @rdname accessors
#' @export
setGeneric("target_annot", function(x) standardGeneric("target_annot"))

#' @rdname accessors
#' @export
setGeneric("target_annot<-", function(x, value) standardGeneric("target_annot<-"))

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
setGeneric("sequence_template", function(x) standardGeneric("sequence_template"))

#' @rdname accessors
#' @export
setGeneric("tapseq_primers", function(x) standardGeneric("tapseq_primers"))

#' @keywords internal
setGeneric("tapseq_primers<-", function(x, value) standardGeneric("tapseq_primers<-"))

#' @rdname accessors
#' @export
setGeneric("pcr_products", function(x) standardGeneric("pcr_products"))

#' @keywords internal
setGeneric("pcr_products<-", function(x, value) standardGeneric("pcr_products<-"))


## METHODS =========================================================================================

# TsIO objects -------------------------------------------------------------------------------------

#' @describeIn TsIO Get sequence_id
#' @export
setMethod("sequence_id", "TsIO", function(x) x@sequence_id)

#' @describeIn TsIO Set sequence_id
#' @export
setReplaceMethod("sequence_id", "TsIO", function(x, value) {
  x@sequence_id <- as.character(value)
  validObject(x)
  x
})

#' @describeIn TsIO Get target_sequence
#' @export
setMethod("target_sequence", "TsIO", function(x) x@target_sequence)

#' @describeIn TsIO Set target_sequence
#' @export
setReplaceMethod("target_sequence", "TsIO", function(x, value) {
  value <- as(value, "DNAString")
  x@target_sequence <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get beads_oligo
#' @export
setMethod("beads_oligo", "TsIO", function(x) x@beads_oligo)

#' @describeIn TsIO Set beads_oligo
#' @export
setReplaceMethod("beads_oligo", "TsIO", function(x, value) {
  value <- as(value, "DNAString")
  x@beads_oligo <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get reverse_primer
#' @export
setMethod("reverse_primer", "TsIO", function(x) x@reverse_primer)

#' @describeIn TsIO Set reverse_primer
#' @export
setReplaceMethod("reverse_primer", "TsIO", function(x, value) {
  value <- as(value, "DNAString")
  x@reverse_primer <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get target_annot
#' @export
setMethod("target_annot", "TsIO", function(x) x@target_annot)

#' @describeIn TsIO Set target_annot
#' @export
setReplaceMethod("target_annot", "TsIO", function(x, value) {
  x@target_annot <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get product_size_range
#' @export
setMethod("product_size_range", "TsIO", function(x) x@product_size_range)

#' @describeIn TsIO Set product_size_range
#' @export
setReplaceMethod("product_size_range", "TsIO", function(x, value) {
  value <- sort(as(value, "integer"), na.last = TRUE)
  x@product_size_range <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_num_return
#' @export
setMethod("primer_num_return", "TsIO", function(x) x@primer_num_return)

#' @describeIn TsIO Set primer_num_return
#' @export
setReplaceMethod("primer_num_return", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_num_return <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get min_primer_region
#' @export
setMethod("min_primer_region", "TsIO", function(x) x@min_primer_region)

#' @describeIn TsIO Set min_primer_region
#' @export
setReplaceMethod("min_primer_region", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@min_primer_region <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_opt_tm
#' @export
setMethod("primer_opt_tm", "TsIO", function(x) x@primer_opt_tm)

#' @describeIn TsIO Set primer_opt_tm
#' @export
setReplaceMethod("primer_opt_tm", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_opt_tm <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_min_tm
#' @export
setMethod("primer_min_tm", "TsIO", function(x) x@primer_min_tm)

#' @describeIn TsIO Set primer_min_tm
#' @export
setReplaceMethod("primer_min_tm", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_min_tm <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get primer_max_tm
#' @export
setMethod("primer_max_tm", "TsIO", function(x) x@primer_max_tm)

#' @describeIn TsIO Set primer_max_tm
#' @export
setReplaceMethod("primer_max_tm", "TsIO", function(x, value) {
  value <- as(value, "integer")
  x@primer_max_tm <- value
  validObject(x)
  x
})

#' @describeIn TsIO Create sequence_template
#' @export
setMethod("sequence_template", "TsIO", function(x) {
  target_sequence <- x@target_sequence
  beads_oligo_revcomp <- reverseComplement(DNAString(x@beads_oligo))
  Biostrings::xscat(target_sequence, beads_oligo_revcomp)
})

#' @describeIn TsIO Get tapseq_primers
#' @export
setMethod("tapseq_primers", "TsIO", function(x) x@tapseq_primers)

#' @keywords internal
setReplaceMethod("tapseq_primers", "TsIO", function(x, value) {
  x@tapseq_primers <- value
  validObject(x)
  x
})

#' @describeIn TsIO Get pcr_products
#' @export
setMethod("pcr_products", "TsIO", function(x) x@pcr_products)

#' @keywords internal
setReplaceMethod("pcr_products", "TsIO", function(x, value) {
  x@pcr_products <- value
  validObject(x)
  x
})


# TsIOList objects ---------------------------------------------------------------------------------

#' @describeIn TsIOList Get sequence_id
#' @export
setMethod("sequence_id", "TsIOList", function(x) {
  vapply(x, FUN = sequence_id, FUN.VALUE = character(1))
})

#' @describeIn TsIOList Get target_sequence
#' @export
setMethod("target_sequence", "TsIOList", function(x) {
  seq_templates <- lapply(x, FUN = target_sequence)
  Biostrings::DNAStringSet(seq_templates)
})

#' @describeIn TsIOList Create sequence_template
#' @export
setMethod("sequence_template", "TsIOList", function(x) {
  seq_templates <- lapply(x, FUN = sequence_template)
  Biostrings::DNAStringSet(seq_templates)
})

#' @describeIn TsIOList Get target_annot
#' @export
setMethod("target_annot", "TsIOList", function(x) {
  annots <- lapply(x, FUN = target_annot)
  GenomicRanges::GRangesList(annots)
})

#' @describeIn TsIOList Get tapseq_primers
#' @export
setMethod("tapseq_primers", "TsIOList", function(x) {
  primers <- lapply(x, FUN = tapseq_primers)
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


## MISCELLANEOUS FUNCTIONS =========================================================================

# show method for TsIO objects to print information about the object to console
setMethod("show", "TsIO", function(object) {

  beads_oligo <- beads_oligo(object)
  if (length(beads_oligo) == 0) beads_oligo <- NA

  reverse_primer <- reverse_primer(object)
  if (length(reverse_primer) == 0) reverse_primer <- NA

  prod_size <- product_size_range(object)

  primers <- length(tapseq_primers(object))
  if (primers == 0) primers <- "None"

  cat(is(object)[[1]], " instance", "\n",
      "  ", length(sequence_template(object)), " bp sequence template", "\n",
      "  Sequence ID: ", sequence_id(object), "\n",
      "  Beads oligo: ", as.character(beads_oligo), "\n",
      "  Right primer: ", as.character(reverse_primer), "\n",
      "  Specified product size range: ", prod_size[1], "-", prod_size[2], " basepairs", "\n",
      "  Designed left primers: ", primers, "\n",
      sep = "")
})
