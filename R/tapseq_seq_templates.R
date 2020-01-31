#' Create TAPseq sequence templates
#'
#' Methods to create sequence templates for TAP-seq primer design. This extracts the transcript
#' sequence from the provided genome and adds the reverse complement of the "Beads-oligo-dT"
#' sequence to the 3' end.
#'
#' @param transcripts A \code{\link[GenomicRanges]{GRanges}} or
#'   \code{\link[GenomicRanges]{GRangesList}} object containing exons of transcripts for which
#'   sequence templates should be created. All exons in a \code{\link[GenomicRanges]{GRanges}}
#'   object are assumed to belong to the same transcript. Multiple transcripts can be provided in a
#'   \code{\link[GenomicRanges]{GRangesList}} object.
#' @param  genome A \code{\link[BSgenome]{BSgenome}} or \code{\link[Biostrings]{DNAStringSet}}
#'   object containing chromosome sequences which should be used to extract transcript sequences.
#' @param beads_oligo (character) The Beads-oligo-dT sequence used for droplet sequencing. If
#'   nothing is specified (\code{beads_oligo = NA}), the 10x V3 Beads-oligo-dT sequence is used.
#'   Can be changed if primers are for instance designed for Drop-seq. Any barcode bases need to be
#'   replaced by \code{N}.
#' @return A \code{\link[Biostrings]{DNAString}} or \code{\link[Biostrings]{DNAStringSet}} object
#'   containing the sequence template(s).
#' @examples
#' \dontrun{
#' library(TAPseq)
#' library(Biostrings)
#' library(BSgenome)
#'
#' # truncated transcripts for chr11 target genes
#' data("chr11_truncated_txs")
#'
#' # human genome (hg38) BSgenome object (needs to be istalled separately from Bioconductor)
#' # the genome sequence could also be loaded from a .fasta file with
#' hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
#'
#' # create sequence templates for all target transcripts on chr11
#' seq_templates <- TAPseqSeqTemplates(chr11_truncated_txs, genome = hg38)
#'
#' # create sequence templates using the Drop-seq Beads-oligo-dT sequence
#' ds_oligo <- "TTTTTTTAAGCAGTGGTATCAACGCAGAGTACNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
#' seq_templates <- TAPseqSeqTemplates(chr11_truncated_txs, genome = hg38, beads_oligo = ds_oligo)
#'}
#' @export
setGeneric("TAPseqSeqTemplates",
           function(transcripts, genome, beads_oligo = NA) standardGeneric("TAPseqSeqTemplates")
)

#' @describeIn TAPseqSeqTemplates Create sequence templates from \code{GRanges} input
#' @export
setMethod("TAPseqSeqTemplates", "GRanges", function(transcripts, genome, beads_oligo) {

  # extract transcript sequences
  tx_seqs <- getTxsSeq(transcripts, genome = genome)

  # create reverse complement of Beads-oligo-dT to add to 3' end of transcripts
  if (is.null(beads_oligo)) {
    stop("beads_oligo cannot be NULL", call. = FALSE)
  }else if (is.na(beads_oligo)) {
    beads_oligo <- reverseComplement(DNAString(get_beads_oligo()))
  }else{
    beads_oligo <- reverseComplement(DNAString(beads_oligo))
  }

  # add bead_seq to transcript sequences to create sequence templates
  seq_templates <- Biostrings::xscat(tx_seqs, beads_oligo)
  names(seq_templates) <- names(tx_seqs)

  return(seq_templates)

}
)

#' @describeIn TAPseqSeqTemplates Create sequence templates from \code{GRangesList} input
#' @export
setMethod("TAPseqSeqTemplates", "GRangesList", function(transcripts, genome, beads_oligo) {

  # extract transcript sequences
  tx_seqs <- getTxsSeq(transcripts, genome = genome)

  # create reverse complement of Beads-oligo-dT to add to 3' end of transcripts
  if (is.null(beads_oligo)) {
    stop("beads_oligo cannot be NULL", call. = FALSE)
  }else if (is.na(beads_oligo)) {
    beads_oligo <- reverseComplement(DNAString(get_beads_oligo()))
  }else{
    beads_oligo <- reverseComplement(DNAString(beads_oligo))
  }

  # add bead_seq to transcript sequences to create sequence templates
  seq_templates <- Biostrings::xscat(tx_seqs, beads_oligo)
  names(seq_templates) <- names(tx_seqs)

  return(seq_templates)

}
)

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
