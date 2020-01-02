#' Get transcript sequences
#'
#' Extract the DNA sequences of all exons of transcript models and concatenate
#' to one sequence per transcript. This is basically a wrapper for
#' \code{\link[GenomicFeatures]{extractTranscriptSeqs}}, which makes sure that the exons are
#' correctly sorted according to their position in the transcript (3' to 5').
#'
#' @param transcripts A \code{\link[GenomicRanges]{GRanges}} or
#'   \code{\link[GenomicRanges]{GRangesList}} object containing exons of transcripts for which
#'   sequences should be extracted. All exons in a \code{\link[GenomicRanges]{GRanges}} object are
#'   assumed to belong to the same transcript. Multiple transcripts can be provided in a
#'   \code{\link[GenomicRanges]{GRangesList}} object.
#' @param  genome A \code{\link[BSgenome]{BSgenome}} or \code{\link[Biostrings]{DNAStringSet}}
#'   object containing chromosome sequences which should be used to extract transcript sequences.
#' @return A \code{\link[Biostrings]{DNAString}} or \code{\link[Biostrings]{DNAStringSet}} object
#'   containing the transcript sequence(s).
#' @examples
#' library(TAPseq)
#' library(Biostrings)
#'
#' # protein-coding exons of transcripts within chr11 region
#' data("chr11_genes")
#' target_txs <- split(chr11_genes, f = chr11_genes$transcript_id)
#'
#' # load human chr11 sequence from fasta file
#' chr11_seq_fasta <- system.file("extdata", "chr11_sequence.fasta.gz", package = "TAPseq")
#' chr11_seq <- readDNAStringSet(chr11_seq_fasta)
#'
#' # get sequences for all target transcripts in chr11 region
#' tx_seqs <- getTxsSeq(target_txs, genome = chr11_seq)
#'
#' \dontrun{
#' # a probably more practical approach is to use Bioconductors BSgenome package to get the
#' # genome sequence
#' library(BSgenome)
#'
#' # human genome (hg38) BSgenome object (needs to be istalled separately from Bioconductor)
#' hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
#'
#' # get sequences for all target transcripts on chr11
#' tx_seqs <- getTxsSeq(target_txs, genome = hg38)
#'}
#' @export
setGeneric("getTxsSeq", function(transcripts, genome) standardGeneric("getTxsSeq") )

#' @describeIn getTxsSeq Obtain transcript sequence from \code{GRangesList} input
#' @export
setMethod("getTxsSeq", "GRangesList", function(transcripts, genome) {

  # abort if genome is not a BSgenome or DNAStringSet object
  if (class(genome) != "BSgenome" & class(genome) != "DNAStringSet") {
    stop("genome must be of class BSgenome or DNAStringSet!", call. = FALSE)
  }

  # make sure that all chromosomes in transcripts are found in genome
  txs_chrs <- as.character(unique(seqnames(unlist(transcripts))))
  if (length(setdiff(txs_chrs, names(genome))) > 0) {
    stop("Not all chromosomes in transcripts found in genome object!", call. = FALSE)
  }

  # get indices of transcripts on positive and negative strand
  txs_pos <- which(all(strand(transcripts) == "+"))
  txs_neg <- which(all(strand(transcripts) == "-"))

  # get transcripts with incorrect or conflicting strand information
  bad_txs <- setdiff(1:length(transcripts), sort(c(txs_pos, txs_neg)))

  # abort if any are found
  if (length(bad_txs) > 0) {
    if (!is.null(names(transcripts))) {
      stop("Incorrect strand information in transcripts! ",
           "Strand of all exons per transcript must be '+' or '-' ('*' is not allowed).\n",
           "Check strand for transcript(s): ",
           paste(names(transcripts)[bad_txs], collapse = ", "), call. = FALSE)
    }else{
      stop("Incorrect strand information in transcripts! ",
           "Strand of all exons per transcript must be '+' or '-' ('*' is not allowed).",
           call. = FALSE)
    }
  }

  # order exons of each transcript according to order in transcript (5' -> 3')
  transcripts[txs_pos] <- sort(transcripts[txs_pos], decreasing = FALSE)
  transcripts[txs_neg] <- sort(transcripts[txs_neg], decreasing = TRUE)

  # get sequences of transcripts
  GenomicFeatures::extractTranscriptSeqs(genome, transcripts = transcripts)

}
)

#' @describeIn getTxsSeq Obtain transcript sequence from \code{GRanges} input
#' @export
setMethod("getTxsSeq", "GRanges", function(transcripts, genome) {

  # transform transcripts to GRangesList
  transcripts <- GenomicRanges::GRangesList(transcripts)

  # obtain transcript sequence
  seq <- getTxsSeq(transcripts, genome = genome)

  # return transcript sequence as DNAString
  seq[[1]]

}
)
