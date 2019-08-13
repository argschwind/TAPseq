#' Truncate transcripts at polyA sites
#'
#' Truncate transcripts at overlapping polyadenylation (polyA) sites to infer  likely 3' ends of
#' transcripts. This is crucial to correctly design TASC-seq primers that amplify fragments of
#' specific lengths. Typically the exons of all annotated transcripts per gene are provided as
#' input. If a polyA site overlaps a single transcript of a given gene, this transcript is truncated
#' and returned. In case a polyA site overlaps multiple transcripts of the same gene, a
#' "metatranscript" consisting of all annotated exons of the overlapping transcripts is generated
#' and truncated. No statements about expressed transcripts can be made if no overlapping polyA
#' sites are found for any transcripts of a gene. In that case a "metatranscript" consisting of
#' the merged exons of that gene is generated and returned.
#'
#' @param transcripts A \code{GRanges} or \code{GRangesList} object containing  exons of the
#'   transcripts to be truncated.
#' @param polyA_sites A \code{GRanges} object containing the polyA sites. This needs to contain a
#'   metadata entry names "score" if the option \code{polyA_select = "score"} is used. PolyA sites
#'   can be either obtained via running \code{\link[TASCseq]{inferPolyASites}} or imported from
#'   an existing bed file via \code{\link[rtracklayer]{import}}.
#' @param extend_3prime_end Specifies how far (bp) 3' ends of transcripts should be extended when
#'   looking for overlapping polyA sites (default = 0). This enables capturing of polyA sites that
#'   occur downstream of annotated 3' ends.
#' @param polyA_select Specifies which eurisic should be used to select the polyA site used to
#'   truncate the transcripts if multiple overlapping polyA sites are found. By default
#'   \code{"downstream"} is used which choses the most downstream polyA site. \code{"score"} selects
#'   the polyA site with the highest score, which correspons to the read coverage when using
#'   \code{\link[TASCseq]{inferPolyASites}} to estimate polyA sites.
#' @param ignore_strand (logical) Specifies whether the strand of polyA sites should be ignored when
#'   looking for overlapping polyA sites. Default is \code{FALSE} and therefore only polyA sites on
#'   the same strand as the transcripts are considered. PolyA sites with strand equals to \code{*}
#'   has the same effect as \code{ignore_strand = TRUE}.
#' @param parallel (logical) Triggers parallel computing using the \code{BiocParallel} package.
#'   This requires that a parallel back-end was registered prior to executing the function.
#'   (default: FALSE).
#'
#' @return Either a \code{GRanges} or \code{GRangesList} object containing the
#'   truncated transcripts.
#'
#' @examples
#' library(TASCseq)
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
                    ignore_strand = FALSE, parallel = FALSE)
             standardGeneric("truncateTxsPolyA")
)

#' @describeIn truncateTxsPolyA Truncate transcripts of one gene provided as \code{GRanges} object
#' @export
setMethod("truncateTxsPolyA", "GRanges",
          function(transcripts, polyA_sites, extend_3prime_end = 0,
                   polyA_select = c("downstream", "upstream", "score"),
                   ignore_strand = FALSE) {

    # abort if input data has wrong format
    if (is.null(transcripts$transcript_id)) {
      stop("transcript_id metadata field required!", call. = FALSE)
    }

    if (class(polyA_sites) != "GRanges"){
      stop("polyA_sites needs to be of class GRanges!", call. = FALSE)
    }

    # get polyA_select argument
    polyA_select <- match.arg(polyA_select)

    # make sure that transcripts are sorted correctly
    transcripts <- sort(transcripts)

    # split transcripts by transcript id
    txs <- split(transcripts, f = transcripts$transcript_id)

    # extend 3' ends of the terminal exons
    txs_ext <- GenomicRanges::GRangesList(
      lapply(X = txs, FUN = extend_tx, downstream = extend_3prime_end)
    )

    # find transcripts and polyA sites that overlap
    overlaps <- as.matrix(
      GenomicRanges::findOverlaps(txs_ext, polyA_sites, ignore.strand = ignore_strand)
    )

    # truncate transcripts based on overlapping polyA sites
    if (nrow(overlaps) > 0) {

      # get overlapping polyA sites
      overlapping_sites <- sort(polyA_sites[unique(overlaps[,2])])

      # select polyA site based on polyA_select
      if (polyA_select == "score") {

        # select polyA site with highest score
        max_score <- max(overlapping_sites$score)
        site <- overlapping_sites[overlapping_sites$score == max_score]

        # if there are multiple sites, only retain the most upstream site
        if (all(BiocGenerics::strand(transcripts) == "+")) {
          site <- site[1]
        }else{
          site <- site[length(site)]
        }

      # select most downstream polyA site
      }else if (polyA_select == "downstream") {
        if (all(BiocGenerics::strand(transcripts) == "+")) {
          site <- overlapping_sites[length(overlapping_sites)]
        }else{
          site <- overlapping_sites[1]
        }

      # select most upstream polyA site
      }else{
        if (all(BiocGenerics::strand(transcripts) == "+")){
          site <- overlapping_sites[1]
        }else{
          site <- overlapping_sites[length(overlapping_sites)]
        }
      }

      # select transcripts that overlap this site
      txs_overlapping <- IRanges::subsetByOverlaps(txs_ext, site,
                                                   ignore.strand = ignore_strand
      )

      # if more than 1 transcript overlap the site it can't be specified which transcript overlaps
      # the site. therefore overlapping exons of all transcripts overlapping the site are merged to
      # derive a consensus transcript model.
      if (length(txs_overlapping) > 1) {

        # merge overlapping exons
        tx_overlapping <- reduce_gene(unlist(txs_overlapping))

      }else{

        # else pick the single transcript overlapping the site
        tx_overlapping <- txs_overlapping[[1]]

      }

      # truncate overlapping transcript
      output <- truncate_tx(tx = tx_overlapping, site = site, ignore_strand = ignore_strand)

      # if no overlaps are found, no statement about the expressed transcript can
      # be made and therefore the consensus (merge) of all transcripts is
      # returned
      }else{

      # merge overlapping exons
      output <- reduce_gene(transcripts)

    }

    # return output
    return(output)
  }
)

#' @describeIn truncateTxsPolyA Truncate transcripts of multiple genes provided as
#'   \code{GRangesList}
#' @export
setMethod("truncateTxsPolyA", "GRangesList",
          function(transcripts, polyA_sites, extend_3prime_end = 0,
                   polyA_select = c("downstream", "upstream", "score"), ignore_strand = FALSE,
                   parallel = FALSE) {

     # apply truncateTxsPolyA.GRanges to each GRanges object in GRangesList
     if (parallel == TRUE){

       output <- BiocParallel::bplapply(X = transcripts, FUN = truncateTxsPolyA,
                                        polyA_sites = polyA_sites,
                                        extend_3prime_end = extend_3prime_end,
                                        polyA_select = polyA_select, ignore_strand = ignore_strand)

     }else{

       # truncate transcripts
       output <- lapply(X = transcripts, FUN = truncateTxsPolyA, polyA_sites = polyA_sites,
                        extend_3prime_end = extend_3prime_end, polyA_select = polyA_select,
                        ignore_strand = ignore_strand)

     }

     # transform output to GRangesList
     GenomicRanges::GRangesList(output)

   }
)
