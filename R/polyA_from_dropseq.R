#' Infer polyA sites from droplet sequencing data
#'
#' Infer polyA sites from 10X, Drop-seq or similar 3' enriched sequencing data. Simple function
#' that looks for peaks in read coverage to estimate potential polyA sites. Default parameters are
#' chosen because they work reasonably well with the example data, but they should typically be
#' empirically selected by verifying the output.
#'
#' @param genes \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} object containing
#'   annotations of genes for which polyA sites are to be estimated.
#' @param bam Path to .bam file containing aligned reads used for polyA site estimation.
#' @param polyA_downstream (numeric) How far downstream of a peak in coverage are polyA sites
#'   expected? Somewhat depends on input DNA fragment size. (default: 100).
#' @param min_cvrg (numeric) Minimal coverage for peaks to be considered for polyA site estimation
#'   (default: 0).
#' @param wdsize (numeric) Window size to estimate sequencing coverage along transcripts (default:
#'   200).
#' @param by (numeric) Steps in basepairs in which the sliding window should be moved along
#'   transcripts to estimate smooth coverage (default: 1).
#' @param extend_downstream (numeric) To which amount should transcript annotations be extended
#'   downstream when estimating polyA sites (default: 0). A reasonable value (e.g. 100-200 bp)
#'   allows to account for polyA sites that fall a few basepairs downstream of terminal exons.
#' @param perc_threshold (numeric) Only sequencing coverage peaks within \code{perc_threshold}
#'   percentile of coverage are considered for polyA site estimation (default: 0.9). Avoids that
#'   small peaks that in coverage are considered, resulting in manby false polyA sites.
#' @param parallel (logical) Triggers parallel computing using the
#'   \code{\link[BiocParallel]{BiocParallel-package}} package. This requires that a parallel
#'   back-end was registered prior to executing the function. (default: FALSE).
#' @return A \code{\link[GenomicRanges:GRanges-class]{GRanges}} object containing coordinates of estimated polyadenylation sites.
#' @examples
#' library(GenomicRanges)
#'
#' # protein-coding exons of genes within chr11 region
#' data("chr11_genes")
#' target_genes <- split(chr11_genes, f = chr11_genes$gene_name)
#'
#' # subset of target genes for quick example
#' target_genes <- target_genes[18:27]
#'
#' # bam file containing aligned Drop-seq reads
#' dropseq_bam <- system.file("extdata", "chr11_k562_dropseq.bam", package = "TAPseq")
#'
#' # infer polyA sites for all target genes with adjusted parameters. parameter values depend on the
#' # input data and at this stage it's best to try different settings and check the results
#' polyA_sites <- inferPolyASites(target_genes, bam = dropseq_bam, polyA_downstream = 50,
#'                                wdsize = 100, min_cvrg = 1, parallel = TRUE)
#' @export
inferPolyASites <- function(genes, bam, polyA_downstream = 100, min_cvrg = 0, wdsize = 200,
                            by = 1, extend_downstream = 0, perc_threshold = 0.9, parallel = FALSE) {

  # verify that genes is a GRangesList
  if (!is(genes , "GRangesList")) {
    stop("genes needs to be a GRangesList containing annotations per gene!", call. = FALSE)
  }

  # calculate per base coverage for all genes ------------------------------------------------------

  ## TO DO: load bam only for chromosomes of target genes

  # load all reads in bam file
  reads <- GenomicAlignments::readGAlignments(bam)

  # calculate coverage per base
  read_cvrg <- GenomicAlignments::coverage(reads)

  # find putative polyA sites ----------------------------------------------------------------------

  # get polyA sites for all genes
  if (parallel == TRUE) {

    polyA_sites <- BiocParallel::bplapply(X = genes, FUN = polyA_site_gene, coverage = read_cvrg,
                          polyA_downstream = polyA_downstream, wdsize = wdsize, by = by,
                          extend_downstream = extend_downstream, perc_threshold = perc_threshold)
  } else {

    polyA_sites <- lapply(X = genes, FUN = polyA_site_gene, coverage = read_cvrg,
                          polyA_downstream = polyA_downstream, wdsize = wdsize, by = by,
                          extend_downstream = extend_downstream, perc_threshold = perc_threshold)

  }

  # transform output to GRanges object
  polyA_sites <- unlist(GRangesList(polyA_sites))

  # filter polyA sites according to min_cvrg
  polyA_sites <- polyA_sites[polyA_sites$score >= min_cvrg]

  # sort according to position
  polyA_sites <- sort(polyA_sites)

  # return output
  return(polyA_sites)

}

# helper functions =================================================================================

# infer polyA sites for one gene
polyA_site_gene <- function(gene, coverage, polyA_downstream, wdsize, by, extend_downstream,
                            perc_threshold) {

  # process annotation and coverage data -----------------------------------------------------------

  # get chromosome and strand of gene
  chr <- as.character(unique(seqnames(gene)))
  strand <- as.character(unique(strand(gene)))

  # merge exons of gene
  exons <- reduce(gene)

  # find putative polyA site if gene is on the '+' strand
  if (strand == "+") {

    # extend terminal exon downstream
    end(exons[length(exons)]) <- end(exons[length(exons)]) + extend_downstream

    # exon ends and starts in spliced (concatenated) transcript
    exon_ends <- cumsum(width(exons))
    exon_starts <- c(1, exon_ends[-length(exon_ends)] + 1)
    tx_exons <- IRanges(exon_starts, exon_ends)

    # calculate coverage in sliding windows along gene ---------------------------------------------

    # coverage of gene
    gene_cvrg <- Views(coverage[[chr]], ranges(exons))

    # coverage of spliced transcript (concatenated)
    tx_cvrg <- unlist(gene_cvrg)

    # define windows based on strand of gene, so that the 3' end is for sure
    # covered by a window
    wd_ends <- seq(from = sum(width(exons)), to = wdsize, by = -by)
    wd_starts <- wd_ends - wdsize + 1
    wds <- IRanges(wd_starts, wd_ends)

    # calculate coverage in each window
    wd_cvrg <- sum(Views(tx_cvrg, wds))

    # define putative polyA site based on local maxima (peaks) -------------------------------------

    # calculate local maxima in smoothed coverage
    peaks <- local_maxima(wd_cvrg)
    peaks_cvrg <- wd_cvrg[peaks]

    # retain peaks which coverage is >= perc_threshold
    min_cvrg <- as.numeric(stats::quantile(wd_cvrg, probs = perc_threshold))
    peaks <- peaks[peaks_cvrg >= min_cvrg]

    # get windows for these peaks
    peak_wds <- wds[peaks]

    # get center of peaks
    peak_centers <- round(start(peak_wds) + wdsize / 2 - 1)

    # define putative polyA site
    polyA_sites <- peak_centers + polyA_downstream

    # if gene is on '-' strand
  } else if (strand == "-") {

    start(exons[1]) <- start(exons[1]) - extend_downstream

    exon_ends <- cumsum(width(exons))
    exon_starts <- c(1, exon_ends[-length(exon_ends)] + 1)
    tx_exons <- IRanges(exon_starts, exon_ends)

    # calculate coverage in sliding windows along gene ---------------------------------------------

    gene_cvrg <- Views(coverage[[chr]], ranges(exons))

    tx_cvrg <- unlist(gene_cvrg)

    wd_starts <- seq(from = 1, to = (sum(width(exons)) - wdsize + 1), by = by)
    wd_ends <- wd_starts + wdsize - 1
    wds <-IRanges(wd_starts, wd_ends)

    wd_cvrg <- sum(Views(tx_cvrg, wds))

    # define putative polyA site based on local maxima (peaks) -------------------------------------

    peaks <- local_maxima(wd_cvrg)
    peaks_cvrg <- wd_cvrg[peaks]

    min_cvrg <- as.numeric(stats::quantile(wd_cvrg, probs = perc_threshold))
    peaks <- peaks[peaks_cvrg >= min_cvrg]

    peak_wds <- wds[peaks]

    peak_centers <- round(start(peak_wds) + wdsize / 2 - 1)

    polyA_sites <- peak_centers - polyA_downstream + 1

  } else {

    stop("Strand not '+' or '-' for at least 1 gene!", call. = FALSE)

  }

  # translate polyA sites to genomic coordinates
  if (length(polyA_sites) > 0) {

    # transform to IRanges object
    polyA_sites <- IRanges(start = polyA_sites, end = polyA_sites)

    # add metadata
    cov_bp <- wd_cvrg[peaks] / wdsize # coverage per bp of peaks
    mcols(polyA_sites) <- data.frame("name" = "polyA_site",
                                     "score" = cov_bp)

    # sort and remove possible duplicate polyA sites due to rounding of coordinates
    polyA_sites <- unique(sort(polyA_sites))

    # get exon overlapping with polyA sites
    polyA_overlaps <- findOverlaps(query = tx_exons, subject = polyA_sites)
    polyA_exons <- S4Vectors::queryHits(polyA_overlaps)

    # only retain polyA sites that overlap any exons
    polyA_sites <- polyA_sites[S4Vectors::subjectHits(polyA_overlaps)]

    # calculate distance of polyA site to the start of the overlapping exon
    pos_diff <- start(polyA_sites) - start(tx_exons[polyA_exons])

    # translate polyA site coordinates to genomic coordinates
    polyA_start <- start(exons[polyA_exons]) + pos_diff

    # create GRanges object with genomic coordinates of polyA site
    polyA_coords <- GRanges(seqnames = chr, strand = strand,
                            ranges = IRanges(polyA_start, polyA_start),
                            mcols(polyA_sites))

  } else {

    # if no polyA sites are found, return empty GRanges object
    polyA_coords <- GRanges()

  }

  # return output
  return(polyA_coords)

}

## FUN: find local maxima in an ordered data series (after post in:
#  https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima)
local_maxima <- function(x) {

  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-1, x)) > 0L
  # y <- diff(c(-.Machine$integer.max, x)) > 0L

  rle(y)$lengths

  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]

  if (x[[1]] == x[[2]]) {

    y <- y[-1]

  }

  # return vector with indices of local maxima
  return(y)

}
