#' Export TAP-seq primers
#'
#' A set of functions for TAP-seq primer export. Convert primers stored in
#' \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} objects to a simple
#' \code{data.frame} for easier export. Or create BED format tracks for primers and write
#' them to files for viewing in a genome browser (e.g. IGV).
#'
#' @param object A \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object containing
#'   designed primers.
#' @param color Color used for the track (Default: black). Can be any of the three kinds of R color
#'   specifications.
#' @return For \code{createPrimerTrack} a \code{data.frame} with the primer track in BED format
#' @examples
#' library(TAPseq)
#'
#' # chr11 primers example data
#' data("chr11_primers")
#'
#' # pick best primers based on predicted off-targets
#' best_primers <- pickPrimers(chr11_primers, n = 1, by = "off_targets")
#'
#' # primers data can be exported to a simple data.frame to e.g. write them to a .csv file
#' primers_df <- primerDataFrame(best_primers)
#' head(primers_df)
#'
#'
#' # primer binding sites in transcript sequences can be converted to genomic coordinates to create
#' # a BED track to visualize primers in a genome browser (e.g. IGV)
#'
#' # create primer BED track with a fancy color
#' track <- createPrimerTrack(best_primers[1:5], color = "steelblue3")
#'
#' # tracks can be written to .bed files using a little helper function (replace con = "" by a file)
#' exportPrimerTrack(track, con = "")
#'
#' \dontrun{
#' # one can easily export primer tracks for multiple TsIO or TsIOList objects (e.g. inner and
#' # outer nested primers) to one .bed file using different colors for each object. see vignette for
#' # a practical example:
#' vignette("tapseq_primer_design", package = "TAPseq")
#'
#' obj1 <- best_primers[1:5]
#' obj2 <- best_primers[6:10]
#' exportPrimerTrack(createPrimerTrack(obj1, color = "steelblue3"),
#'                   createPrimerTrack(obj2, color = "goldenrod1"),
#'                   con = "path/to/file.bed")
#'
#' }
#' @name exportPrimers
NULL

#' @rdname exportPrimers
#' @export
setGeneric("createPrimerTrack", function(object, color = 1) standardGeneric("createPrimerTrack") )

#' @describeIn exportPrimers Create primer BED track from \code{TsIO} objects
#' @export
setMethod("createPrimerTrack", "TsIO", function(object, color) {

  # check that object has valid target annotations
  annot <- target_annot(object)
  if (length(annot) == 0) {
    stop("Input object does not contain valid target annotations", call. = FALSE)
  }

  # parse color argument
  color <- paste0(grDevices::col2rgb(color)[, 1], collapse = ",")

  # get genomic coordinates of all primers
  primer_coords <- calc_coords_primers(object)

  # create bed tracks
  bed_tracks <- lapply(primer_coords, FUN = create_bed_track, rgb_col = color)
  bed_tracks <- bind_rows(bed_tracks)
  as.data.frame(bed_tracks)

})

#' @describeIn exportPrimers Create primer BED track from \code{TsIOList} objects
#' @export
setMethod("createPrimerTrack", "TsIOList", function(object, color) {

  # check that all objects have valid target annotations
  annot <- target_annot(object)
  n_exons <- vapply(annot, FUN = length, FUN.VALUE = integer(1))
  if (any(n_exons == 0)) {
    stop("Not all TsIO objects in input valid contain target annotations", call. = FALSE)
  }

  # parse color argument
  color <- paste0(grDevices::col2rgb(color)[, 1], collapse = ",")

  # get genomic coordinates of all primers
  primer_coords <- lapply(object, FUN = calc_coords_primers)
  primer_coords <- unlist(primer_coords)

  # create bed tracks
  bed_tracks <- lapply(primer_coords, FUN = create_bed_track, rgb_col = color)
  bed_tracks <- bind_rows(bed_tracks)
  as.data.frame(bed_tracks)

})

#' @describeIn exportPrimers Export primer BED tracks files
#'
#' @param con Connection to which tracks are written. Typically a .bed file.
#' @param ... One or more primer BED tracks created by \code{\link[TAPseq]{createPrimerTrack}}.
#' @export
exportPrimerTrack <- function(..., con) {

  # process arguments
  if (missing(con)) {
    args <- list(...)
    if (length(args) == 0) {
      stop("Nothing to export...", call. = FALSE)
    }else if (length(args) == 1) {
      stop('Argument "con" is missing, with no default', call. = FALSE)
    }else{
      con <- args[[length(args)]]
      tracks <- args[-length(args)]
    }
  }else{
    tracks <- list(...)
    if (length(tracks) == 0) {
      stop("Nothing to export...", call. = FALSE)
    }
  }

  # write header line to bed file
  header <- paste('track name="Tardrop primers" description="TAP-seq primers for targeted',
                  'single-cell transcriptomics" itemRgb="On"')
  write(header, file = con)

  # write tracks to file
  tracks <- bind_rows(tracks)
  utils::write.table(tracks, file = con, append = TRUE, quote = FALSE, row.names = FALSE,
                     col.names = FALSE, sep = "\t")

}

#' @rdname exportPrimers
#'
#' @param object A \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object containing
#'   designed primers.
#' @export
setGeneric("primerDataFrame", function(object) standardGeneric("primerDataFrame") )

#' @describeIn exportPrimers Create a \code{data.frame} with primer data from \code{TsIO}
#' @export
setMethod("primerDataFrame", "TsIO", function(object) {
  create_primer_df(object)
})

#' @describeIn exportPrimers Create a \code{data.frame} with primer data from \code{TsIOList}
#' @export
setMethod("primerDataFrame", "TsIOList", function(object) {
  output <- lapply(object, FUN = create_primer_df)
  bind_rows(output)
})

# HELPER FUNCTIONS =================================================================================

# get primer coordinates for primers in one TsIO object
calc_coords_primers <- function(object) {

  # get designed primers and target sequence and annotations
  primers <- tapseq_primers(object)
  annot <- target_annot(object)
  seq <- target_sequence(object)

  # the start and end positions of the genes exons within the transcript are inferred based on the
  # provided gene annotation. the cumulative sum of exon lengths is used to get end position of
  # exon junctions within the transcript
  exon_ends <- cumsum(width(annot))
  exon_starts <- c(1, (exon_ends[-length(exon_ends)] + 1))
  exons <- IRanges(start = exon_starts, end = exon_ends)

  # infer genomic coordinates of primer binding sites
  primers_list <- split(primers, f = names(primers))
  lapply(primers_list, FUN = calc_coords_one_primer, exons = exons, annot = annot)

}

# calculate genomic coordinates of primer binding sites based on provided transcript annotations
calc_coords_one_primer <- function(primer, exons, annot) {

  # find exons that overlap with the primers
  overlaps <- IRanges::findOverlaps(query = primer, subject = exons)
  primer_exons <- subjectHits(overlaps)

  # split primer site by overlapping exons in case the primer overlaps an exon-exon junction
  primer_site <- IRanges::pintersect(IRanges::findOverlapPairs(primer, exons))

  # the primer binding site position(s) have to be transformed to genomic coordinates. to do this,
  # the positions of the primer site(s) relative to the start of the exon where the primer binds
  # are added to the genomic coordinates of the exon start
  if (all(strand(annot) == "-")) {
    primer_start <- start(annot[primer_exons]) + end(exons[primer_exons]) - end(primer_site)
    primer_end <- start(annot[primer_exons]) + end(exons[primer_exons]) - start(primer_site)
  }else{
    primer_start <-  start(annot[primer_exons]) + start(primer_site) - start(exons[primer_exons])
    primer_end <- start(annot[primer_exons]) + end(primer_site) - start(exons[primer_exons])
  }

  # create GRanges object for genomic primer binding site
  primer_coords <- GRanges(seqnames = as.character(unique(seqnames(annot))),
                           ranges = IRanges(primer_start, primer_end),
                           strand = unique(strand(annot)))

  # add metadata and return output
  mcols(primer_coords) <- cbind(primer_id = names(primer), mcols(primer))
  return(primer_coords)

}

# create bed track for one primer binding site
create_bed_track <- function(primer_coords, rgb_col) {

  # sort coordinates in case there are multiple sites (should already be sorted)
  primer_coords <- sort(primer_coords)

  # extract metadata
  metadat <- unique(mcols(primer_coords))

  # create .bed fields as a data.frame with one row
  data.frame(stringsAsFactors = FALSE,

    # genomic coordinates of the primer
    chrom = as.character(unique(seqnames(primer_coords))),
    start = start(primer_coords)[1] - 1,
    end = end(primer_coords)[length(primer_coords)],

    # attributes of the primer
    name =  metadat$primer_id,
    score = metadat$penalty,
    strand = as.character(unique(strand(primer_coords))),

    # display the whole primer as thick line (shows '>>' for strand)
    thickStart = start(primer_coords)[1] - 1,
    thickEnd = end(primer_coords)[length(primer_coords)],

    # color of the track
    itemRgb = rgb_col,

    # the regions of the track where the primer binds are displayed as 'blocks'. some primers span
    # exon junctions (i.e. introns in the genomic annotation), and this allows these junctions to be
    # displayed as thin lines similar to introns of gene annotations.
    blockCount = length(primer_coords),
    blockSizes = paste(width(primer_coords), collapse = ","),
    blockStarts = paste(start(primer_coords) - start(primer_coords)[1], collapse = ",")

  )

}

# create a data.frame containing primers for export as .csv file from a TsIO object
create_primer_df <- function(object) {

  # get data stored in object
  seq_id <- sequence_id(object)
  target_seq <- target_sequence(object)
  primers <- tapseq_primers(object)
  pcr_prods <- pcr_products(object)

  # return empty data.frame if no primers are found
  if (length(primers) == 0) return(data.frame())

  # convert primers to data.frame
  primers_df <- as.data.frame(primers)

  # get target sequence length and expected pcr product sizes
  seq_len <-length(target_seq)
  pcr_len <- data.frame(primer_id = names(pcr_prods), pcr_product_size = width(pcr_prods),
                        stringsAsFactors = FALSE)

  # create output data.frame
  output <- data.frame(seq_id = seq_id, seq_len = seq_len, primers_df[, 1:12],
                       stringsAsFactors = FALSE)
  output <- dplyr::rename(output, primer_id = names, primer_len = width)
  output <- left_join(output, pcr_len, by = "primer_id")

  # add off-targets to output if provided
  off_target_cols <- c("intergenic_off_targets", "intronic_off_targets", "exonic_off_targets")
  off_targets <- colnames(primers_df) %in% off_target_cols
  cbind.data.frame(output, primers_df[, off_targets], stringsAsFactors = FALSE)

}
