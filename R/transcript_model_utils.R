# --------------------------------------------------------------------------------------------------
# This file contains various utility functions for truncating transcript models by
# truncate_txs_polyA(). These functions are not exported to the NAMESPACE and are intended to be
# used as internal functions. Serious errors can occur when these functions are used out of their
# context. YOU HAVE BEEN WARNED!
# --------------------------------------------------------------------------------------------------

# truncate transcripts of a gene at an overlapping polyA site (if any are found). polyA_select is
# used to control how the polyA site for truncation is selected if multiple polyA sites overlap with
# the gene. if the selected polyA site overlaps with multiple annotated transcripts, truncation is
# performed on a meta-transcript created by merging all transcripts overlapping with the site. if no
# polyA site overlaps with the gene, a meta-transcript containing all merged transcripts is returned
truncate_tx_polyA <- function(transcripts, polyA_sites, extend_3prime_end = 0,
                              polyA_select = c("downstream", "upstream", "score"),
                              transcript_id = "transcript_id", gene_id = "gene_id",
                              exon_number = "exon_number", ignore_strand = FALSE) {

  # get polyA_select argument
  polyA_select <- match.arg(polyA_select)

  # split transcripts by transcript id
  txs <- split(transcripts, f = mcols(transcripts)[[transcript_id]])

  # extend 3' ends of the terminal exons
  txs_ext <- endoapply(txs, FUN = extend_tx, downstream = extend_3prime_end)

  # find polyA sites overlapping transcripts
  overlaps_sites <- findOverlaps(txs_ext, polyA_sites, ignore.strand = ignore_strand)

  # extract overlaping sites
  site_ids <- unique(subjectHits(overlaps_sites))
  sites <- polyA_sites[site_ids]
  sites$id <- site_ids

  # truncate transcripts if at least one overlapping site is found
  if (length(sites) > 0) {

    # order sites according to strand of transcripts
    if (all(strand(transcripts) == "-")) {
      sites <- sort(sites, decreasing = TRUE)
    }else{
      sites <- sort(sites)
    }

    # select polyA site for truncation based on polyA_select
    if (polyA_select == "downstream") {
      site <- sites[length(sites)]
    }else if (polyA_select == "upstream") {
      site <- sites[1]
    }else{
      site <- sites[which.max(sites$score)]
    }

    # select transcripts overlapping with the selected site
    overlapping_txs <- txs_ext[queryHits(overlaps_sites)[subjectHits(overlaps_sites) == site$id]]

    # if more than 1 transcript overlap the site it can't be specified which transcript overlaps
    # the site. therefore overlapping exons of all transcripts overlapping the site are merged to
    # derive a consensus transcript model.
    if (length(overlapping_txs) > 1) {
      tx <- reduce_gene(unlist(overlapping_txs), transcript_id = transcript_id, gene_id = gene_id,
                        exon_number = exon_number)
    }else{
      tx <- sort(overlapping_txs[[1]])
    }

    # truncate overlapping transcript at selected polyA site
    if (all(strand(tx) == "-")) {
      restrict(tx, start = end(site)[length(site)] + 1, end = end(tx)[length(tx)])
    }else{
      restrict(tx, start = start(tx)[1], end = start(site)[1] - 1)
    }

  # if no overlapping polyA sites are found, no statement about the expressed transcript can be made
  # and therefore the consensus (merge) of all transcripts is returned
  }else if (length(unique(mcols(transcripts)[[transcript_id]])) == 1) {
     return(transcripts)
  }else{
    reduce_gene(transcripts, transcript_id = transcript_id, gene_id = gene_id,
                exon_number = exon_number)
  }
}

# extend a provided transcript model at the 3' (downstream) or 5' (upstream) by the specified number
# of base pairs. Only extends the terminal exons of the transcript
extend_tx <- function(tx, upstream = 0, downstream = 0) {

  # make sure that tx is sorted correctly
  tx <- sort(tx)

  # extend terminal exons depending on the strand
  if (all(strand(tx) == "-")) {
    start(tx[1]) <- start(tx[1]) - downstream
    end(tx[length(tx)]) <- end(tx[length(tx)]) + upstream
  }else{
    end(tx[length(tx)]) <- end(tx[length(tx)]) + downstream
    start(tx[1]) <- start(tx[1]) - upstream
  }

  # return modifed tx
  return(tx)

}

# Merge overlapping exons of a gene in GRanges format. Basically calls GenomicRanges::reduce, but
# keeps metadata (id unique for all exons). gene_id, transcript_id and exon_number are names of
# expected columns in metadata for these features. If found in input gene's metadata, they will be
# used to create a new transcript id and exon numbers. If not found in input gene's metadata, these
# columns simply won't appear in the output metadata.
reduce_gene <- function(gene, gene_id = "gene_id", transcript_id = "transcript_id",
                        exon_number = "exon_number") {

  # reduce ranges of gene to create merged annotations
  gene_red <- reduce(gene)

  # get DataFrame containing metadata columns for input gene annotations
  metadat_df <- mcols(gene)

  # get columns with unique values and transform to list
  metadat <- lapply(X = metadat_df, FUN = unique)

  # set columns with non-unique values to NA
  metadat[vapply(metadat, FUN = length, FUN.VALUE = integer(1)) > 1] <- as.integer(NA)

  # add new transcript_id (gene_id + "-MERGED") to identify merged transcripts
  gid <- as.character(metadat[[gene_id]])
  metadat[[transcript_id]] <- paste(c(gid, "MERGED"), collapse = "-")

  # create exon number and add to metadata
  exon_num <- seq_len(length(gene_red))
  if (all(strand(gene_red) == "-")) {
    exon_num <- sort(exon_num, decreasing = TRUE)
  }
  metadat[[exon_number]] <- exon_num

  # make sure that metadat only contains data found in input and as same classes (except for tx id)
  metadat_df_classes <- vapply(metadat_df, FUN = class, FUN.VALUE = character(1))
  metadat_df_classes[[transcript_id]] <- "character"
  metadat_df_classes <- metadat_df_classes[names(metadat_df)]
  metadat <- metadat[names(metadat_df)]
  metadat <- mapply(FUN = as, metadat, metadat_df_classes, SIMPLIFY = FALSE)

  # add metadat to gene_red and return output
  mcols(gene_red) <- DataFrame(metadat)
  return(gene_red)

}
