# --------------------------------------------------------------------------------------------------
# This file contains various utility functions for truncating transcript models by
# truncate_txs_polyA(). These functions are not exported to the NAMESPACE and are intended to be
# used as internal functions.
# --------------------------------------------------------------------------------------------------

# Extend a provided transcript model at the 3' (downstream) or 5' (upstream)
# by the specified number of base pairs. Only operates on the terminal exons of
# the transcript.
extend_tx <- function(tx, upstream = 0, downstream = 0){

  # make sure that all exons are on the same chromosome
  if (length(unique(GenomeInfoDb::seqnames(tx))) > 1){

    stop("Transcript maps to multiple sequences, unable to extend correctly!")

  }

  # make sure that exons do not overlap
  if (length(GenomicRanges::reduce(tx)) < length(tx)){

    stop("Overlaps between exons found! Verify transcript model.")

  }

  # make sure that tx is sorted correctly
  tx <- sort(tx)

  # extend terminal exons depending on the strand
  if (all(BiocGenerics::strand(tx) == "+")){

    BiocGenerics::end(tx[length(tx)]) <- BiocGenerics::end(tx[length(tx)]) + downstream
    BiocGenerics::start(tx[1]) <- BiocGenerics::start(tx[1]) - upstream

  }else if (all(BiocGenerics::strand(tx) == "-")){

    BiocGenerics::start(tx[1]) <- BiocGenerics::start(tx[1]) - downstream
    BiocGenerics::end(tx[length(tx)]) <- BiocGenerics::end(tx[length(tx)]) + upstream

  }else{

    stop("Incorrect strand information, cannot extend correctly!")

  }

  # return modifed tx
  return(tx)

}

# Overlap a transcript model with a genomic site (GRanges) and truncate it at
# the site if it overlaps any of the exons.
truncate_tx <- function(tx, site, ignore_strand = FALSE){

  # make sure that all exons are on the same chromosome
  if (length(unique(GenomeInfoDb::seqnames(tx))) > 1){

    stop("Transcript maps to multiple sequences, unable to extend correctly!")

  }

  # make sure that exons do not overlap
  if (length(GenomicRanges::reduce(tx)) < length(tx)){

    stop("Overlaps between exons found! Verify transcript model.")

  }

  # make sure that tx is sorted correctly
  tx <- sort(tx)

  # get overlap site with tx
  overlaps <- as.matrix(
    GenomicRanges::findOverlaps(tx, site, ignore.strand = ignore_strand)
    )

  # truncate transcript if an overlap is found
  if (length(overlaps)){

    # extract overlapping exon
    exon <- tx[overlaps[1, 1]]

    # calculate truncated part of cut_exon and add truncated part plus all
    # upstream exons to output
    if (all(BiocGenerics::strand(tx) == "+")){

      # get all exons upstream of truncated exon (if any)
      if (overlaps[1, 1] > 1){

        # upstream exons
        output <- tx[1:(overlaps[1, 1] - 1)]

      }else{

        # transcript truncated in first exon
        output <- GenomicRanges::GRanges()

      }

      # get the fraction of the overlapping exon upstream of the site (if there
      # is any upstream part), and add to output
      if (BiocGenerics::start(site) > BiocGenerics::start(exon)){

        # split exon at truncating site and get upstream fraction
        split_exon <- GenomicRanges::setdiff(exon, site, ignore.strand = TRUE)
        trunc_exon <- split_exon[1]

        # add meta data and strand
        S4Vectors::mcols(trunc_exon) <- S4Vectors::mcols(exon)
        BiocGenerics::strand(trunc_exon) <- BiocGenerics::strand(exon)

        # add truncated exon to output
        output <- c(output, trunc_exon)

      }

    }else if (all(BiocGenerics::strand(tx) == "-")){

      # get all exons upstream of truncated exon (if any)
      if (length(tx) > overlaps[1, 1]){

        # upstream exons
        output <- tx[(overlaps[1, 1] + 1):length(tx)]

      }else{

        # transcript truncated in first exons
        output <- GenomicRanges::GRanges()

      }

      # get the fraction of the overlapping exon upstream of the site (if there
      # is any upstream part), and add to output
      if (BiocGenerics::end(site) < BiocGenerics::end(exon)){

        # split exon at truncating site and get upstream fraction
        split_exon <- GenomicRanges::setdiff(exon, site, ignore.strand = TRUE)
        trunc_exon <- split_exon[length(split_exon)]

        # add meta data and strand
        S4Vectors::mcols(trunc_exon) <- S4Vectors::mcols(exon)
        BiocGenerics::strand(trunc_exon) <- BiocGenerics::strand(exon)

        # add truncated exon to output
        output <- c(trunc_exon, output)

      }

    }else{

      stop("Incorrect strand information, cannot extend correctly!")

    }

  }else{

    # don't alter transcript in case no overlaps were found
    message("No overlaps found, returning non-truncated tx...")
    output <- tx

    }

  # return truncated transcript
  return(output)

}

# Merge overlapping exons of a gene in GRanges format. Currently only works
# with genes formatted as loaded from GTF files. keep_source specifies whether
# the GTF source field should be kept or replaced by "merged_gene" (default).
reduce_gene <- function(gene, keep_source = FALSE){

  # get gene metadata
  metadat <- S4Vectors::mcols(gene)

  # check for mandatory GTF fields in metadata
  mandat <- c("source", "type", "score", "phase", "gene_id", "transcript_id")

  if (!all(mandat %in% colnames(metadat))){

    stop("Not all mandatory GTF fields found! Check gene input.")

  }

  # stop if gene contains multiple types (e.g. exons and cds)
  if (length(unique(gene$type)) > 1){

    stop("Gene contains multiple types! ",
         "Cannot merge features of different types.")

  }

  # reduce ranges of gene
  gene_red <- GenomicRanges::reduce(gene)

  # if no exons can be merged, return the original gene, else merge exons
  if (length(gene_red) == length(gene)){

    return(gene)

  }else{

    # get columns with unique values and transform to list
    metadat <- lapply(X = metadat, FUN = unique)

    # set columns with non-unique values to NA
    metadat[sapply(X = metadat, FUN = length) > 1] <- as.character(NA)

    ##### mandatory fields in GTF files #####

    # keep original source or replace by "merged_gene"
    if (keep_source == TRUE){

      metadat$source <- as.factor(metadat$source)

    }else{

      metadat$source <- as.factor("merged_gene")

    }

    # make sure that score is NA, because it cannot be merged
    metadat$score <- as.numeric(NA)

    # make sure that phase is integer
    metadat$phase <- as.integer(metadat$phase)

    # create new transcript id based on gene id
    metadat$transcript_id <- paste0(metadat$gene_id, "-MERGED")

    ##### common fields in GTF files #####

    # transcript name
    gene_name <- metadat$gene_name

    if (!is.null(metadat$transcript_name) & !is.null(gene_name)){

      # set new transcript name based on gene name
      metadat$transcript_name <- paste0(gene_name, "-MERGED")

    }

    # exon number
    if (!is.null(metadat$exon_number)){

      # calculate new exon number based on strand
      if (all(BiocGenerics::strand(gene_red) == "+")){

        metadat$exon_number <- as.character(seq(from = 1,
                                                to = length(gene_red), by = 1))

      }else{

        metadat$exon_number <- as.character(seq(from = length(gene_red),
                                                to = 1,by = -1))

      }

    }

    ##### finalize output #####

    # create and add new metadata
    S4Vectors::mcols(gene_red) <- metadat

    # return output
    return(gene_red)

  }

}
