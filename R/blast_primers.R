## Functions and methods to blast primer sequences to identify potential off-target priming

#' Create genome and transcriptome fasta file
#'
#' Create a fasta file containing the DNA sequences of a genome and its transcripts. This function
#' is used to create the input for \code{\link[TAPseq]{createBLASTDb}}. Currently only tested with
#' human and mouse genomes.
#'
#' @param genome The \code{\link[BSgenome]{BSgenome}} object of the genome to be used to obtain the
#'   genome and (optional) transcript sequences.
#' @param output_fasta Path to the output fasta file.
#' @param annot A \code{\link[GenomicRanges]{GRanges}} object containing all exons of transcripts to
#'   be considered. Needs to contain \code{transcript_id} and \code{gene_name} metadata columns. To
#'   correctly identify off-target primer binding events, \code{gene_name} should be the same as the
#'   names of the sequence templates used for primer design! If not specified, no transcript
#'   sequences will be included in the output fasta file.
#' @param include_genome (logical) Specifies whether the genome sequence should be included in the
#'   output fasta file.
#' @param exclude_chr_pattern Pattern used to exclude any unwanted chromosomes from genome, such as
#'   scaffolds with "_" by default. Set to NULL to include all chromosomes in genome.
#' @param compress (logical) Create a gzipped output fasta file.
#' @export
createFasta <- function(genome, output_fasta, annot = NULL, include_genome = TRUE,
                        exclude_chr_pattern = "_", compress = FALSE){

  # transcript sequences ---------------------------------------------------------------------------

  if (!is.null(annot)) {
    message("Obtaining transcript sequences...")

    # make sure that annot is a GRanges object
    if (!is(annot, "GRanges")) stop("annot needs to be a GRanges object!", call. = FALSE)

    # create a new identifier for each exon consisting of the transcript id, the type and the gene
    # name. This id will later be used as the title of each entry in the fasta file.
    tx_ids <- paste(annot$transcript_id, "transcript", annot$gene_name)

    # split exons by this new tx_ids into GRangesList with each element
    # containing all exons of a given transcript
    txs <- split(annot, f = tx_ids)

    # split transcripts according to strand
    txs_pos <- txs[all(strand(txs) == "+")]
    txs_neg <- txs[all(strand(txs) == "-")]

    # get sequences of all transcripts
    txs_seqs_pos <- getSeq(x = genome, names = txs_pos)
    txs_seqs_neg <- getSeq(x = genome, names = txs_neg)

    # concatenate (unlist) sequences of all exons per transcript
    txs_seqs_pos_cat <- DNAStringSet(lapply(X = txs_seqs_pos, FUN = unlist))
    txs_seqs_neg_cat <- DNAStringSet(lapply(X = txs_seqs_neg, FUN = function(x){
      unlist(rev(x))
    }))

    # merge sequences from positive and negative strand and sort according to name
    txs_seqs <- c(txs_seqs_pos_cat, txs_seqs_neg_cat)
    txs_seqs <- txs_seqs[order(names(txs_seqs))]

  }else{
    txs_seqs <- DNAStringSet()
  }

  # genome sequence --------------------------------------------------------------------------------

  if (include_genome == TRUE){
    message("Obtaining genome sequence...")

    # get chromosome names for which sequences should be extracted
    if (!is.null(exclude_chr_pattern)) {
      chroms <- grep(names(genome), pattern = exclude_chr_pattern, invert = TRUE, value = TRUE)
    }else{
      chroms <- names(genome)
    }

    # get sequences
    genome_seq <- getSeq(genome, names = chroms)

    # add chromosome names to sequences
    chrom_names <- paste(chroms, "chromosome", sub(chroms, pattern = "chr", replacement = ""))
    names(genome_seq) <- chrom_names

  }else{
    genome_seq <- DNAStringSet()
  }

  # finalize output --------------------------------------------------------------------------------

  # combine genome and transcriptome sequences
  gt_seqs <- c(genome_seq, txs_seqs)

  # write sequences to output fasta file
  if (length(gt_seqs) > 0) {
    message("Writing to output...")
    Biostrings::writeXStringSet(gt_seqs, filepath = output_fasta, format = "fasta",
                                compress = compress)
  }else{
    stop("Nothing to write to output!", call. = FALSE)
  }

}

#' Create BLAST database
#'
#' Create a database which can be used to estimate potential off-target priming of TAP-seq primers.
#'
#' @param blastdb A file path to the blastdb that should be created.
#' @param fasta A fasta file containing all sequences to be included in the database. Can be created
#'   with \code{\link[TAPseq]{createFasta}}.
#' @param makeblastdb Path (character) to the \code{makeblastdb} executable. Usually this is
#'   inferred when loading/attaching the package.
#' @param compression What compression was used when generating the fasta file? Default is "auto",
#'   which tries to infer compression from fasta filename. Currently only allows for gzip
#'   compression (.gz file extension).
#' @param title Optional title for BLAST database.
#' @export
createBLASTDb <- function(blastdb, fasta, makeblastdb = getOption("TAPseq.makeblastdb"),
                          compression = c("auto", "gzip", "none"), title = "None") {

  # get compression argument
  compression <- match.arg(compression)

  # infer compression from fasta file extension if compression is set to "auto"
  if (compression == "auto"){
    fasta_ext <- tools::file_ext(fasta)
    if (fasta_ext == "gz") compression <- "gzip"
    else compression <- "none"
  }

  # create blast data base based on fasta input compression
  if (compression == "gzip"){

    # parameters for blast data base generation
    db_args <- paste("-in - -input_type fasta -dbtype nucl -title", title, "-out", blastdb)

    # create full command
    cmd <- paste("gzip -dc", fasta, "|", makeblastdb, db_args)

    # run command
    system(cmd)

  }else{

    # parameters for blast data base generation
    db_args <- paste("-in", fasta, "-input_type fasta -dbtype nucl -title", title, "-out", blastdb)

    # create blast database
    system2(command = makeblastdb, args = db_args)

  }
}

#' BLAST TAP-seq primers
#'
#' Use BLAST to align designed TAP-seq primers against a genome and transcriptome database to
#' estimate off-target priming potential. Only hits where at least a specified portion of the
#' sequence involving the 3' end of the primer aligns with not more than a certain number of
#' mismatches are considered.
#'
#' @param object A \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object containing
#'   designed primers.
#' @param blastdb BLAST database to which primers should be aligned. Can be created with
#'   \code{\link[TAPseq]{createBLASTDb}}.
#' @param annot A \code{\link[GenomicRanges]{GRanges}} object containing exons of transcripts, which
#'   will be used to annotate blast hits on the genome sequence. The same object as used in
#'   \code{\link[TAPseq]{createFasta}}.
#' @param max_mismatch Maximum number of mismatches allowed for off-target hits (default: 0).
#' @param min_aligned Minimum portion of the primer sequence starting from the 3' end that must
#'   align for off-target hits (default: 0.75).
#' @param tmpdir Directory needed to store temporary files.
#' @param blastn Path (character) to the \code{blastn} executable. Usually this is inferred when
#'   loading/attaching the package.
#' @return A \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object with the number
#'   of potential off-target priming hits added to the TAP-seq primer metadata.
#' @export
setGeneric("blastPrimers",
           function(object, blastdb, annot, max_mismatch = 0, min_aligned = 0.75,
                    tmpdir = tempdir(), blastn = getOption("TAPseq.blastn"))
             standardGeneric("blastPrimers")
)

#' @describeIn blastPrimers Blast primers in a \code{TsIO} object
#' @export
setMethod("blastPrimers", "TsIO", function(object, blastdb, annot, max_mismatch, min_aligned,
                                           tmpdir, blastn) {

  # extract designed tap-seq primers
  primers <- tapseq_primers(object)

  if (length(primers) > 0) {

    # blast primers and process hits to count potential off-target primer binding
    primers <- blast_primers(primers, blastdb = blastdb, annot = annot, max_mismatch = max_mismatch,
                             min_aligned = min_aligned, tmpdir = tmpdir, blastn = blastn)

    # add annotated primers back to input object
    tapseq_primers(object) <- primers
    return(object)

  }else{
    message("No primers found in TsIO object!")
    return(object)
  }

})

#' @describeIn blastPrimers Blast primers in a \code{TsIOList} object
#' @export
setMethod("blastPrimers", "TsIOList", function(object, blastdb, annot, max_mismatch, min_aligned,
                                               tmpdir, blastn) {

  # extract designed tap-seq primers
  primers <- tapseq_primers(object)

  if (length(primers) > 0) {

    # blast primers and process hits to count potential off-target primer binding
    primers <- blast_primers(primers, blastdb = blastdb, annot = annot, max_mismatch = max_mismatch,
                             min_aligned = min_aligned, tmpdir = tmpdir, blastn = blastn)

    # split primers by target
    targets <- sub("\\.primer_left_\\d+$", "", names(primers))
    primers_split <- S4Vectors::split(primers, f = targets)

    # add primers to every TsIO object
    mendoapply(FUN = `tapseq_primers<-`, object, primers_split)

  }else{
    message("No primers found in TsIOList object!")
    return(object)
  }

})


# HELPER FUNCTIONS =================================================================================

# blast primers in IRanges format, process hits and add off-target hits to metadata
blast_primers <- function(primers, blastdb, annot, max_mismatch, min_aligned, tmpdir, blastn) {

  # extract primer sequence and store in DNAStringSet
  primer_seqs <- structure(mcols(primers)$sequence, names = names(primers))
  primer_seqs <- DNAStringSet(primer_seqs)

  # blast primers and get hits
  message("Running BLAST...")
  hits <- blast_primer_seqs(primer_seqs, blastdb = blastdb, tmpdir = tmpdir, blastn = blastn)

  # process blast hits
  message("Processing hits...")
  hits <- process_blast_hits(hits, annot = annot, tx_hits_as_exonic = TRUE)

  # filter hits for aligned sequence portion and max mismatches
  hits <- filter_3p_hits(hits, max_mismatch = max_mismatch, min_aligned = min_aligned)

  # count the number of off-target hits per primer
  off_hits <- count_off_targets(hits)

  # get existing primer meta data and remove any blast hits if there were already some.
  primer_meta <- mcols(primers)
  primer_meta <- primer_meta[, !colnames(primer_meta) %in% colnames(off_hits)]

  # add 0 in off_hits for weird cases where no hits were found by blast
  off_hits <- off_hits[rownames(primer_meta), ]
  off_hits[is.na(off_hits)] <- 0

  # add off-target hits to primer meta data
  mcols(primers) <- cbind(primer_meta, off_hits)
  return(primers)

}

# blast primers against an existing blast database
blast_primer_seqs <- function(primer_seqs, blastdb, tmpdir, blastn) {

  # save to temporary fasta file
  primers_fasta <- file.path(tmpdir, "primer_seqs.fasta")
  Biostrings::writeXStringSet(primer_seqs, filepath = primers_fasta)

  # set blastn parameters
  blast_out <- file.path(tmpdir, "blast_out.txt")
  blast_args <- paste("-query", primers_fasta, "-out", blast_out, "-db", blastdb, "-task blastn-short",
                      "-outfmt '6 qaccver qlen saccver stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore'")

  # run blastn for short query sequences
  system2(command = blastn, args = blast_args)

  # read blast ouptut
  hits <- read.table(blast_out, stringsAsFactor = FALSE, sep = "\t", quote = "")

  # set colnames
  colnames(hits) <- c("qaccver", "qlen", "saccver", "stitle", "pident", "length", "mismatch",
                      "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  # delete temporary fasta and blast output files
  unlink(c(primers_fasta, blast_out))

  # return found hits
  return(hits)

}

# process blast output and overlap genomic hits with annotations to get their target. annot has to
# be GRanges object containing the exon annotations of all genes that should be used for the
# annotation of genomic hits. tx_hits_as_exonic defines whether blast hits in transcripts should be
# named "exonic" (TRUE) or "transcript".
process_blast_hits <- function(hits, annot, tx_hits_as_exonic = TRUE) {

  # parse stitle field -----------------------------------------------------------------------------

  # split stitle into type and name
  stitle_split <- strsplit(hits$stitle, split = " ")
  types <- vapply(X = stitle_split, FUN = "[[", 2, FUN.VALUE = character(1))
  names <- vapply(X = stitle_split, FUN = "[[", 3, FUN.VALUE = character(1))

  # add types and names to hits (replace stitle)
  hits <- cbind.data.frame(hits[, 1:3], "name" = names, "type" = types, hits[, 5:ncol(hits)],
                           stringsAsFactors = FALSE)

  # prepare annotation data ------------------------------------------------------------------------

  # set gene_names as element names for annot
  names(annot) <- annot$gene_name

  # get gene locus coordinates (1st to last exon) of each gene
  genes <- split(annot, f = annot$gene_name)  # get all exons per gene
  genes <- unlist(range(genes))

  # process genomic hits ---------------------------------------------------------------------------

  # split hits into genomic and transcript his
  thits <- hits[hits$type == "transcript", ]
  ghits <- hits[hits$type == "chromosome", ]

  # get chromosome ids
  chr <- ghits$saccver

  # get start and end of primer binding sites
  start_end <- as.matrix(ghits[, c("sstart", "send")])

  # order start and stop for each site (i.e. each row)
  start_end <- t(apply(X = start_end, MARGIN = 1, FUN = sort))

  # create GRanges with all primer binding sites
  hit_coords <- IRanges(start = start_end[, 1], end = start_end[, 2])
  hit_coords <- GRanges(seqnames = chr, ranges = hit_coords)

  # get gene loci and exons overlapping the primer binding sites
  gene_overlaps <- find_hit_overlaps(hit_coords = hit_coords, subjects = genes)
  exon_overlaps <- find_hit_overlaps(hit_coords = hit_coords, subjects = annot)

  # add gene overlaps as name to genomic hits (i.e. the gene lous a given hit
  # overlaps with)
  ghits$name <- gene_overlaps

  # set type according to overlaps found: exon, intron or intergenic
  ghits$type[!is.na(exon_overlaps)] <- "exonic"
  ghits$type[!is.na(gene_overlaps) & is.na(exon_overlaps)] <- "intronic"
  ghits$type[ghits$type == "chromosome"] <- "intergenic"  # remaining hits

  # combine genomic and transcript hits ------------------------------------------------------------

  # change type of transcript hits to "exonic" if specified
  if (tx_hits_as_exonic == "TRUE") thits$type <- "exonic"

  # combine genome and transcriptome hits
  output <- rbind.data.frame(thits, ghits)

  # order according to primer id
  output[order(output$qaccver), ]

}

# function to find overlapping subjects for each hit
find_hit_overlaps <- function(hit_coords, subjects) {

  # find overlaps of of primer binding sites with subjects
  overlaps <- GenomicRanges::findOverlaps(query = hit_coords, subject = subjects)

  # create data.frame containing the subject name for each overlap
  overlaps <- data.frame("primer" = from(overlaps), "subject" = names(subjects[to(overlaps)]),
                         stringsAsFactors = FALSE)

  # split overlapping subjects by primer index
  overlaps_split <- split(overlaps$subject, f = overlaps$primer)

  # paste unique overlapping subjects for each hit together, if multiple are found
  overlaps_split <- lapply(X = overlaps_split, FUN = unique)
  overlaps <- sapply(X = overlaps_split, FUN = paste, collapse = ";")

  # create vector for output containing NA for each primer hit
  output <- rep(NA, times = length(hit_coords))

  # fill in any overlapping subjects for each primer hit
  output[as.numeric(names(overlaps))] <- overlaps
  return(output)

}

# filter hits for hits involving the 3' end of primers. max_mismatch specifies the number of
# mismatched bases allowed in alignments to be considered. min_aligned (0-1) specifies the minimal
# fraction of the primer sequence that needs to align to the hit. (e.g. max_mismatch = 1 and
# min_aligned = 0.75 means that 75% of the the primer sequence is required to align to the target
# with a maximum of 1 mismatch.
filter_3p_hits <- function(hits, max_mismatch = 0, min_aligned = 0.75) {

  # remove hits with more mismatches than max_mismatch
  hits <- hits[hits$mismatch <= max_mismatch, ]

  # remove hits that cover less than min_match_perc of the primer sequence
  min_length <- round(hits$qlen * min_aligned)
  hits <- hits[hits$length >= min_length, ]

  # retain hits around the 3' end of the primer (end of the primer = qlen)
  hits[hits$qend == hits$qlen, ]

}

# count the number of off-targets per primer and off-target type (exonic, intronic, intergenic).
count_off_targets <- function(hits) {

  # get hits in genic regions and only retain primer id, type and taget name
  genic_hits <- hits[hits$type != "intergenic", c("qaccver", "name", "type")]
  genic_hits <- unique(genic_hits)

  # get intergenic hits
  intergenic_hits <- hits[hits$type == "intergenic", c("qaccver", "name", "type")]

  # combine unique genic and intergenic hits and convert query id and type to factors
  unique_hits <- rbind.data.frame(genic_hits, intergenic_hits)
  unique_hits$qaccver <- as.factor(unique_hits$qaccver)
  unique_hits$type <- factor(unique_hits$type, levels = c("intergenic", "intronic", "exonic"))

  # define on- and off-target hits based on intended primer target from primer id
  unique_hits$target <- sub(unique_hits$qaccver, pattern = "\\..*", replacement = "")
  unique_hits$on_target <- unique_hits$name == unique_hits$target
  unique_hits[unique_hits$type == "intergenic", "on_target"] <- FALSE

  # only retain off-target hits
  off_target_hits <- unique_hits[unique_hits$on_target == FALSE, ]

  # split hits by primer id
  off_target_hits <- split(off_target_hits, f = off_target_hits$qaccver, drop = FALSE)

  # count number of off-target hits per primer
  n_off_hits <- lapply(off_target_hits, FUN = function(x) {
    off_counts <- table(x$type)
    t(as.matrix(off_counts))
  })

  # convert into one data.frame
  output <- data.frame(do.call(rbind, n_off_hits), row.names = names(n_off_hits))
  colnames(output) <- paste0(colnames(output), "_off_targets")
  return(output)

}
