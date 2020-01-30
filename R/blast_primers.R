#' Create BLAST database
#'
#' Create a BLAST database containing genome and transcriptome sequences. Used by
#' \code{\link[TAPseq]{blastPrimers}} to estimate potential off-target priming of TAP-seq primers.
#' The created database contains both sequence files for BLAST and annotations to process the
#' results.
#'
#' @param genome A \code{\link[BSgenome]{BSgenome}} (or \code{\link[Biostrings]{DNAStringSet}})
#'   object containing the sequences of all chromosomes to obtain genome and transcript sequences.
#' @param annot A \code{\link[GenomicRanges]{GRanges}} object containing all exons of transcripts to
#'   be considered.
#' @param blastdb Path to the directory where the database should be created. Will be created if not
#'   existing. If the directory already exist, the function will raise a warning, but overwrite any
#'   previous blast database files (other files stay untouched).
#' @param standard_chromosomes (logical) Specifies whether only standard chromosomes should be
#'   included in output genome sequences (e.g. chr1-22, chrX, chrY, chrM for homo sapiens).
#' @param tx_id,tx_name,gene_name,gene_id (character) Column names in annot metadata containing
#'   transcript id, transcript name, gene name and gene id information.
#' @param verbose (logical) If \code{TRUE}, additional information from \code{makeblastdb} is
#'   printed to the console. Default: \code{FALSE}.
#' @param makeblastdb Path to the \code{makeblastdb} executable. Usually this is inferred when
#'   loading/attaching the package.
#' @param title Optional title for BLAST database.
#' @examples
#' \dontrun{
#' library(TAPseq)
#' library(BSgenome)
#'
#' # human genome (hg38) BSgenome object
#' hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
#'
#' # get annotations for BLAST
#' annot_url <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
#'                     "gencode.v32.annotation.gtf.gz")
#' annot <- import(annot_url, format = "gtf")
#' blast_exons <- annot[annot$type == "exon" & annot$gene_type == "protein_coding"]
#'
#' # build BLAST database
#' blastdb <- file.path(tempdir(), "blastdb")
#' createBLASTDb(genome = hg38, annot = blast_exons, blastdb = blastdb)
#' }
#' @export
createBLASTDb <- function(genome, annot, blastdb, standard_chromosomes = TRUE,
                          tx_id = "transcript_id", tx_name = "transcript_name",
                          gene_name = "gene_name", gene_id = "gene_id", title = "TAP-seq_GT_DB",
                          verbose = FALSE, makeblastdb = getOption("TAPseq.makeblastdb")) {

  # make sure that genome is a BSgenome object or DNAStringSet object
  if (!any(is(genome, "BSgenome"), is(genome, "DNAStringSet"))) {
    stop("genome must be of class BSgenome or DNAStringSet!", call. = FALSE)
  }

  # create blast database directory if required
  dir.create(blastdb, recursive = TRUE, showWarnings = TRUE)

  # basename for all files in database
  db_base <- file.path(blastdb, "gt_blast")

  # process annotations ----------------------------------------------------------------------------

  # make sure that annot is a GRanges object
  message("Processing annotations...")
  if (is(annot, "GRangesList")) {
    warning("Coercing annot from 'GRangesList' to 'GRanges' (i.e. unlist(annot)).", call. = FALSE)
    annot <- unlist(annot)
  } else if (!is(annot, "GRanges")) {
    stop("annot needs to be a 'GRanges' object.", call. = FALSE)
  }

  # check if any specified identifiers are missing from annotation metadata
  ids <- c(tx_id, tx_name, gene_id, gene_name)
  missing_ids <- !ids %in% colnames(mcols(annot))
  if (any(missing_ids)) {
    stop("'", paste(ids[missing_ids], collapse = ", "), "' not found in annot metadata!",
         call. = FALSE)
  }

  # process annotations for tap-seq blast database
  mcols(annot) <- mcols(annot)[ids]
  colnames(mcols(annot)) <- c("transcript_id", "transcript_name", "gene_id", "gene_name")

  # save annotations to blast database directory
  saveRDS(annot, file = paste0(db_base, ".annot.rds"), compress = TRUE)

  # create genome and transcript sequences ---------------------------------------------------------

  # extract genome and transcriptome sequences
  gt_seqs <- get_gt_sequences(genome = genome, annot = annot,
                             standard_chromosomes = standard_chromosomes)

  # save to fasta file in blast database directory
  message("Writing to fasta file...")
  seq_fasta <- tempfile(pattern = "gt_seqs_", tmpdir = blastdb, fileext = ".fasta")
  Biostrings::writeXStringSet(gt_seqs, filepath = seq_fasta, format = "fasta", compress = FALSE)

  # create blast sequence database -----------------------------------------------------------------

  # process verbose argument
  if (verbose == TRUE) {
    stderr_print <- ""
  } else {
    stderr_print <- FALSE
  }

  # parameters for blast data base generation
  db_arg <- paste("-in", seq_fasta, "-input_type fasta -dbtype nucl -title", title, "-out", db_base)

  # create blast database
  message("Creating BLAST database...")
  system2(command = makeblastdb, args = db_arg, stdout = stderr_print)

  # delete temporary fasta file
  unlink(seq_fasta)
  message("Done!")

}

#' BLAST TAP-seq primers
#'
#' Use BLAST to align designed TAP-seq primers against a genome and transcriptome database to
#' estimate off-target priming potential. Only hits where at least a specified portion of the
#' sequence involving the 3' end of the primer aligns with not more than a certain number of
#' mismatches are considered.
#'
#' \code{blastPrimers} counts the number of genes in which a primer has 1) exonic hits or 2)
#' intronic hits, or 3) the number of hits in intergenic regions of the genome. The exonic and
#' intronic counts should be interptreted as: "In how many genes does a primer have exonic
#' (or intronic) hits?".
#'
#' If a BLAST hit falls in both intronic and exonic regions of a given gene (i.e. exonic for one
#' transcript, intronic for another transcript), only the exonic hit is counted for that gene. If a
#' primer has for instance 3 BLAST hits in one gene, 2 exonic and 1 intronic, then one exonic hit
#' and one intronic hit is counted for that gene.
#'
#' If seqID of the designed primers (\code{\link[TAPseq]{sequence_id}}) refer to the target
#' gene/transcripts and can be found in the BLAST database annotations via \code{primer_targets},
#' then only off-target hits are counted. This is usually the case if input for primer design was
#' produced from target gene annotations.
#'
#' @param object A \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object containing
#'   designed primers.
#' @param blastdb TAP-seq BLAST database created with \code{\link[TAPseq]{createBLASTDb}}.
#' @param max_mismatch Maximum number of mismatches allowed for off-target hits (default: 0).
#' @param min_aligned Minimum portion of the primer sequence starting from the 3' end that must
#'   align for off-target hits (default: 0.75).
#' @param primer_targets Specifies what should be used to identify primer targets for off-target
#'   identification. I.e. to what does the 'seqID' in TsIO objects refer? Can be a subset of
#'   \code{transcript_id}, \code{transcript_name}, \code{gene_id} or \code{gene_name}. By default
#'   all 4 are checked. Set to \code{NULL} to disable any off-target identification. See Details for
#'   more information.
#' @param tmpdir Directory needed to store temporary files.
#' @param blastn Path (character) to the \code{blastn} executable. Usually this is inferred when
#'   loading/attaching the package.
#' @return A \code{\link[TAPseq]{TsIO}} or \code{\link[TAPseq]{TsIOList}} object with the number
#'   of potential off-target priming hits added to the TAP-seq primer metadata.
#' @examples
#' \dontrun{
#' library(TAPseq)
#' library(BSgenome)
#'
#' # human genome (hg38) BSgenome object
#' hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
#'
#' # get annotations for BLAST
#' annot_url <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
#'                     "gencode.v32.annotation.gtf.gz")
#' annot <- import(annot_url, format = "gtf")
#' blast_exons <- annot[annot$type == "exon" & annot$gene_type == "protein_coding"]
#'
#' # build BLAST database
#' blastdb <- file.path(tempdir(), "blastdb")
#' createBLASTDb(genome = hg38, annot = blast_exons, blastdb = blastdb)
#'
#' # chr11 primers example data (already contains off-targets, but we can overwrite them)
#' data("chr11_primers")
#' chr11_primers <- chr11_primers[1:3]  # only use a small subset for this example
#'
#' # run blast to identify potential off-targets
#' chr11_primers <- blastPrimers(chr11_primers, blastdb = blastdb)
#' tapseq_primers(chr11_primers)
#'
#' # allow 1 mismatch between primer and off-target
#' chr11_primers <- blastPrimers(chr11_primers, blastdb = blastdb, max_mismatch = 1)
#' tapseq_primers(chr11_primers)
#' }
#' @export
setGeneric("blastPrimers",
           function(object, blastdb, max_mismatch = 0, min_aligned = 0.75,
                    primer_targets = c("transcript_id", "transcript_name", "gene_id", "gene_name"),
                    tmpdir = tempdir(), blastn = getOption("TAPseq.blastn"))
             standardGeneric("blastPrimers")
)

#' @describeIn blastPrimers BLAST primers in a \code{TsIO} object
#' @export
setMethod("blastPrimers", "TsIO", function(object, blastdb, max_mismatch, min_aligned,
                                           primer_targets, tmpdir, blastn) {

  # check primer_targets argument for valid choices
  if (!is.null(primer_targets)) {
    primer_targets <- match.arg(primer_targets, several.ok = TRUE)
  }

  # extract designed tap-seq primers
  primers <- tapseq_primers(object)

  if (length(primers) > 0) {

    # blast primers and process hits to count potential off-target primer binding
    primers <- blast_primers(primers, blastdb = blastdb, max_mismatch = max_mismatch,
                             min_aligned = min_aligned, primer_targets = primer_targets,
                             tmpdir = tmpdir, blastn = blastn)

    # add annotated primers back to input object
    tapseq_primers(object) <- primers
    return(object)

  }else{
    message("No primers found in TsIO object!")
    return(object)
  }

})

#' @describeIn blastPrimers BLAST primers in a \code{TsIOList} object
#' @export
setMethod("blastPrimers", "TsIOList", function(object, blastdb, max_mismatch, min_aligned,
                                               primer_targets, tmpdir, blastn) {

  # check primer_targets argument for valid choices
  if (!is.null(primer_targets)) {
    primer_targets <- match.arg(primer_targets, several.ok = TRUE)
  }

  # extract designed tap-seq primers
  primers <- tapseq_primers(object)

  if (length(primers) > 0) {

    # blast primers and process hits to count potential off-target primer binding
    primers <- blast_primers(primers, blastdb = blastdb, max_mismatch = max_mismatch,
                             min_aligned = min_aligned, primer_targets = primer_targets,
                             tmpdir = tmpdir, blastn = blastn)

    # split primers by target
    targets <- sub("\\.primer_left_\\d+$", "", names(primers))
    primers_split <- S4Vectors::split(primers, f = targets)

    # add empty IRanges objects for TsIO objects without any primers
    missing <- setdiff(names(object), names(primers_split))
    missing_ranges <- IRangesList(lapply(missing, FUN = function(x) IRanges() ))
    names(missing_ranges) <- missing
    primers_split <- c(primers_split, missing_ranges)[names(object)]

    # add primers to every TsIO object
    mendoapply(FUN = `tapseq_primers<-`, object, primers_split)

  }else{
    message("No primers found in TsIOList object!")
    return(object)
  }

})

# HELPER FUNCTIONS =================================================================================

#' Get genome and transcriptome sequences
#'
#' Get DNA sequences of all chromosomes and all annotated transcripts of a genome. This function is
#' used to create the sequences in \code{\link[TAPseq]{createBLASTDb}}.
#'
#' @param genome A \code{\link[BSgenome]{BSgenome}} (or \code{\link[Biostrings]{DNAStringSet}})
#'   object containing the chromosome sequences to obtain genome and / or transcript sequences.
#' @param annot A \code{\link[GenomicRanges]{GRanges}} object containing all exons of transcripts to
#'   be considered. If not specified, no transcript sequences will be included in the output fasta
#'   file.
#' @param tx_id,tx_name,gene_name,gene_id (character) Column names in annot metadata containing
#'   transcript id, transcript name, gene name and gene id information. These column are mandatory,
#'   but can contain internal names (e.g. "transcript-1" or "gene-1").
#' @param include_genome (logical) Specifies whether the genome sequence should be included in the
#'   output fasta file.
#' @param standard_chromosomes (logical) Specifies whether only standard chromosomes should be
#'   included in output genome sequences (e.g. chr1-22, chrX, chrY, chrM for homo sapiens).
#' @param compress (logical) Create a gzipped output fasta file.
#' @return A \code{\link[Biostrings]{DNAStringSet}} object containing the genome and transcriptome
#'   sequences.
#' @keywords internal
get_gt_sequences <- function(genome, annot = NULL, tx_id = "transcript_id",
                            tx_name = "transcript_name", gene_name = "gene_name",
                            gene_id = "gene_id", include_genome = TRUE,
                            standard_chromosomes = TRUE) {

  # transcript sequences ---------------------------------------------------------------------------

  if (!is.null(annot)) {
    message("Obtaining transcript sequences...")

    # check if any specified identifiers are not found in annotation metadata
    metadat <- mcols(annot)
    ids <- c(tx_id, tx_name, gene_id, gene_name)
    req_ids <- ids[!is.na(ids)]
    missing_ids <- !req_ids %in% colnames(metadat)
    if (any(missing_ids)) {
      stop("'", paste(req_ids[missing_ids], collapse = ", "), "' not found in annot metadata!",
           call. = FALSE)
    }

    # create a new identifier for each exon consisting of the transcript id & name, and gene id &
    # name. This id will later be used as the title of each entry in the fasta file.
    tx_ids <- paste("lcl|transcript", metadat[[tx_id]], metadat[[tx_name]], metadat[[gene_id]],
                    metadat[[gene_name]], sep = ";")

    # split exons by this new tx_ids into GRangesList with each element
    # containing all exons of a given transcript
    txs <- split(annot, f = tx_ids)

    # get sequences of all transcripts
    txs_seqs <- getTxsSeq(txs, genome = genome)

  }else{
    txs_seqs <- DNAStringSet()
  }

  # genome sequence --------------------------------------------------------------------------------

  if (include_genome == TRUE){
    message("Obtaining genome sequence...")

    # get chromosome names for which sequences should be extracted
    if (standard_chromosomes ==  TRUE) {
      chroms <- GenomeInfoDb::standardChromosomes(genome)
    }else{
      chroms <- names(genome)
    }

    # get sequences
    genome_seq <- getSeq(genome, names = chroms)

    # add "sequence type" to names
    names(genome_seq) <- paste("lcl|chromosome", names(genome_seq), sep = ";")

  }else{
    genome_seq <- DNAStringSet()
  }

  # combine and return genome and transcriptome sequences
  c(genome_seq, txs_seqs)

}

# blast primers in IRanges format, process hits and add off-target hits to metadata
blast_primers <- function(primers, blastdb, max_mismatch, min_aligned, primer_targets, tmpdir,
                          blastn) {

  # load annot from blast database
  annot <- readRDS(file.path(blastdb, "gt_blast.annot.rds"))

  # extract primer sequence and store in DNAStringSet
  primer_seqs <- structure(mcols(primers)[["sequence"]], names = names(primers))
  primer_seqs <- DNAStringSet(primer_seqs)

  # blast primers and get hits
  message("Running BLAST...")
  hits <- blast_primer_seqs(primer_seqs, blastdb = blastdb, tmpdir = tmpdir, blastn = blastn)

  # annotate blast hits
  message("Processing hits...")
  hits <- annotate_blast_hits(hits, annot = annot)

  # filter for primer 3' end hits with minimum aligned sequence portion and max mismatches
  hits <- filter_3p_hits(hits, max_mismatch = max_mismatch, min_aligned = min_aligned)

  # find primer target genes and filter for off-target hits
  if (!is.null(primer_targets)) {
    hits <- get_primer_target_genes(hits, annot = annot, primer_targets = primer_targets)
    off_target_filter <- hits$gene_id != hits$target_gene_id | is.na(hits$gene_id)
    hits <- hits[off_target_filter, ]
  }

  # count the number of hits per primer
  hit_counts <- count_primer_hits(hits)

  # get existing primer meta data and remove any blast hits if there were already some added
  primer_meta <- mcols(primers)
  primer_meta <- primer_meta[, !colnames(primer_meta) %in% colnames(hit_counts)]

  # add 0 in hit_counts for cases where no hits are found (e.g. no off-targets)
  hit_counts <- hit_counts[rownames(primer_meta), ]
  hit_counts[is.na(hit_counts)] <- 0

  # add off-target hits to primer meta data
  mcols(primers) <- cbind(primer_meta, hit_counts)
  return(primers)

}

# blast primers against an existing blast database
blast_primer_seqs <- function(primer_seqs, blastdb, tmpdir, blastn) {

  # save to temporary fasta file
  primers_fasta <- file.path(tmpdir, "primer_seqs.fasta")
  Biostrings::writeXStringSet(primer_seqs, filepath = primers_fasta)

  # create blastn arguments
  blastdb_bn <- file.path(blastdb, "gt_blast")
  outfmt <- c("qaccver", "qlen", "saccver", "length", "mismatch", "qend", "sstart", "send")
  blast_args <- paste0("-query ", primers_fasta, " -db ", blastdb_bn, " -task blastn-short ",
                       "-outfmt '6 ", paste(outfmt, collapse = " "), "'")

  # run blastn for short query sequences
  hits <- system2(command = blastn, args = blast_args, stdout = TRUE)

  # process output text into data.frame
  hits <- utils::read.table(text = hits, col.names = outfmt, sep = "\t", stringsAsFactor = FALSE)

  # delete temporary fasta
  unlink(primers_fasta)

  # return found hits
  return(hits)

}

# process blast output and overlap genomic hits with annotations to get their target. annot contains
# exons of all genes used to annotate hits and must contain the metadata columns gene_id and
# gene_name. typically the annot file stored in a blast database created by createBLASTDb is used.
annotate_blast_hits <- function(hits, annot) {

  # add hit identifier for every hit
  hits$hit_id <- paste0("hit_", seq_len(nrow(hits)))

  # split saccver into type and name
  hits$saccver <- sub("lcl\\|", "", hits$saccver)
  hits <- separate(hits, col = "saccver", sep = ";", fill = "right", remove = TRUE,
                   into = c("subject_type", "id", "name", "gene_id", "gene_name"))

  # set gene_ids as names for annot
  gene_identifiers <- as.data.frame(mcols(annot)[, c("gene_id", "gene_name")])
  names(annot) <- gene_identifiers$gene_id

  # get merged exons for every gene and gene locus coordinates
  exons <- split(annot, f = names(annot))
  exons_red <- reduce(exons)
  genes <- range(exons_red)

  # split hits into genomic and transcript his
  ghits <- hits[hits$subject_type == "chromosome", ]
  thits <- hits[hits$subject_type == "transcript", ]

  # process genomic hits ---------------------------------------------------------------------------

  # create GRanges object with primer binding sites for all hits in ghits
  chr <- ghits$id
  primer_site <- as.matrix(ghits[, c("sstart", "send")])
  primer_site <- t(apply(X = primer_site, MARGIN = 1, FUN = sort))  # order starts and stops
  hit_coords <- GRanges(seqnames = chr,
                        ranges = IRanges(start = primer_site[, 1], end = primer_site[, 2]))
  names(hit_coords) <- ghits$hit_id

  # get gene loci and exons overlapping the primer binding sites
  gene_overlaps <- find_hit_overlaps(hit_coords = hit_coords, subjects = genes)
  exon_overlaps <- find_hit_overlaps(hit_coords = hit_coords, subjects = exons_red)
  overlaps <- bind_rows(gene = gene_overlaps, exon = exon_overlaps, .id = "type")

  # convert to wide format to add hit type (intronic or exonic)
  overlaps <- pivot_wider(overlaps, names_from = "type", values_from = "hit",
                          values_fill = list(hit = 0))
  overlaps$type <- c("intronic", "exonic")[rowSums(overlaps[, c("gene", "exon")])]
  overlaps <- overlaps[, c("hit_id", "gene_id", "type")]

  # add gene name for every gene id
  overlaps <- left_join(x = overlaps, y = distinct(gene_identifiers), by = "gene_id")

  # merge with data in ghits
  ghits <- ghits[, !colnames(ghits) %in% c("gene_id", "gene_name")]
  ghits_processed <- left_join(x = ghits, y = overlaps, by = "hit_id")

  # set type to intergenic for all hits that did not overlap any exons or gene loci
  ghits_processed$type[is.na(ghits_processed$type)] <- "intergenic"

  # process transcript hits ------------------------------------------------------------------------

  # add type to thits and extract the same columns as for genomic hits
  thits$type <- "exonic"
  thits_processed <- thits[, colnames(ghits_processed)]

  # combined processed ghits and thits
  bind_rows(ghits_processed, thits_processed)

}

# function to find overlapping subjects for each hit
find_hit_overlaps <- function(hit_coords, subjects) {

  # find overlaps of of primer binding sites with subjects
  overlaps <- findOverlaps(query = hit_coords, subject = subjects)

  # create data.frame containing the subject name for each overlap
  data.frame(hit_id = names(hit_coords)[queryHits(overlaps)], hit = 1,
             gene_id = names(subjects)[subjectHits(overlaps)], stringsAsFactors = FALSE)
}

# filter hits for hits involving the 3' end of primers. max_mismatch specifies the number of
# mismatched bases allowed in alignments to be considered. min_aligned (0-1) specifies the minimal
# fraction of the primer sequence that needs to align to the hit. (e.g. max_mismatch = 1 and
# min_aligned = 0.75 means that 75% of the the primer sequence is required to align to the target
# with a maximum of 1 mismatch.
filter_3p_hits <- function(hits, max_mismatch, min_aligned) {

  # remove hits with more mismatches than max_mismatch
  hits <- hits[hits$mismatch <= max_mismatch, ]

  # remove hits that cover less than min_match_perc of the primer sequence
  min_length <- round(hits$qlen * min_aligned)
  hits <- hits[hits$length >= min_length, ]

  # retain hits around 3' end of primers (end of query (primer) alignment == end of the primer)
  hits[hits$qend == hits$qlen, ]

}

# get the intended target gene for every primer in annotated blast hits. annot is the same object as
# used in annotate_blast_hits.
get_primer_target_genes <- function(hits, annot, primer_targets) {

  # get sequence id from primer id
  primers_ids <- unique(hits$qaccver)
  primer_seq_ids <- structure(sub(".primer_left_\\d+", "", primers_ids), names = primers_ids)

  # create data.frame containing transcript_ids, transcript_names, gene_ids and gene_names
  gene_names <- distinct(as.data.frame(mcols(annot)))
  gene_names$target_gene <- gene_names$gene_id

  # only retain gene_names data on specified primer_targets
  gene_names <- gene_names[, c(primer_targets, "target_gene")]

  # find primer targets in data provided in annot
  target_genes <- lapply(primer_seq_ids, FUN = find_primer_target_genes, gene_names = gene_names)

  # count number of inferred target_genes per primer
  n_match <- unlist(lapply(target_genes, FUN = length))
  if (!all(n_match == 1)) {
    warning("Can't find primer targets for all primers with primer_targets set to '",
            paste(primer_targets, collapse= ", "), "'!", call. = FALSE)
  }

  # create data.frame with target gene for every primer
  target_gene_ids <- lapply(target_genes[n_match == 1], FUN = function(x) {
    data.frame(target_gene_id = x, stringsAsFactors = FALSE)
  })
  target_gene_ids <- bind_rows(target_gene_ids, .id = "qaccver")

  # add primer target genes to hits
  if (nrow(target_gene_ids) > 0) {
    hits <- left_join(x = hits, y = target_gene_ids, by = "qaccver")
  } else {
    hits$target_gene_id <- as.character(NA)
  }

  # return hits with primer targets
  return(hits)

}

# find all target gene ids associated with a given primer_seq_id
find_primer_target_genes <- function(primer_seq_id, gene_names) {
  primer_match <- gene_names[, -5] == primer_seq_id
  gene_names_primer <- gene_names[rowSums(primer_match) > 0, ]
  unique(gene_names_primer$target_gene)
}

# count the number of blast hits
count_primer_hits <- function(hits) {

  # get unique hits in genic regions and only retain primer id, type and target name
  genic_hits <- hits[hits$type != "intergenic", c("qaccver", "gene_id", "type")]
  genic_hits <- distinct(genic_hits)

  # get intergenic hits
  intergenic_hits <- hits[hits$type == "intergenic", c("qaccver", "gene_id", "type")]

  # combine genic and intergenic hits and convert query id and type to factors
  unique_hits <- bind_rows(genic_hits, intergenic_hits)
  unique_hits$qaccver <- as.factor(unique_hits$qaccver)
  unique_hits$type <- factor(unique_hits$type, levels = c("intergenic", "intronic", "exonic"))

  # split hits by primer id
  unique_hits <- split(unique_hits, f = unique_hits$qaccver, drop = FALSE)

  # count number of off-target hits per primer
  n_hits <- lapply(unique_hits, FUN = function(x) {
    hit_counts <- table(x$type)
    as.data.frame(t(as.matrix(hit_counts)))
  })

  # convert into one data.frame
  output <- bind_rows(n_hits)
  rownames(output) <- names(n_hits)
  colnames(output) <- paste0(colnames(output), "_off_targets")
  return(output)

}
