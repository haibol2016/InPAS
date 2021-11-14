#' Helper function to label the last component of a genomic feature for each
#' transcript
#'
#' @param gr A tibble converted from an object of [GenomicRanges::GRanges-class]
#' @param feature_alt A character(1) vector, specifying the type of genomic
#'   features, such as "CDS", "exon", "utr3", "utr5".
#'
#' @importFrom dplyr as_tibble mutate filter arrange bind_rows group_by
#'   left_join summarise n
#' @author Haibo Liu
#' @return An object of [GenomicRanges::GRanges-class]
#' @keywords internal

assign_feature <- function(gr, feature_alt = "utr3") {
  gr <- gr %>%
    dplyr::arrange(
      seqnames, transcript,
      ifelse(strand == "+", -start, start)) %>%
    dplyr::mutate(feature = ifelse(!duplicated(paste0(seqnames, "|", transcript)),
    feature_alt, feature))
  gr
}

#' Extract gene models from a TxDb object
#'
#' Extract gene models from a TxDb object and annotate last 3' UTR exons and
#' the last CDSs
#'
#' @details A good practice is to perform read alignment using a
#' reference genome from Ensembl/GenCode including only the primary
#' assembly and build a TxDb using the GTF/GFF files downloaded
#' from the same source as the reference genome, such as
#' BioMart/Ensembl/GenCode. For instruction, see Vignette of the
#' GenomicFeatures. The UCSC reference genomes and their
#' annotation can be very cumbersome.
#'
#' @param TxDb An object of [GenomicFeatures::TxDb-class]
#' @param edb An object of [ensembldb::EnsDb-class]
#' @param removeScaffolds A logical(1) vector, whether the scaffolds should be
#'   removed from the genome If you use a TxDb containing alternative
#'   scaffolds, you'd better to remove the scaffolds.
#' @return A [GenomicRanges::GRanges-class] object for gene models
#'
#' @import GenomicFeatures
#' @importFrom plyranges as_granges complement_ranges disjoin_ranges filter
#'   group_by mutate reduce_ranges reduce_ranges_directed remove_names select
#'   set_genome_info shift_downstream summarise
#' @importFrom dplyr as_tibble mutate select pull filter arrange bind_rows
#'   group_by left_join summarise n
#' @importFrom magrittr %>%
#' @importFrom AnnotationDbi mapIds
#' @author Haibo Liu
#' @export
#'
#' @examples
#' library("EnsDb.Hsapiens.v86")
#' library("GenomicFeatures")
#' samplefile <- system.file("extdata",
#'                          "hg19_knownGene_sample.sqlite",
#'                           package = "GenomicFeatures")
#' TxDb <- loadDb(samplefile)
#' edb <- EnsDb.Hsapiens.v86
#'
#' parsed_Txdb <- parse_TxDb(TxDb, edb, 
#'                           removeScaffolds = TRUE)

parse_TxDb <- function(TxDb = NULL, edb = NULL,
                      removeScaffolds = FALSE) {
  if (is.null(TxDb) || !is(TxDb, "TxDb")) {
    stop("TxDb must be an object of TxDb!")
  }

  if (is.null(edb) || !is(edb, "EnsDb")) {
    stop("edb must be an object of EnsDb!")
  }

  if (removeScaffolds) {
    seqlevels(TxDb) <- seqlevels(TxDb)[grepl("^(chr)?(\\d+|[XY])$",
      seqlevels(TxDb),
      perl = TRUE
    )]
  }

  ## get all transcripts
  tx <- unlist(transcriptsBy(TxDb, by = "gene")) %>%
    plyranges::mutate(gene = names(.)) %>%
    plyranges::remove_names() %>%
    plyranges::select(-c("tx_id")) %>%
    data.frame() %>%
    as_tibble()

  ## recover CDS and UTRs for coding transcripts
  cds <- unlist(cdsBy(TxDb, by = "tx", use.names = TRUE)) %>%
    plyranges::mutate(transcript = names(.), feature = "CDS") %>%
    plyranges::remove_names() %>%
    plyranges::select(-c("cds_id", "cds_name")) %>%
    data.frame() %>%
    as_tibble()
  cds <- assign_feature(cds, feature_alt = "lastCDS")


  utr5 <- unlist(fiveUTRsByTranscript(TxDb,
    use.names = TRUE
  )) %>%
    plyranges::mutate(transcript = names(.), feature = "utr5") %>%
    plyranges::remove_names() %>%
    plyranges::select(-c("exon_id", "exon_name")) %>%
    data.frame() %>%
    as_tibble()

  utr3 <- unlist(threeUTRsByTranscript(TxDb,
    use.names = TRUE
  )) %>%
    plyranges::mutate(transcript = names(.), feature = "utr3") %>%
    plyranges::remove_names() %>%
    plyranges::select(-c("exon_id", "exon_name")) %>%
    data.frame() %>%
    as_tibble()
  utr3 <- assign_feature(utr3, feature_alt = "lastutr3")

  ## exons for non-coding transcripts
  exons <- unlist(exonsBy(TxDb,
    by = "tx",
    use.names = TRUE
  )) %>%
    plyranges::mutate(transcript = names(.)) %>%
    plyranges::remove_names() %>%
    plyranges::select(-c("exon_id", "exon_name"))

  noncoding_exons <- exons %>%
    plyranges::filter(!transcript %in% cds$transcript) %>%
    data.frame() %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(seqnames, transcript, -exon_rank)

  nexons <- noncoding_exons %>%
    dplyr::group_by(transcript) %>%
    summarise(num_exon = n())

  # noncoding RNA with single exon
  singleton <- nexons$transcript[nexons$num_exon == 1]
  noncoding_exons_singleton <- noncoding_exons %>%
    dplyr::filter(transcript %in% singleton) %>%
    dplyr::mutate(feature = "ncRNA")

  # noncoding RNA with mulitple exon
  noncoding_exons_multiplex <- noncoding_exons %>%
    dplyr::filter(!transcript %in% singleton)

  multiplex <- nexons[nexons$num_exon != 1, ] %>%
    as.data.frame()
  rownames(multiplex) <- multiplex[, 1]
  multiplex <- multiplex[, -1, drop = FALSE]
  multiplex <- 
    multiplex[unique(noncoding_exons_multiplex$transcript), , drop = F]
  feature <- do.call("c", lapply(multiplex[, 1], function(.x) {
    c("lastutr3", "lastCDS", rep("CDS", .x - 2))
  }))
  noncoding_exons_multiplex$feature <- feature
  noncoding_exons <- noncoding_exons_multiplex %>%
    bind_rows(noncoding_exons_singleton)

  ## reconstruct the transcriptome
  tx <- bind_rows(utr5, cds, utr3, noncoding_exons) %>%
    dplyr::left_join(dplyr::select(tx, tx_name, gene),
      by = c("transcript" = "tx_name")
    ) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::mutate(exon = paste(seqnames, feature,
      paste(transcript, exon_rank, sep = "_"),
      sep = ":"
    ))
  
  ## clean up
  rm(utr5, cds, utr3, noncoding_exons, singleton,
    noncoding_exons_multiplex, noncoding_exons_singleton,
    multiplex, nexons, exons)
  gc()
  
  ## Gene ID type
  keytype <- ""
  if (any(grepl("^ENS", tx$gene[1:10], perl = TRUE))) {
    keytype <- "GENEID"
    ensembl_gene_id <- gsub("\\.\\d+$", "", tx$gene, perl = TRUE)
  } else if (any(grepl("^\\d+$", tx$gene[1:10], perl = TRUE))) {
    keytype <- "ENTREZID"
    ensembl_gene_id <- tx$gene
  }

  if (keytype == "GENEID" || keytype == "ENTREZID") {
    symbol <- AnnotationDbi::mapIds(edb,
      keys = unique(ensembl_gene_id),
      column = "GENENAME",
      keytype = keytype,
      multiVals = "list"
    )
    symbol <- sapply(symbol, function(.ele) {
      paste(unique(.ele[!is.na(.ele)]), collapse = ";")
    })
    symbol <- ifelse(symbol[ensembl_gene_id] == "",
      ensembl_gene_id, symbol[ensembl_gene_id]
    )
    tx <- tx %>% dplyr::mutate(symbol = symbol[ensembl_gene_id])
  } else {
    tx <- tx %>% dplyr::mutate(symbol = gene)
  }
  tx <- tx %>% plyranges::as_granges()
}
