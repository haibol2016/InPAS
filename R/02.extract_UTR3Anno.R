#' extract 3' UTR information from a [GenomicFeatures::TxDb-class] object
#'
#' extract 3' UTR information from a [GenomicFeatures::TxDb-class] object. The
#' 3'UTR is defined as the last 3'UTR fragment for each transcript and it will
#' be cut if there is any overlaps with other exons.
#'
#' @details A good practice is to perform read alignment using a reference
#'   genome from Ensembl/GenCode including only the primary assembly and build a
#'   TxDb using the GTF/GFF files downloaded from the same source as the
#'   reference genome, such as BioMart/Ensembl/GenCode. For instruction, see
#'   Vignette of the GenomicFeatures. The UCSC reference genomes and their
#'   annotation can be very cubersome.
#'
#' @param TxDb an object of [GenomicFeatures::TxDb-class]
#' @param edb An object of [ensembldb::EnsDb-class]
#' @param removeScaffolds A logical(1) vector, whether the scaffolds should be
#'   removed from the genome If you use a TxDb containing alternative
#'   scaffolds, you'd better to remove the scaffolds.
#' @param MAX_EXONS_GAP An integer(1) vector, maximal gap sizes between last
#'   known CP sites to downstream exons
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom plyranges as_granges complement_ranges disjoin_ranges filter group_by mutate reduce_ranges reduce_ranges_directed remove_names select set_genome_info shift_downstream summarise
#' @importFrom dplyr as_tibble mutate filter arrange bind_rows group_by
#'   left_join summarise n
#'
#' @author Jianhong Ou, Haibo Liu
#' @return An object of [GenomicRanges::GRangesList], containing GRanges for
#'   extracted 3' UTRs, and the corresponding last CDSs and next.exon.gap for 
#'   each chromosome/scaffold.
#'
#' @export
#'
#' @examples
#' library("EnsDb.Hsapiens.v86")
#' library("GenomicFeatures")
#' samplefile <- system.file("extdata",
#'                           "hg19_knownGene_sample.sqlite",
#'                            package = "GenomicFeatures")
#' TxDb <- loadDb(samplefile)
#' edb <- EnsDb.Hsapiens.v86
#' utr3 <- extract_UTR3Anno(TxDb, edb,
#'                  removeScaffolds = TRUE,
#'                  MAX_EXONS_GAP = 10000)

extract_UTR3Anno <- function(TxDb = NULL, edb = NULL,
                           removeScaffolds = FALSE,
                           MAX_EXONS_GAP = 10000) {
  if (is.null(TxDb) || !is(TxDb, "TxDb")) {
    stop("TxDb is required and must be an object of GenomicFeatures::TxDb-class!")
  }
  if (is.null(edb) || !is(edb, "EnsDb")) {
    stop("edb is required and must be an object of ensembldb::EnsDb-class!")
  }

  tx <- parse_TxDb(TxDb, edb, removeScaffolds = removeScaffolds)
  ## extract utr3 sharing same start positions from same gene
  if (any(is.na(tx$gene) | tx$gene == "")) {
    stop("unexpected things happend at checking gene IDs", 
         " gene should not contain NA or empty strings")
  }

  
  utr3 <- tx %>%
    plyranges::filter(feature != "ncRNA" &
      feature %in% c("utr3", "lastutr3")) %>%
    unique()

  utr3.last <- utr3 %>% plyranges::filter(feature == "lastutr3")

  exons.not.utr3.last <- tx %>%
    plyranges::filter(feature != "lastutr3") %>%
    unique()

  utr3.last.block <- utr3.last

  ## reduce GRanges to cover from min_start to max_end
  getRange <- function(gr, f) {
    gr <- gr %>%
      plyranges::group_by(seqnames, !!!syms(f), strand) %>%
      plyranges::summarise(start = min(start), end = max(end)) %>%
      plyranges::as_granges()
    gr
  }

  ## mask 5UTR-CDS and any other region from last 3UTR
  # non.utr3 <- getRange(exons.not.utr3.last, "transcript") %>%
  #    plyranges::filter(width < 1000) ##why 10000? because some gene is unbelivable huge.
  # non.utr3 <- reduce(c(non.utr3, reduce(exons.not.utr3.last)))
  non.utr3 <- reduce(exons.not.utr3.last)

  removeMask <- function(gr, mask) {
    # mask is already reduced
    # mask.reduce <- reduce(mask)
    mask <- non.utr3
    gr <- utr3.last.block
    ol.mask <- findOverlaps(mask, gr,
      ignore.strand = TRUE, minoverlap = 1
    )
    if (length(ol.mask) > 0) {
      gr.ol <- gr[subjectHits(ol.mask)]
      gr.nol <- gr[-unique(subjectHits(ol.mask))]
      mask.ol <- mask[queryHits(ol.mask)]
      gr.ol$rel <- "ol"
      gr.ol$rel[end(mask.ol) >= start(gr.ol) &
        start(mask.ol) <= start(gr.ol)] <- "left"
      gr.ol$rel[end(mask.ol) >= end(gr.ol) &
        start(mask.ol) <= end(gr.ol)] <- "right"
      # remove covered gr by mask
      gr.ol$rel[end(mask.ol) >= end(gr.ol) &
        start(mask.ol) <= start(gr.ol)] <- "cover"
      # Haibo add "-1" and " +1"
      end(gr.ol)[gr.ol$rel == "right"] <-
        start(mask.ol)[gr.ol$rel == "right"] - 1
      start(gr.ol)[gr.ol$rel == "left"] <-
        end(mask.ol)[gr.ol$rel == "left"] + 1

      # mask in bewteen gr
      start(gr.ol)[gr.ol$rel == "ol" &
        as.character(strand(gr.ol)) == "+"] <-
        end(mask.ol)[gr.ol$rel == "ol" &
          as.character(strand(gr.ol)) == "+"] + 1
      end(gr.ol)[gr.ol$rel == "ol" &
        as.character(strand(gr.ol)) == "-"] <-
        start(mask.ol)[gr.ol$rel == "ol" &
          as.character(strand(gr.ol)) == "-"] - 1
      gr.ol <- gr.ol[gr.ol$rel != "cover"]
      gr.ol$rel <- NULL

      gr.ol <- unique(gr.ol)
      exons.dup <- gr.ol[gr.ol$exon %in%
        gr.ol$exon[duplicated(gr.ol$exon)]]
      if (length(exons.dup) > 0) {
        exons.nd <- gr.ol[!gr.ol$exon %in%
          gr.ol$exon[duplicated(gr.ol$exon)]]
        exons.dup <- split(exons.dup, exons.dup$exon)
        exons.dup <- lapply(exons.dup, function(.ele) {
          # .ele <- exons.dup[["chr1|ENST00000667728.1_2|lastutr3"]]
          ## get disjoined parts
          .disj <- disjoin(.ele)
          .cnt <- countOverlaps(.disj, .ele)
          .disj <- .disj[order(-.cnt)]
          ## get the last 3' UTR segment not masked by any non-utr frags
          .disj <- .disj[1]
          mcols(.disj) <- mcols(.ele[1])
          .disj
        })
        exons.dup <- unlist(GRangesList(exons.dup))
        names(exons.dup) <- NULL
        if (length(exons.nd) > 0) {
          gr.ol <- c(exons.dup, exons.nd)
        } else {
          gr.ol <- exons.dup
        }
      }
      gr <- c(gr.ol, gr.nol)
    } else {
      gr
    }
  }
  ####################################################################
  #               MASK ONE by non-UTR exonic regions                 #
  ####################################################################

  utr3.last <- removeMask(utr3.last, non.utr3)
  utr3.last <- unique(utr3.last)
  ####################################################################
  #           MASK TWO by 3' UTR on the opposite strand              #
  ####################################################################
  utr3.reverse.strand <- utr3.last %>%
    plyranges::mutate(strand = ifelse(strand == "+", "-", "+")) %>%
    plyranges::reduce_ranges_directed()

  ol.utr3.rev <- findOverlaps(utr3.last, utr3.reverse.strand,
    ignore.strand = FALSE, minoverlap = 1
  )
  if (length(ol.utr3.rev) > 0) {
    utr3.last.ol <- utr3.last[queryHits(ol.utr3.rev)]
    utr3.last.nol <- utr3.last[-unique(queryHits(ol.utr3.rev))]
    utr3.rev.ol <- utr3.reverse.strand[subjectHits(ol.utr3.rev)]
    strand <- as.character(strand(utr3.last.ol)) == "+"
    idx <- (start(utr3.last.ol) < start(utr3.rev.ol) & strand) |
      (end(utr3.last.ol) > end(utr3.rev.ol) & !strand)
    utr3.last.ol <- utr3.last.ol[idx]
    utr3.rev.ol <- utr3.rev.ol[idx]
    strand <- strand[idx]
    end(utr3.last.ol)[strand] <- start(utr3.rev.ol)[strand] - 1
    start(utr3.last.ol)[!strand] <- end(utr3.rev.ol)[!strand] + 1
    utr3.last <- c(utr3.last.nol, utr3.last.ol[width(utr3.last.ol) > 1])
  }

  ## keep the last 3' utr for each transcripts
  utr3.last$feature <- "utr3"
  utr3.last <- assign_feature(utr3.last, feature_alt = "lastutr3") %>%
    plyranges::filter(feature == "lastutr3", width > 1) %>%
    unique()

  #####################################################################
  #      Collapsing last UTRs of the same start/end coordinates       #
  #                on the +/- strand, respective per gene             #
  #####################################################################

  collapse_same_start_utr <- function(gr) {
    utr3.last <- gr
    utr3.last <- utr3.last %>%
      unique() %>%
      plyranges::mutate(
        start.utr3.last = paste(seqnames, start,
          gene,
          sep = ":"
        ),
        end.utr3.last = paste(seqnames, end,
          gene,
          sep = ":"
        )
      )

    start.dup <- with(
      utr3.last,
      unique(start.utr3.last[duplicated(start.utr3.last) &
        as.character(strand(utr3.last)) == "+"])
    )
    end.dup <- with(
      utr3.last,
      unique(end.utr3.last[duplicated(end.utr3.last) &
        as.character(strand(utr3.last)) == "-"])
    )
    utr3.last.dup <- utr3.last %>%
      plyranges::filter(start.utr3.last %in% start.dup |
        end.utr3.last %in% end.dup)

    if (length(utr3.last.dup) > 0) {
      utr3.last.nd <- utr3.last %>%
        plyranges::filter(!(start.utr3.last %in% start.dup |
          end.utr3.last %in% end.dup)) %>%
        plyranges::select(-c(start.utr3.last, end.utr3.last)) %>%
        dplyr::mutate(annotatedProximalCP = "unknown")

      utr3.last.dup <- utr3.last.dup %>%
        plyranges::mutate(
          dup.group =
            ifelse(strand == "+",
              start.utr3.last, end.utr3.last
            )
        ) %>%
        plyranges::select(-c(start.utr3.last, end.utr3.last)) %>%
        data.frame() %>%
        dplyr::as_tibble() %>%
        dplyr::arrange(dup.group, desc(width))

      utr3_end <- ifelse(utr3.last.dup$strand == "+",
        utr3.last.dup$end, utr3.last.dup$start
      )
      utr3_end <- split(utr3_end, utr3.last.dup$dup.group)
      annotation <- lapply(utr3_end, function(.x) {
        label <- paste("proximalCP",
          paste(.x[-1], collapse = "_"),
          sep = "_"
        )
      })
      annotation <- do.call("c", annotation)
      utr3.last.dup <- utr3.last.dup %>%
        dplyr::filter(!duplicated(dup.group)) %>%
        dplyr::mutate(
          annotatedProximalCP =
            annotation[dup.group]
        ) %>%
        dplyr::select(-dup.group) %>%
        as_granges()

      utr3.last <- c(utr3.last.nd, utr3.last.dup)
    } else {
      utr3.last <- gr %>%
        plyranges::mutate(annotatedProximalCP = "unknown")
    }
    utr3.last
  }
  utr3.last <- collapse_same_start_utr(utr3.last)

  ######################################################################
  #               Mark truncated 3UTR (3821 truncated utr3)            #
  ######################################################################

  utr3.last$truncated <- FALSE
  utr3.last$truncated[as.character(strand(utr3.last)) == "+" &
    end(utr3.last) <
      end(utr3.last.block[match(
        utr3.last$exon,
        utr3.last.block$exon
      )])] <- TRUE
  utr3.last$truncated[as.character(strand(utr3.last)) == "-" &
    start(utr3.last) >
      start(utr3.last.block[match(
        utr3.last$exon,
        utr3.last.block$exon
      )])] <- TRUE

  #####################################################################
  #            split the overlapping UTRs of the same genes           #
  #####################################################################
  utr3.last_reduced <- reduce(utr3.last, ignore.strand = FALSE)
  count_overlap <- countOverlaps(utr3.last_reduced, utr3.last)
  utr3.last_reduced <- utr3.last_reduced[count_overlap > 1]

  utr3.last_ol <- GenomicRanges::findOverlaps(utr3.last_reduced,
    utr3.last,
    ignore.strand = TRUE, minoverlap = 1
  )
  # non-overlapping utr3
  utr3.last_isolated <- utr3.last[-unique(subjectHits(utr3.last_ol))]

  # split overlapping utr3.last of the same genes into
  # compatible segaments for coverage calculation
  # !!! overlapping utr3.last of different genes are removed
  uniq_pairs <- data.frame(utr3.last_ol)

  uniq_pairs <- split(uniq_pairs[, 2], uniq_pairs[1])
  utr3.last_overlap_split <- lapply(uniq_pairs, function(.x) {
    gr <- utr3.last[.x]
    if (length(unique(gr$gene)) > 1) {
      gr <- NULL
    } else {
      gr_reduced <- reduce(gr)

      strand <- strand(gr_reduced)
      feature <- gr$feature[1]
      gene <- gr$gene[1]
      symbol <- gr$symbol[1]

      if (as.character(strand(gr_reduced)) == "+") {
        gr <- gr[order(start(gr))]
        st <- start(gr)
        segment_start <- st
        segment_end <- c(st[-1] - 1, end(gr_reduced))
        end_sort <- gr[order(-end(gr))]
        transcript <- gr$transcript
        exon <- gr$exon
        exon_rank <- gr$exon_rank
        truncated <- c(
          rep(TRUE, length(st) - 1),
          end_sort$truncated[1]
        )
        CP <- grep("^proximalCP_", gr$annotatedProximalCP,
          perl = TRUE, value = TRUE
        )
        if (length(CP) >= 1) {
          CP <- gsub(
            "_proximalCP", "",
            paste(CP, collapse = "_")
          )
          CP <- paste(CP, paste(end(end_sort)[-1],
            collapse = "_"
          ),
          sep = "_"
          )
        } else {
          CP <- paste("proximalCP", paste(end(end_sort)[-1],
            collapse = "_"
          ),
          sep = "_"
          )
        }
      } else {
        gr <- gr[order(end(gr))]
        end <- end(gr)
        segment_end <- end
        segment_start <- c(start(gr_reduced), end[-length(end)] + 1)
        transcript <- gr$transcript
        exon <- gr$exon
        exon_rank <- gr$exon_rank
        start_sort <- gr[order(start(gr))]
        truncated <- c(
          start_sort$truncated[1],
          rep(TRUE, length(end) - 1)
        )
        CP <- grep("^proximalCP_", gr$annotatedProximalCP,
          perl = TRUE, value = TRUE
        )
        if (length(CP) >= 1) {
          CP <- gsub(
            "_proximalCP", "",
            paste(CP, collapse = "_")
          )
          CP <- paste(CP, paste(start(start_sort)[-1],
            collapse = "_"
          ),
          sep = "_"
          )
        } else {
          CP <- paste("proximalCP",
            paste(start(start_sort)[-1],
              collapse = "_"
            ),
            sep = "_"
          )
        }
      }
      gr <- GRanges(data.frame(
        seqnames =
          seqnames(gr_reduced),
        start = segment_start,
        end = segment_end,
        strand = strand(gr_reduced)
      ),
      exon_rank = exon_rank,
      transcript = transcript,
      feature = feature,
      gene = gene,
      exon = exon,
      symbol = symbol,
      annotatedProximalCP = CP,
      truncated = truncated
      )
    }
    gr
  })
  rmv_idx <- do.call("c", lapply(utr3.last_overlap_split, is.null))
  utr3.last_overlap_split <- 
    unlist(GRangesList(utr3.last_overlap_split[!rmv_idx]))
  utr3.last <- c(utr3.last_isolated, utr3.last_overlap_split)

  ###################################################################
  #    reprocess 3' UTRs with the same start                        #
  ###################################################################
  utr3.last.sameStart <- utr3.last[grepl(
    "proximalCP",
    utr3.last$annotatedProximalCP
  )]
  utr3.last.unknown <- utr3.last[utr3.last$annotatedProximalCP == "unknown"]
  utr3.last.sameStart.CP <-
    strsplit(gsub("proximalCP_", "",
                  utr3.last.sameStart$annotatedProximalCP), "_")
  utr3.last.sameStart.Start <- start(utr3.last.sameStart)
  utr3.last.sameStart.End <- end(utr3.last.sameStart)
  utr3.last.sameStart$annotatedProximalCP <- mapply(function(cp, st, en) {
    cp <- as.integer(cp)
    cp <- cp[cp > st & cp < en]
    if (length(cp) == 0) {
      return("unknown")
    } else {
      return(paste("proximalCP", paste(cp, collapse = "_"), sep = "_"))
    }
  }, utr3.last.sameStart.CP,
  utr3.last.sameStart.Start,
  utr3.last.sameStart.End,
  SIMPLIFY = TRUE
  )

  utr3.clean <- c(utr3.last.unknown, utr3.last.sameStart)
  utr3.clean <-
    utr3.clean[order(
      as.character(seqnames(utr3.clean)),
      start(utr3.clean)
    )]

  ####################################################################
  #       Determine the gap between the untracated 3' UTRs end       #
  #                 and the next downstream  exon                    #
  ####################################################################
  genome.info <- as.data.frame(seqinfo(TxDb))[, 1:2]
  gaps <- tx %>%
    plyranges::reduce_ranges() %>%
    plyranges::set_genome_info(
      seqnames = rownames(genome.info),
      seqlengths = genome.info$seqlengths,
      is_circular = genome.info$isCircular
    ) %>%
    plyranges::complement_ranges()

  utr3.clean.ext1 <- utr3.clean %>%
    plyranges::filter(!truncated) %>%
    plyranges::shift_downstream(shift = 1L)

  ol <- findOverlaps(utr3.clean.ext1, gaps, ignore.strand = TRUE)
  ol.utr3.clean <- utr3.clean[queryHits(ol)]

  next.exons.gap <- gaps[subjectHits(ol)] %>%
    plyranges::mutate(strand = strand(ol.utr3.clean))
  mcols(next.exons.gap) <- mcols(ol.utr3.clean)

  wid <- width(next.exons.gap) > MAX_EXONS_GAP
  width(next.exons.gap)[wid & as.character(strand(next.exons.gap)) == "+"] <-
    MAX_EXONS_GAP
  start(next.exons.gap)[wid & as.character(strand(next.exons.gap)) == "-"] <-
    end(next.exons.gap)[
      wid & as.character(strand(next.exons.gap)) == "-"
    ] - MAX_EXONS_GAP + 1
  next.exons.gap <- next.exons.gap %>%
    plyranges::mutate(feature = "next.exon.gap")
  utr3.clean <- utr3.clean %>%
    plyranges::mutate(feature = "utr3")

  ###################################################################
  #             Get last CDS for coverage compensation              #
  ###################################################################

  CDS.last <- tx %>%
    plyranges::filter(feature == "lastCDS") %>%
    plyranges::filter(transcript %in% utr3.clean$transcript)
  CDS.last <- CDS.last[match(
    utr3.clean$transcript,
    CDS.last$transcript
  )]
  mcols(CDS.last) <- mcols(utr3.clean)
  CDS.last$feature <- "CDS"

  utr3.fixed <- c(utr3.clean, next.exons.gap, CDS.last)

  names(utr3.fixed) <-
    paste(utr3.fixed$exon, utr3.fixed$symbol,
      utr3.fixed$feature,
      sep = "|"
    )
  utr3.fixed <- split(utr3.fixed, seqnames(utr3.fixed))
  utr3.fixed
}
