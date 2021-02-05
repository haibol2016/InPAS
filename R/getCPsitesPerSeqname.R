#' Estimate the CP sites
#'
#' Estimate the CP sites for UTRs on a given chromosome
#'
#' @param chr.cov Coverage list for one chromosome
#' @param utr3 output of [utr3Annotation()]
#' @param background how to get the local background
#' @param z2s output of [zScoreThreshold()]
#' @param depth.weight output of [depthWeight()]
#' @param genome  A [BSgenome::BSgenome-class] object
#' @param MINSIZE min size of short form
#' @param window_size window size
#' @param search_point_START search start point
#' @param search_point_END search end point
#' @param cutStart  cut from start
#' @param cutEnd cut from end
#' @param adjust_distal_polyA_end adjust distal site or not
#' @param coverage_threshold cutoff value for coverage
#' @param long_coverage_threshold cutoff value for long form
#' @param PolyA_PWM  polyA PWM
#' @param classifier  classifier
#' @param classifier_cutoff classifier cutoff
#' @param shift_range shift range
#' @param step Adjust step, default 1, means adjust by each base by
#'   cleanUpdTSeq.
#' @param two_way Search the proximal site from both direction or not
#' @param hugeData logical(1). Is this dataset consume too much memory? if it is
#'   TRUE, the coverage will be saved into temporary files.
#' @param tmpfolder temp folder could save and reload the analysis data for
#'   resume analysis
#' @param silence logical(1), indicating whether progress is reported or not. By
#'   default, FALSE
#'
#' @return a data frame
#' @seealso [CPsites()],[searchProximalCPs()], [proximalAdj()],
#'   [proximalAdjByPWM()], [proximalAdjByCleanUpdTSeq()],[PAscore()],
#'   [PAscore2()]
#' @import GenomicRanges
#' @importFrom BSgenome getSeq matchPWM
#' @keywords internal

estimateCPsites <- function(chr.cov, 
                            utr3,
                            background, 
                            z2s,
                            depth.weight,
                            genome, 
                            MINSIZE = 10, 
                            window_size,
                            search_point_START,
                            search_point_END,
                            cutStart, cutEnd,
                            adjust_distal_polyA_end,
                            coverage_threshold,
                            long_coverage_threshold,
                            PolyA_PWM, 
                            classifier,
                            classifier_cutoff,
                            shift_range,
                            step = 1,
                            two_way = FALSE,
                            hugeData = TRUE,
                            tmpfolder = NULL, 
                            silence = FALSE) {
  if (!is.na(PolyA_PWM)[1]) {
    if (!is(PolyA_PWM, "matrix")) stop("PolyA_PWM must be matrix")
    if (any(rownames(PolyA_PWM) != c("A", "C", "G", "T"))) {
      stop("rownames of PolyA_PWM must be c('A', 'C', 'G', 'T')")
    }
  }
  if (missing(chr.cov) || missing(genome) || missing(utr3)) {
    stop("coverage for a chromosome/sc, genome and utr3 are required.")
  }
  if (!is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome.")
  }
  if (!is(utr3, "GRanges") |
      !all(utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
    stop("utr3 must be output of function of utr3Annotation")
  }
  
  if (seqlevelsStyle(utr3) != seqlevelsStyle(genome)) {
    stop("the seqlevelsStyle of utr3 must be same as genome")
  }
  if (hugeData && is.null(tmpfolder)){
    stop("An explicit tempfolder is required for huge data")
  }
  
  if (is.character(chr.cov)) {
    chr.cov <- readRDS(file = chr.cov)
    if (!is.list(chr.cov)) {
      stop(
        "Something wrong when load big data. ",
        "Maybe the tempfile is broken!"
      )
    }
  }
  chr.cov <- chr.cov[sapply(chr.cov, mean) > 0]
  if (length(chr.cov) == 0) {
    return(NULL)
  }
  curr_UTR.gr <- utr3[names(chr.cov)]
  seqn <- as.character(seqnames(curr_UTR.gr))[1]

  utr3.utr <- curr_UTR.gr[curr_UTR.gr$feature == "utr3"]
  utr3.gap <- curr_UTR.gr[curr_UTR.gr$feature == "next.exon.gap"]
  co <- countOverlaps(utr3.gap, utr3.utr,
    maxgap = 1, ignore.strand = TRUE
  )
  utr3.gap <- utr3.gap[co > 1]
  curr_UTR.gr$conn_next_utr3 <-
    curr_UTR.gr$transcript %in% utr3.gap$transcript
  curr_UTR <- split(curr_UTR.gr, curr_UTR.gr$transcript)
  conn_next_utr3 <- sapply(curr_UTR, function(.UTR) {
    .UTR$conn_next_utr3[1]
  })
  chr.cov.merge <- lapply(curr_UTR, function(.UTR) {
    .UTR <- .UTR[order(start(.UTR))]
    chr.utr3TotalCov <- chr.cov[names(.UTR)]
    chr.utr3TotalCov <-
      mapply(function(.covList, .start, .end, .property) {
        # set names for each position
        .posList <- .start:.end
        if (length(dim(.covList)) == 0) .covList <- t(.covList)
        rownames(.covList) <-
          paste(.property, .posList, sep = "_SEP_")
        .covList
      }, chr.utr3TotalCov,
      start(.UTR), end(.UTR),
      .UTR$feature,
      SIMPLIFY = FALSE
      )

    chr.utr3TotalCov <- do.call(rbind, chr.utr3TotalCov)
    if (as.character(strand(.UTR))[1] == "-") { ## reverse the negative strand
      chr.utr3TotalCov <-
        chr.utr3TotalCov[rev(rownames(chr.utr3TotalCov)), ,
          drop = FALSE
        ]
    }

    # if the range of "cutstart" is (0, 1), percentage; otherwise absolute bases
    if (!is.na(cutStart)) {
      if (cutStart < 1) {
        cutStart <- floor(length(chr.utr3TotalCov) * cutStart)
      }
      if (cutStart > 0) {
        chr.utr3TotalCov <-
          chr.utr3TotalCov[-(1:cutStart), , drop = FALSE]
      }
    }
    chr.utr3TotalCov
  })
  if (!silence) message("chromsome", seqn, "coverage merged.")
  ## chr.cov.merge should be a list with named numeric
  ## filter the coverage
  coverage_quality <- sapply(chr.cov.merge, function(.ele) {
    if (nrow(.ele[grepl("utr3_SEP_", rownames(.ele)), ,
      drop = FALSE]) > MINSIZE) {
      any(colMeans(.ele[1:min(nrow(.ele), 100), ,
        drop = FALSE]) > coverage_threshold)
    } else {
      FALSE
    }
  })

  chr.cov.merge <- chr.cov.merge[coverage_quality]
  conn_next_utr3 <- conn_next_utr3[coverage_quality]
  if (!silence) {
    message(
      "chromsome", seqn,
      "quality filtered by first 100nt."
    )
  }

  if (!is.null(tmpfolder)) {
    if (file.exists(file.path(tmpfolder, seqn))) {
      load(file.path(tmpfolder, seqn))
      if (!silence) {
        message(
          "chromsome",
          seqn, "coverage loaded."
        )
      }
    }
  }
  if (!exists("CPsite_search_step")) {
    CPsite_search_step <- 0
  }

  if (length(chr.cov.merge) > 0) {
    ## step1 search distal cp sites
    if (!silence) message("chromsome", seqn, "distal search ... start")
    curr_UTR <- curr_UTR[names(chr.cov.merge)]
    if (CPsite_search_step < 1) {
      CPsite_search_step <- 1
      chr.abun <- searchDistalCPs(
        chr.cov.merge,
        conn_next_utr3,
        curr_UTR,
        window_size,
        depth.weight,
        long_coverage_threshold,
        background,
        z2s
      )
      if (!is.null(tmpfolder)) {
        save(
          list = c("chr.abun", "CPsite_search_step"),
          file = file.path(tmpfolder, seqn)
        )
      }
      if (!silence) message("chromsome", seqn, "distal search ... done")
    }

    ## step2 adjust distal CP sites
    if (!silence) message("chromsome", seqn, "distal adjust ... start")
    if (adjust_distal_polyA_end && is(classifier, "PASclassifier")) {
      if (CPsite_search_step < 2) {
        CPsite_search_step <- 2
        chr.abun <- distalAdj(
          chr.abun, classifier, classifier_cutoff,
          shift_range, genome, step
        )
        if (!is.null(tmpfolder)) {
          save(
            list = c("chr.abun", "CPsite_search_step"),
            file = file.path(tmpfolder, seqn)
          )
        }
        if (!silence) {
          message("chromsome", seqn, "distal adjust ... done")
        }
      }
    }

    ## step3 search proximal CP sites
    if (!silence) message("chromsome", seqn, "proximal search ... start")
    if (CPsite_search_step < 3) {
      CPsite_search_step <- 3
      chr.abun <- searchProximalCPs(
        chr.abun, curr_UTR,
        window_size, MINSIZE,
        cutEnd,
        search_point_START,
        search_point_END,
        two_way
      )
      if (!is.null(tmpfolder)) {
        save(
          list = c("chr.abun", "CPsite_search_step"),
          file = file.path(tmpfolder, seqn)
        )
      }
      if (!silence) {
        message("chromsome", seqn, "proximal searched ... done")
      }
    }

    ## step4 adjust proximal CP sites
    if (!silence) message("chromsome", seqn, "proximal adjust ... start")
    if (is(PolyA_PWM, "matrix") || is(classifier, "PASclassifier")) {
      if (CPsite_search_step < 4) {
        CPsite_search_step <- 4
        chr.abun <- proximalAdj(
          chr.abun,
          MINSIZE,
          PolyA_PWM,
          genome,
          classifier,
          classifier_cutoff,
          shift_range,
          search_point_START,
          step
        )
        if (!is.null(tmpfolder)) {
          save(
            list = c("chr.abun", "CPsite_search_step"),
            file = file.path(tmpfolder, seqn)
          )
        }
        if (!silence) {
          message("chromsome", seqn, "proximal adjust ... done")
        }
      }
    }
    chr.abun <- polishCPs(chr.abun)
  } else {
    chr.abun <- NULL
  }
  
  if (!is.null(chr.abun)){
    utr3.shorten.UTR <- utr3 %>%
      plyranges::filter(feature == "utr3") %>%
      plyranges::select(-c(feature)) %>%
      plyranges::filter(transcript %in%
                          rownames(chr.abun))
    
    chr.abun <-
      chr.abun[utr3.shorten.UTR$transcript, , drop = FALSE]
    
    utr3.shorten.UTR$fit_value <- unlist(chr.abun[, "fit_value"])
    utr3.shorten.UTR$Predicted_Proximal_APA <-
      unlist(chr.abun[, "Predicted_Proximal_APA"])
    utr3.shorten.UTR$Predicted_Distal_APA <-
      unlist(chr.abun[, "Predicted_Distal_APA"])
    utr3.shorten.UTR$type <- unlist(chr.abun[, "type"])
    utr3.shorten.UTR$utr3start <- unlist(chr.abun[, "utr3start"])
    utr3.shorten.UTR$utr3end <- unlist(chr.abun[, "utr3end"])
    utr3.shorten.UTR <-
      utr3.shorten.UTR[!is.na(utr3.shorten.UTR$Predicted_Proximal_APA)]
  } else {
    utr3.shorten.UTR <- NULL
  }
  utr3.shorten.UTR
}
