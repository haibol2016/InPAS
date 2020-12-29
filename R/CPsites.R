#' predict cleavage and polyadenylation (CP) sites
#'
#' predict alternative cleavage and polyadenylation (CP or APA) sites.
#'
#' @param coverage coverage for each sample, output from
#'   [coverageFromBedGraph()]
#' @param groupList a list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param utr3 output of [utr3Annotation()]
#' @param window_size window size for novel distal CP site searching and
#'   adjusted CP site searching, default: 100
#' @param search_point_START an positive integer. start point for searching 
#'   proximal CP sites relative to the 5' extremity of a given extracted 3' UTR 
#'   entity.
#' @param search_point_END an positive integer. end point for searching 
#'   searching proximal CP sites relative to the 5' extremity of a given
#'   extracted 3' UTR entity.
#' @param cutStart what percentage or how many nucleotides should be removed
#'   from the start before searching for CP sites? It can be a decimal in 
#'   [0, 1) or non-negative integer. 0.1 means 10 percent, 25 means cut first
#'   25 bases
#' @param cutEnd what percentage or how many nucleotides should be removed from
#'   the end before searching for CP sites. See cutStart
#' @param adjust_distal_polyA_end If true, adjust distal polyA end by
#'   [cleanUpdTSeq::cleanUpdTSeq-package]
#' @param coverage_threshold cutoff  threshold of coverage for first 100
#'   nucleotides. If the coverage of first 100 nucleotides is lower than
#'   coverage_threshold, that transcript will be dropped
#' @param long_coverage_threshold cutoff threshold for coverage in the region of
#'   long 3' UTR form. If the coverage in the region of long form is less than
#'   long_coverage_threshold, that transcript will be dropped
#' @param background the range for calculating cutoff threshold of local
#'   background
#' @param TxDb an object of [GenomicFeatures::TxDb-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be
#' removed from the genome. If you use a TxDb containing alternative
#' scaffolds, you'd better to remove the scaffolds.
#' @param PolyA_PWM position Weight Matrix of polyA signal, such as AAUAAA.
#' @param classifier an object of [cleanUpdTSeq::PASclassifier-class]
#' @param classifier_cutoff floating point number between 0 and 1. This is the
#'   cutoff used to assign whether a putative pA is true or false. For example,
#'   classifier_cutoff = 0.5 will assign an putative pA site with prob.1 > 0.5
#'   to the True class (1), and any putative pA site with prob.1 <= 0.5 as False
#'   (0)
#' @param step integer, adjusting step size, default 1, means adjust by each
#'   base using the cleanUpdTSeq algorithm
#' @param two_way logical, Search the proximal site from both direction or not
#' @param shift_range integer, the shift range for polyA site searching
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply.
#' @param tmpfolder temp folder which could be used to save and reload the
#'   analysis data for resuming analysis
#' @param silence report progress or not. By default it doesn't report progress.
#' @references Cheung MS, Down TA, Latorre I, Ahringer J. Systematic bias in
#'   high-throughput sequencing data and its correction by BEADS. Nucleic Acids
#'   Res. 2011 Aug;39(15):e103. doi: 10.1093/nar/gkr425. Mappability could be
#'   calculated by [GEM](http://algorithms.cnag.cat/wiki/Man:gem-mappability)
#'   Derrien T, Estelle J, Marco Sola S, Knowles DG, Raineri E, Guigo R, Ribeca
#'   P. Fast computation and applications of genome mappability. PLoS One.
#'   2012;7(1):e30377. doi: 10.1371/journal.pone.0030377.
#'
#' @return An object of \code{\link[GenomicRanges]{GRanges}}
#' @import S4Vectors Biobase BiocParallel GenomicRanges
#'   GenomicFeatures methods
#' @importFrom BSgenome getSeq matchPWM
#' @importFrom utils object.size
#' @importFrom plyranges as_granges complement_ranges disjoin_ranges filter
#'   group_by mutate reduce_ranges reduce_ranges_directed remove_names select
#'   set_genome_info shift_downstream summarise
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle mapSeqlevels seqlevels
#' @export
#'
#' @examples
#' if (interactive()) {
#'   library(BSgenome.Mmusculus.UCSC.mm10)
#'   bedgraphs <- system.file("extdata",
#'     "Baf3.extract.bedgraph",
#'     package = "InPAS"
#'   )
#'   data(utr3.mm10)
#'   tags <- "Baf3"
#'   genome <- BSgenome.Mmusculus.UCSC.mm10
#'   coverage <- coverageFromBedGraph(bedgraphs,
#'     tags, genome,
#'     removeScaffolds = FALSE,
#'     hugeData = FALSE
#'   )
#'   CP <- CPsites(
#'     coverage = coverage,
#'     groupList = tags,
#'     genome = genome,
#'     utr3 = utr3.mm10,
#'     coverage_threshold = 5,
#'     long_coverage_threshold = 5
#'   )
#' }
CPsites <- function(coverage, groupList = NULL,
                    genome, utr3,
                    window_size = 100,
                    search_point_START = 50,
                    search_point_END = NA,
                    cutStart = window_size,
                    cutEnd = 0,
                    adjust_distal_polyA_end = TRUE,
                    coverage_threshold = 5,
                    long_coverage_threshold = 2,
                    background = c(
                      "same_as_long_coverage_threshold",
                      "1K", "5K", "10K", "50K"
                    ),
                    TxDb = NA,
                    removeScaffolds = FALSE,
                    PolyA_PWM = NA, classifier = NA,
                    classifier_cutoff = .8, step = 1,
                    two_way = FALSE,
                    shift_range = window_size,
                    BPPARAM = NULL,
                    tmpfolder = NULL,
                    silence = TRUE) {
  gcCompensation <- NA
  mappabilityCompensation <- NA
  FFT <- FALSE
  fft.sm.power <- 20
  if (!is.na(PolyA_PWM)[1]) {
    if (!is(PolyA_PWM, "matrix")) stop("PolyA_PWM must be matrix")
    if (any(rownames(PolyA_PWM) != c("A", "C", "G", "T"))) {
      stop("rownames of PolyA_PWM must be c('A', 'C', 'G', 'T')")
    }
  }
  if (missing(coverage) || missing(genome) || missing(utr3)) {
    stop("coverage, genome and utr3 are required.")
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
  if (!is.null(tmpfolder)) dir.create(tmpfolder, showWarnings = FALSE)
  background <- match.arg(background)
  introns <- GRanges()
  if (background != "same_as_long_coverage_threshold") {
    if (!is(TxDb, "TxDb")) {
      stop("TxDb is missing when you want local background")
    }
    if (!identical(unname(genome(TxDb))[1], unname(genome(utr3))[1])) {
      stop("genome of TxDb should be same as genome of utr3")
    }
    if (removeScaffolds) {
      TxDb_seqlevels <- seqlevels(TxDb)
      seqlevels(TxDb) <- TxDb_seqlevels[grepl("^(chr)?(\\d+|[XY])$", TxDb_seqlevels)]
    }
    introns <- unlist(intronsByTranscript(TxDb, use.names = TRUE)) %>%
      plyranges::disjoin_ranges() ## strand-unaware
    exons <- exons(TxDb) %>% plyranges::disjoin_ranges()

    intron_exons <- disjoin(c(exons, introns))
    ol.exons <- findOverlaps(exons, intron_exons)

    ## pure intronic regions
    introns <- intron_exons[-unique(subjectHits(ol.exons))]
  }
  if (!silence) message("total backgroud ... done.")
  utr3 <- utr3[utr3$feature != "CDS"]
  MINSIZE <- 10
  hugeData <- is.character(coverage[[1]])
  if (hugeData && !is.null(names(groupList))[1]) {
    cov.load <- FALSE
    if (!is.null(tmpfolder)) {
      if (file.exists(file.path(tmpfolder, "cov.rds"))) {
        load(file.path(tmpfolder, "cov.rds"))
        cov.load <- TRUE
      }
    }
    if (!cov.load) {
      coverage <-
        mergeCoverage(coverage, groupList,
          genome,
          BPPARAM = BPPARAM
        )
    }

    if (!is.null(tmpfolder) && !cov.load) {
      save(list = "coverage", file = file.path(tmpfolder, "cov.rds"))
    }
    hugeData <- FALSE
    depth.weight <- depthWeight(coverage, hugeData = hugeData)
    totalCov <- totalCoverage(coverage, genome, hugeData = hugeData)
  } else {
    depth.weight <- depthWeight(coverage, hugeData = hugeData, groupList)
    totalCov <- totalCoverage(coverage, genome, hugeData = hugeData, groupList)
  }
  if (!silence) message("total coverage ... done.")
  z2s <- zScoreThreshold(background, introns, totalCov, utr3)
  if (!silence) message("backgroud around 3utr ... done.")
  if (!is.null(tmpfolder) &&
    file.exists(file.path(tmpfolder, "utr3TotalCov.rds"))) {
    utr3TotalCov <- readRDS(file.path(tmpfolder, "utr3TotalCov.rds"))
  } else {
    utr3TotalCov <-
      UTR3TotalCoverage(utr3, totalCov,
        gcCompensation,
        mappabilityCompensation,
        FFT = FFT,
        fft.sm.power = fft.sm.power
      )
    objSize <- as.numeric(object.size(utr3TotalCov)) / (1024^3)
    if (objSize > 4) {
      ## huge data, try to save the data and load later
      utr3TotalCov <- mapply(function(.ele, .name) {
        if (!is.null(tmpfolder)) {
          utr3TotalCov.tmpfile <-
            file.path(tmpfolder, paste0("utr3TotalCov.", .name))
        } else {
          utr3TotalCov.tmpfile <- tempfile()
        }

        ## save actual coverage to tempfile
        saveRDS(.ele, utr3TotalCov.tmpfile, compress = "xz")
        if (!silence) {
          message(
            "save utr3TotalCov", .name,
            "at", utr3TotalCov.tmpfile,
            " ... done."
          )
        }
        utr3TotalCov.tmpfile
      }, utr3TotalCov, names(utr3TotalCov), SIMPLIFY = FALSE)

      if (!is.null(tmpfolder)) {
        saveRDS(utr3TotalCov, file.path(tmpfolder, "utr3TotalCov.rds"))
      }
    }
  }
  if (!silence) message("coverage in 3utr ... done.")

  rm(list = c("coverage", "totalCov"))
  gc()

  if (!is.null(BPPARAM)) {
    shorten_UTR_estimation <-
      bplapply(utr3TotalCov, estimateCPsites,
        BPPARAM = BPPARAM, utr3 = utr3,
        MINSIZE = MINSIZE,
        window_size = window_size,
        search_point_START = search_point_START,
        search_point_END = search_point_END,
        cutStart = cutStart, cutEnd = cutEnd,
        adjust_distal_polyA_end = adjust_distal_polyA_end,
        background = background,
        z2s = z2s,
        coverage_threshold = coverage_threshold,
        long_coverage_threshold = long_coverage_threshold,
        PolyA_PWM = PolyA_PWM,
        classifier = classifier,
        classifier_cutoff = classifier_cutoff,
        shift_range = shift_range,
        depth.weight = depth.weight,
        genome = genome,
        step = step,
        two_way = two_way,
        tmpfolder = tmpfolder,
        silence = silence
      )
  } else {
    shorten_UTR_estimation <-
      lapply(utr3TotalCov, estimateCPsites,
        utr3 = utr3, MINSIZE = MINSIZE,
        window_size = window_size,
        search_point_START = search_point_START,
        search_point_END = search_point_END,
        cutStart = cutStart,
        cutEnd = cutEnd,
        adjust_distal_polyA_end = adjust_distal_polyA_end,
        background = background,
        z2s = z2s,
        coverage_threshold = coverage_threshold,
        long_coverage_threshold = long_coverage_threshold,
        PolyA_PWM = PolyA_PWM,
        classifier = classifier,
        classifier_cutoff = classifier_cutoff,
        shift_range = shift_range,
        depth.weight = depth.weight,
        genome = genome,
        step = step,
        two_way = two_way,
        tmpfolder = tmpfolder,
        silence = silence
      )
  }
  shorten_UTR_estimation <-
    do.call(
      rbind,
      unname(shorten_UTR_estimation[!sapply(
        shorten_UTR_estimation,
        is.null
      )])
    )
  utr3.shorten.UTR <- utr3 %>%
    plyranges::filter(feature == "utr3") %>%
    plyranges::select(-c(feature)) %>%
    plyranges::filter(transcript %in%
      rownames(shorten_UTR_estimation))

  shorten_UTR_estimation <-
    shorten_UTR_estimation[utr3.shorten.UTR$transcript, , drop = FALSE]

  utr3.shorten.UTR$fit_value <- unlist(shorten_UTR_estimation[, "fit_value"])
  utr3.shorten.UTR$Predicted_Proximal_APA <-
    unlist(shorten_UTR_estimation[, "Predicted_Proximal_APA"])
  utr3.shorten.UTR$Predicted_Distal_APA <-
    unlist(shorten_UTR_estimation[, "Predicted_Distal_APA"])
  utr3.shorten.UTR$type <- unlist(shorten_UTR_estimation[, "type"])
  utr3.shorten.UTR$utr3start <- unlist(shorten_UTR_estimation[, "utr3start"])
  utr3.shorten.UTR$utr3end <- unlist(shorten_UTR_estimation[, "utr3end"])
  utr3.shorten.UTR <-
    utr3.shorten.UTR[!is.na(utr3.shorten.UTR$Predicted_Proximal_APA)]
}
