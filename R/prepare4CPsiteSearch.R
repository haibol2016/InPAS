#' prepare data for predicting cleavage and polyadenylation (CP) sites
#'
#' @param coverage coverage for each sample, output from
#'   [coverageFromBedGraph()]
#' @param groupList a list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param utr3 output of [utr3Annotation()]
#' @param background the range for calculating cutoff threshold of local
#'   background
#' @param TxDb an object of [GenomicFeatures::TxDb-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be
#' removed from the genome. If you use a TxDb containing alternative
#' scaffolds, you'd better to remove the scaffolds.
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

prepare4CPsiteSearch <- function(coverage, 
                        groupList = NULL,
                        genome, 
                        utr3,
                        background = c(
                            "same_as_long_coverage_threshold",
                            "1K", "5K", "10K", "50K"
                        ),
                        TxDb = NA,
                        removeScaffolds = FALSE,
                        BPPARAM = NULL,
                        tmpfolder = NULL,
                        silence = TRUE){
    
    gcCompensation <- NA
    mappabilityCompensation <- NA
    FFT <- FALSE
    fft.sm.power <- 20

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
    if (!is.null(tmpfolder))
    {
        dir.create(tmpfolder, recursive = TRUE, showWarnings = FALSE)
    }
    
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
            
            # merge the coverage data of each sample in a given group together
            coverage <-
                mergeCoverage(coverage, groupList,
                              genome,
                              removeScaffolds = removeScaffolds,
                              BPPARAM = BPPARAM
                )
        }
        
        if (!is.null(tmpfolder) && !cov.load) {
            save(list = "coverage", file = file.path(tmpfolder, "cov.rds"))
        }
        
        # only the number of groups of merged coverage
        hugeData <- FALSE
        depth.weight <- depthWeight(coverage, hugeData = hugeData)
        totalCov <- totalCoverage(coverage, genome, 
                                  hugeData = hugeData,
                                  removeScaffolds = removeScaffolds)
    } else {
        depth.weight <- depthWeight(coverage, hugeData = hugeData, groupList)
        totalCov <- totalCoverage(coverage = coverage, 
                                  genome = genome, 
                                  hugeData = hugeData, 
                                  removeScaffolds = removeScaffolds, 
                                  groupList)
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
    list(utr3TotalCov = utr3TotalCov, background = background,
         z2s = z2s, depth.weight = depth.weight)
}
