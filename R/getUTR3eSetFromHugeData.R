#' prepare 3' UTR coverage data for usage test
#'
#' generate a UTR3eSet object with PDUI information for statistic tests
#'
#' @param CPsites outputs of [CPsites()]
#' @param coverage coverage for each sample, outputs of [coverageFromBedGraph()],
#'   [getCovFileList()], or [getCov4SmallExperiment()]
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param utr3 output of [utr3Annotation()]
#' @param UTR3CDS.cov an object of [GenomicRanges::GRanges-class], output of
#'   [integrate3UTRUsage()]
#' @param normalize normalization method
#' @param ... parameter can be passed into
#'   [preprocessCore::normalize.quantiles.robust()]
#' @param BPPARAM  An optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply.
#' @param singleSample prepare data for singleSample analysis? default is FALSE
#'
#' @return An object of [UTR3eSet-class] which contains following elements:
#'   usage: an [GenomicRanges::GRanges-class] object with CP sites info. PDUI: a
#'   matrix of PDUI PDUI.log2: log2 transformed PDUI matrix short: a matrix of
#'   usage of short form long: a matrix of usage of long form if singleSample is
#'   TRUE, one more element, signals, will be included.
#' @export
#' @import GenomicRanges
#' @importFrom preprocessCore colSummarizeAvg colSummarizeMedian
#'   normalize.quantiles normalize.quantiles.robust
#'
#' @examples
#' path <- system.file("extdata", package = "InPAS")
#' load(file.path(path, "CPs.MAQC.rda"))
#' load(file.path(path, "coverage.MAQC.rda"))
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' data(utr3.hg19)
#' getUTR3eSet(
#'   CPsites = CPs,
#'   coverage = coverage,
#'   genome = BSgenome.Hsapiens.UCSC.hg19,
#'   utr3 = utr3.hg19
#' )
getUTR3eSetFromHugeData <- function(CPsites, coverage,
                        genome, utr3,
                        UTR3CDS.cov,
                        normalize = c(
                            "none", "quantiles",
                            "quantiles.robust",
                            "mean", "median"
                        ),
                        ...,
                        BPPARAM = NULL, singleSample = FALSE) {
    if (missing(coverage) || missing(CPsites)) {
        stop("CPsites and coverage are required.")
    }
    if (missing(utr3) || missing(genome)) {
        stop("utr3 and genome are required.")
    }
    if (missing(UTR3CDS.cov) || !is(UTR3CDS.cov, "GRanges")){
        stop("UTR3CDS.cov is required and must be an object of GRanges")
    }
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }
    if (!is(utr3, "GRanges") ||
        !all(utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
        stop("utr3 must be output of function of utr3Annotation")
    }
    normalize <- match.arg(normalize)
    hugeData <- is.character(coverage[[1]])
    if ((!singleSample) && length(coverage) == 1) {
        message("Single sample mode is on")
        singleSample <- TRUE
    }
    if (singleSample) {
        if (length(coverage) > 1){
            message("Only first sample will be used.")
        } 
        UTRusage <- UTR3CDS.cov %>% plyranges::filter(!is.na(source))
    } else {
        UTRusage <- UTR3CDS.cov %>% plyranges::filter(!is.na(source))
    }
    
    UTRusage <- split(UTRusage, UTRusage$transcript)
    UTRusage <- UTRusage[sapply(UTRusage, length) == 2]
    UTRusage <- unlist(UTRusage, use.names = FALSE)
    UTRusage.total <- UTRusage[UTRusage$source == "short"]
    UTRusage.long <- UTRusage[UTRusage$source == "long"]
    UTRusage.total <- UTRusage.total[!duplicated(UTRusage.total$transcript)]
    UTRusage.long <- UTRusage.long[match(
        UTRusage.total$transcript,
        UTRusage.long$transcript
    )]
    PDUItable <- UTRusage.total
    PDUItable$data <- NULL
    start(PDUItable)[as.character(strand(PDUItable)) == "-"] <-
        PDUItable$Predicted_Distal_APA[as.character(strand(PDUItable)) == "-"]
    end(PDUItable)[as.character(strand(PDUItable)) == "+"] <-
        PDUItable$Predicted_Distal_APA[as.character(strand(PDUItable)) == "+"]
    
    UTRusage.total.data <- do.call(rbind, UTRusage.total$data)
    UTRusage.long.data <- do.call(rbind, UTRusage.long$data)
    UTRusage.short.data <- UTRusage.total.data - UTRusage.long.data
    lt0 <- apply(UTRusage.short.data, 1, function(.ele) any(.ele < 0))
    if (any(lt0)) {
        CDSusage <- UTR3CDS.cov %>% 
            plyranges::filter(feature == "CDS" &
                            transcript %in% 
                                unique(PDUItable$transcript[lt0]))
        CDSusage.data <- do.call(rbind, CDSusage$data)
        rownames(CDSusage.data) <- CDSusage$transcript
        idx <- match(PDUItable$transcript[lt0], CDSusage$transcript)
        
        ## why this is reasonable?
        UTRusage.short.data[lt0[!is.na(idx)]] <-
            CDSusage.data[idx[!is.na(idx)], ] - 
            UTRusage.long.data[lt0[!is.na(idx)]]
        UTRusage.short.data[UTRusage.short.data < 0] <- 0
    }
    
    ## normalization?
    normalize.foo <- function(exprs, FUN) {
        avgs <- FUN(exprs)$Estimates
        scaling.factors <- avgs / avgs[1]
        scaling.factors <- matrix(rep(scaling.factors, nrow(exprs)),
                                  ncol = length(scaling.factors), byrow = TRUE
        )
        exprs <- exprs / scaling.factors
    }
    normalize.mean <- function(exprs) {
        normalize.foo(exprs, colSummarizeAvg)
    }
    normalize.median <- function(exprs) {
        normalize.foo(exprs, colSummarizeMedian)
    }
    if (normalize != "none") {
        UTRusage.long.short.data <- rbind(
            UTRusage.long.data,
            UTRusage.short.data
        )
        UTRusage.long.short.data <-
            switch(normalize,
                   quantile = normalize.quantiles(UTRusage.long.short.data),
                   quantile.robust = normalize.quantiles.robust(UTRusage.long.short.data, ...),
                   mean = normalize.mean(UTRusage.long.short.data),
                   median = normalize.median(UTRusage.long.short.data),
                   UTRusage.long.short.data
            ) ## suppose the output should be same order
        UTRusage.long.data <-
            UTRusage.long.short.data[1:nrow(UTRusage.long.data), , drop = FALSE]
        UTRusage.short.data <-
            UTRusage.long.short.data[-(1:nrow(UTRusage.long.data)), , drop = FALSE]
    }
    
    UTRusage.PDUI <- UTRusage.long.data / (UTRusage.long.data + UTRusage.short.data)
    UTRusage.PDUI.log2 <- log2(UTRusage.PDUI + .Machine$double.xmin)
    rownames(UTRusage.PDUI) <-
        rownames(UTRusage.PDUI.log2) <-
        rownames(UTRusage.long.data) <-
        rownames(UTRusage.short.data) <- PDUItable$transcript
    PDUItable$source <- NULL
    if (singleSample) {
        signals.short <- UTRusage.total$data2
        signals.long <- UTRusage.long$data2
        
        ## ??
        cut50 <- function(x, y) {
            z <- c(x, y)
            if (length(z) > 50) {
                xat <- floor(50 * length(x) / length(z))
                z <- tapply(z, cut(1:length(z), 50), mean)
            } else {
                xat <- length(x)
            }
            c(xat, z)
        }
        signals <- mapply(cut50, signals.short, signals.long, SIMPLIFY = FALSE)
        names(signals) <- UTRusage.total$transcript
        new("UTR3eSet",
            usage = PDUItable, PDUI = UTRusage.PDUI,
            PDUI.log2 = UTRusage.PDUI.log2,
            short = UTRusage.short.data,
            long = UTRusage.long.data,
            signals = signals
        )
    } else {
        new("UTR3eSet",
            usage = PDUItable, PDUI = UTRusage.PDUI,
            PDUI.log2 = UTRusage.PDUI.log2,
            short = UTRusage.short.data,
            long = UTRusage.long.data
        )
    }
}
