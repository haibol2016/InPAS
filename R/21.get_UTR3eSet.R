#' prepare 3' UTR coverage data for usage test
#'
#' generate a UTR3eSet object with PDUI information for statistic tests
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param normalize A character(1) vector, spcifying the normalization method. 
#'   It can be "none", "quantiles", "quantiles.robust", "mean", or "median"
#' @param ... parameter can be passed into
#'   [preprocessCore::normalize.quantiles.robust()]
#' @param singleSample A logical(1) vector, indicating whether data is prepared
#'   for analysis in a singleSample mode? Default, FALSE
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
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples
#' if (interactive()) {
#'    library(BSgenome.Mmusculus.UCSC.mm10)
#'    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#'    genome <- BSgenome.Mmusculus.UCSC.mm10
#'    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'    
#'    ## load UTR3 annotation and convert it into a GRangesList
#'    data(utr3.mm10)
#'    utr3 <- split(utr3.mm10, seqnames(utr3.mm10), drop = TRUE)
#'    
#'    bedgraphs <- system.file("extdata",c("Baf3.extract.bedgraph",
#'                                         "UM15.extract.bedgraph"), 
#'                            package = "InPAS")
#'    tags <- c("Baf3", "UM15")
#'    metadata <- data.frame(tag = tags, 
#'                           condition = c("Baf3", "UM15"),
#'                           bedgraph_file = bedgraphs)
#'    outdir = tempdir()
#'    write.table(metadata, file =file.path(outdir, "metadata.txt"), 
#'                sep = "\t", quote = FALSE, row.names = FALSE)
#'    
#'    sqlite_db <- setup_sqlitedb(metadata = file.path(outdir, 
#'                                          "metadata.txt"), outdir)
#'    coverage <- list()
#'    for (i in seq_along(bedgraphs)) {
#'    coverage[[tags[i]]] <- get_ssRleCov(bedgraph = bedgraphs[i],
#'                             tag = tags[i],
#'                             genome = genome,
#'                             sqlite_db = sqlite_db,
#'                             outdir = outdir,
#'                             chr2exclude = "chrM",
#'                             BPPARAM = NULL)}
#'                             
#'    data4CPsSearch <- setup_CPsSearch(sqlite_db,
#'                                      genome,
#'                                      chr.utr3 = utr3[["chr6"]],
#'                                      seqname = "chr6",
#'                                      background = "10K",
#'                                      TxDb = TxDb,
#'                                      chr2exclude = "chrM",
#'                                      hugeData = TRUE,
#'                                      outdir = outdir,
#'                                      minZ = 2,
#'                                      cutStart = 10,
#'                                      MINSIZE = 10,
#'                                      coverage_threshold = 5)
#'    ## polyA_PWM
#'    load(system.file("extdata", "polyA.rda", package = "InPAS"))
#'    
#'    ## load the Naive Bayes classifier model from the cleanUpdTSeq package
#'    library(cleanUpdTSeq)
#'    data(classifier)
#'    
#'    CPs <- search_CPs(seqname = "chr6",
#'                      sqlite_db = sqlite_db, 
#'                      chr.utr3 = utr3[["chr6"]],
#'                      genome = genome, 
#'                      MINSIZE = 10, 
#'                      window_size = 100,
#'                      search_point_START =50,
#'                      search_point_END = NA,
#'                      cutEnd = 0,
#'                      adjust_distal_polyA_end = TRUE,
#'                      long_coverage_threshold = 2,
#'                      PolyA_PWM = pwm, 
#'                      classifier = classifier,
#'                      classifier_cutoff = 0.8,
#'                      shift_range = 100,
#'                      step = 5,
#'                      two_way = FALSE,
#'                      outdir = outdir)
#' utr3_cds_cov <- get_regionCov(chr.utr3 = utr3[["chr6"]],
#'                               sqlite_db,
#'                               outdir,
#'                               phmm = FALSE)
#' eSet <- get_UTR3eSet(sqlite_db,
#'                      normalize ="none", 
#'                      singleSample = FALSE)
#' test_out <- test_dPDUI(eset = eSet, 
#'                        method = "fisher.exact",
#'                        normalize = "none",
#'                        sqlite_db = sqlite_db) 
#' }


get_UTR3eSet <- function(sqlite_db,
                        normalize = c("none", "quantiles",
                                      "quantiles.robust",
                                      "mean", "median"),
                        ...,
                        singleSample = FALSE) {
    if (missing(sqlite_db)|| length(sqlite_db) != 1 || !file.exists(sqlite_db)){
        stop("sqlite_db, a path to the SQLite database is required!")
    }
    if (!is.logical(singleSample) || length(singleSample) != 1)
    {
        stop("singleSample must be a logical(1) vector")
    }
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
    utr3_cov <- dbReadTable(db_conn, "utr3cds_coverage")
    dbDisconnect(db_conn)
    if (nrow(utr3_cov) < 1){
        stop("Coverage data for UTR3 and last CDSs is not available")
    }
    
    ## load the coverage object for UTR3 and last CDSs
    UTR3_CDS.cov <- lapply(utr3_cov$coverage_file, function(.x){
        readRDS(.x)
    })
    UTR3_CDS.cov <- unlist(GRangesList(UTR3_CDS.cov))
    
    normalize <- match.arg(arg = normalize, 
                           choices = c("none", "quantiles",
                                       "quantiles.robust",
                                       "mean", "median"))
    
    if ((!singleSample) && "data2" %in% colnames(mcols(UTR3_CDS.cov))) {
        message("Single sample mode is on")
        singleSample <- TRUE
    }
    UTRusage <- UTR3_CDS.cov %>% plyranges::filter(!is.na(source))
    UTRusage <- split(UTRusage, UTRusage$transcript)
    UTRusage <- UTRusage[sapply(UTRusage, length) == 2]
    UTRusage <- unlist(UTRusage, use.names = FALSE)
    UTRusage.total <- UTRusage[UTRusage$source == "short"]
    UTRusage.long <- UTRusage[UTRusage$source == "long"]
    UTRusage.total <- UTRusage.total[!duplicated(UTRusage.total$transcript)]
    UTRusage.long <- UTRusage.long[match(UTRusage.total$transcript,
                                         UTRusage.long$transcript)]
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
        CDSusage <- UTR3_CDS.cov %>% 
            plyranges::filter(feature == "CDS" &
                             transcript %in% unique(PDUItable$transcript[lt0]))
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
            UTRusage.short.data)
        UTRusage.long.short.data <-
            switch(normalize,
                   quantile = normalize.quantiles(UTRusage.long.short.data),
                   quantile.robust = 
                      normalize.quantiles.robust(UTRusage.long.short.data, ...),
                   mean = normalize.mean(UTRusage.long.short.data),
                   median = normalize.median(UTRusage.long.short.data),
                   UTRusage.long.short.data) 
        ## suppose the output should be same order
        UTRusage.long.data <-
            UTRusage.long.short.data[1:nrow(UTRusage.long.data), , drop = FALSE]
        UTRusage.short.data <-
            UTRusage.long.short.data[-(1:nrow(UTRusage.long.data)), ,
                                     drop = FALSE]
    }
    
    UTRusage.PDUI <- 
        UTRusage.long.data / (UTRusage.long.data + UTRusage.short.data)
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
        signals <- mapply(cut50, signals.short, 
                          signals.long, SIMPLIFY = FALSE)
        names(signals) <- UTRusage.total$transcript
        new("UTR3eSet",
            usage = PDUItable, 
            PDUI = UTRusage.PDUI,
            PDUI.log2 = UTRusage.PDUI.log2,
            short = UTRusage.short.data,
            long = UTRusage.long.data,
            signals = signals)
    } else {new("UTR3eSet",
                usage = PDUItable, 
                PDUI = UTRusage.PDUI,
                PDUI.log2 = UTRusage.PDUI.log2,
                short = UTRusage.short.data,
                long = UTRusage.long.data)
    }
}
