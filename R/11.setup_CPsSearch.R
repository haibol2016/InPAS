#' prepare data for predicting cleavage and polyadenylation (CP) sites
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param utr3 An object of [GenomicRanges::GRangesList-class], output of
#'   [extract_UTR3Anno()]
#' @param background A character(1) vector, the range for calculating cutoff
#'   threshold of local background. It can be "same_as_long_coverage_threshold",
#'   "1K", "5K","10K", or "50K".
#' @param TxDb an object of [GenomicFeatures::TxDb-class]
#' @param removeScaffolds A logical(1) vector, whether the scaffolds should be
#'   removed from the genome If you use a TxDb containing alternative
#'   scaffolds, you'd better to remove the scaffolds.
#' @param hugeData A logical(1) vector, indicating whether it is huge data
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   the coverage data. If it doesn't exist, it will be created.
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam()
#'   
#' @param silence report progress or not. By default it doesn't report progress.
#'
#' @return A list as described below:
#'   \describe{\item{utr3TotalCov}{chromosome-wise 3' UTR coverage in summarized
#'                  View format}
#'              \describe{\item{chr1}{A filename for chr1 3' UTR coverage in
#'              summarized View format}
#'                        \item{chr2}{A filename for chr2 3' UTR coverage in
#'              summarized View format}
#'                        \item{chrN}{A filename for chrN 3' UTR coverage in
#'              summarized View format}}
#'              \item{background}{The type of methods for bckground
#'                    coverage calculation}
#'              \item{z2s}{Z-score cutoff thresholds for each 3' UTRs}
#'              \item{depth.weight}{A named vector containing depth weight}}
#'
#' @import S4Vectors Biobase GenomicRanges GenomicFeatures methods
#' @importFrom BSgenome getSeq matchPWM
#' @importFrom utils object.size
#' @importFrom plyranges as_granges complement_ranges disjoin_ranges filter
#'   group_by mutate reduce_ranges reduce_ranges_directed remove_names select
#'   set_genome_info shift_downstream summarise
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle mapSeqlevels seqlevels
#' @importFrom BiocParallel bptry bpok bplapply bpparam
#' @export
#' @author Jianhong Ou, Haibo Liu
#' 
#'  @examples
#'  if (interactive()) {
#'     library(BSgenome.Mmusculus.UCSC.mm10)
#'     library("TxDb.Mmusculus.UCSC.mm10.knownGene")
#'
#'     genome <- BSgenome.Mmusculus.UCSC.mm10
#'     TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'
#'     ## load UTR3 annotation and convert it into a GRangesList
#'     data(utr3.mm10)
#'     utr3 <- split(utr3.mm10, seqnames(utr3.mm10))
#'
#'     bedgraphs <- system.file("extdata",c("Baf3.extract.bedgraph",
#'                                          "UM15.extract.bedgraph"),
#'                             package = "InPAS")
#'     tags <- c("Baf3", "UM15")
#'     metadata <- data.frame(tag = tags,
#'                            condition = c("Baf3", "UM15"),
#'                            bedgraph_file = bedgraphs)
#'     outdir = tempdir()
#'     write.table(metadata, file =file.path(outdir, "metadata.txt"),
#'                 sep = "\t", quote = FALSE, row.names = FALSE)
#'
#'     sqlite_db <- setup_sqlitedb(metadata = file.path(outdir,
#'                                                      "metadata.txt"),
#'                                 outdir)
#'     coverage <- list()
#'     for (i in seq_along(bedgraphs)){
#'     coverage[[tags[i]]] <- get_ssRleCov(bedgraph = bedgraphs[i],
#'                              tag = tags[i],
#'                              genome = genome,
#'                              sqlite_db = sqlite_db,
#'                              outdir = outdir,
#'                              removeScaffolds = TRUE)
#'     }
#'     coverage_files <- assemble_allCov(sqlite_db,
#'                                      outdir,
#'                                      genome,
#'                                      removeScaffolds = TRUE)
#'     data4CPsitesSearch <- setup_CPsSearch(sqlite_db,
#'                                           genome,
#'                                           utr3,
#'                                          background = "10K",
#'                                          TxDb = TxDb,
#'                                          removeScaffolds = TRUE,
#'                                          hugeData = TRUE,
#'                                          outdir = outdir)
#'  }
#'


setup_CPsSearch <- function(sqlite_db,
                            genome, 
                            utr3,
                            background = c(
                                "same_as_long_coverage_threshold",
                                "1K", "5K", "10K", "50K"),
                            TxDb = NA,
                            removeScaffolds = FALSE,
                            hugeData = TRUE,
                            outdir,
                            BPPARAM = NULL,
                            silence = FALSE) {
    gcCompensation <- NA
    mappabilityCompensation <- NA
    FFT <- FALSE
    fft.sm.power <- 20
    MINSIZE <- 10

    if (missing(outdir)){
        stop("An explicit output directory is required")
    } else {
        outdir <- file.path(outdir, "CPsites.out")
        if (!dir.exists(outdir)){
            dir.create(outdir, recursive = TRUE, 
                       showWarnings = FALSE)
        }
        outdir <- normalizePath(outdir)
    }
    
    if (missing(genome) || missing(utr3)) {
        stop("genome and utr3 are required.")
    }
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }
    if (!is(utr3, "GRangesList") ||
        !all(utr3[[1]]$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
        stop("utr3 must be output of function of extract_UTR3Anno()")
    }
    if (seqlevelsStyle(utr3[[1]]) != seqlevelsStyle(genome)) {
        stop("the seqlevelsStyle of utr3 must be same as genome")
    }
    if (missing(sqlite_db)){
        stop("The path to the sQLite database is required")
    }

    background <- match.arg(arg = background, 
                            choices = c("same_as_long_coverage_threshold",
                                      "1K", "5K", "10K", "50K"))
    introns <- GRanges()
    if (background != "same_as_long_coverage_threshold") {
        if (!is(TxDb, "TxDb")) {
            stop("TxDb is missing when you want local background")
        }
        # if (!identical(unname(genome(TxDb))[1], unname(genome(utr3))[1])) {
        #     stop("genome of TxDb should be same as genome of utr3")
        # }
        if (removeScaffolds) {
            TxDb_seqlevels <- seqlevels(TxDb)
            seqlevels(TxDb) <- 
                TxDb_seqlevels[grepl("^(chr)?(\\d+|[XY])$", TxDb_seqlevels)]
        }
        introns <- unlist(intronsByTranscript(TxDb, use.names = TRUE)) %>%
            plyranges::disjoin_ranges() ## strand-unaware
        exons <- exons(TxDb) %>% plyranges::disjoin_ranges()
        
        intron_exons <- disjoin(c(exons, introns))
        ol.exons <- findOverlaps(exons, intron_exons)
        
        ## pure intronic regions
        introns <- intron_exons[-unique(subjectHits(ol.exons))]
    }
    if (!silence) message("total backgroud ... done.\n")
    
    ## utr3 is GRangesList
    utr3 <- endoapply(utr3, function(.x){
        .x <- .x[.x$feature != "CDS"]
    })
    
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
    metadata <- dbReadTable(db_conn, "metadata")
    chromosome_cov <- dbReadTable(db_conn, "chromosome_coverage")
    dbDisconnect(db_conn)
    
    depth.weight <- get_depthWeight(sqlite_db, hugeData = hugeData)
    
    if (nrow(chromosome_cov) < 1) {
        stop("chromosome_coverage table in the SQLite database is empty")
    }
    totalCov <- lapply(chromosome_cov$chr, function(.x){
        get_totalCov(sqlite_db = sqlite_db, 
                      seqname = .x, hugeData = hugeData)[[1]]
    })
    names(totalCov) <- chromosome_cov$chr
    total_coverage_filename <- file.path(outdir, "TotalCov.RDS")
    saveRDS(totalCov, file = total_coverage_filename)
    
    if (hugeData) { hugeData <- FALSE }

    if (!silence) message("total coverage ... done.\n")
    
    z2s <- get_zScoreCutoff(background, introns, totalCov, utr3)
    if (!silence) message("backgroud around 3utr ... done.\n")
    
    utr3TotalCov <- get_UTR3TotalCov(utr3, totalCov,
                                     BPPARAM = BPPARAM,
                                     gcCompensation,
                                     mappabilityCompensation,
                                     FFT = FFT,
                                     fft.sm.power = fft.sm.power)
    
    ## split utr3TotalCov by chromosome and save individual coverage
    utr3TotalCov <- mapply(function(.ele, .name) {
        utr3TotalCov.tmpfile <- file.path(outdir, 
                                    paste(.name, "utr3TotalCov.RDS",sep = "_"))
        saveRDS(.ele, utr3TotalCov.tmpfile, compress = "xz")
        if (!silence) {
            message("save utr3TotalCov for ", .name, " at ", 
                    utr3TotalCov.tmpfile, " ... done.\n")
        }
        utr3TotalCov.tmpfile
    }, utr3TotalCov, names(utr3TotalCov), SIMPLIFY = FALSE)
    
    filename_df <- data.frame(chr = names(utr3TotalCov),
                              coverage_file = unlist(utr3TotalCov))
    insist_execute_sqlite <- fix_dbLockError()
    
    tryCatch({db_conn <- dbConnect(drv = RSQLite::SQLite(), 
                                   dbname = sqlite_db)
    ## remove existing record in utr3_total_coverage
    res <- dbSendStatement(db_conn,
                           paste0("DELETE FROM utr3_total_coverage;"))
    dbClearResult(res)

    ## insert into database the per-chromosome utr3_total_coverage
    res <- dbSendStatement(db_conn, 
                           "INSERT INTO 
                           utr3_total_coverage (chr, coverage_file) 
                           VALUES (:chr, :coverage_file);", 
                           filename_df)
    dbClearResult(res)

    ## remove existing record in total_coverage
    res <- dbSendStatement(db_conn,
                           paste0("DELETE FROM total_coverage;"))
    dbClearResult(res)

    ## insert into total_coverage
    res <- dbSendStatement(db_conn, 
                           paste0("INSERT INTO 
                           total_coverage (coverage_file) 
                           VALUES ('", total_coverage_filename, "');"))
    dbClearResult(res)
    }, error = function(e) {
        print(paste(conditionMessage(e)))
    },finally = {dbDisconnect(db_conn)})
            
    list(utr3TotalCov = utr3TotalCov, 
         background = background,
         z2s = z2s, 
         depth.weight = depth.weight)
}