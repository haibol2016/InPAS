#' prepare data for predicting cleavage and polyadenylation (CP) sites
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param chr.utr3 An object of [GenomicRanges::GRanges-class], an element of the
#'   output of [extract_UTR3Anno()]
#' @param seqname A character(1), the name of a chromosome/scaffold
#' @param background A character(1) vector, the range for calculating cutoff
#'   threshold of local background. It can be "same_as_long_coverage_threshold",
#'   "1K", "5K","10K", or "50K".
#' @param TxDb an object of [GenomicFeatures::TxDb-class]
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or 
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#' @param hugeData A logical(1) vector, indicating whether it is huge data
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param silence report progress or not. By default it doesn't report progress.
#' @param minZ A numeric(1), a Z score cutoff value
#' @param cutStart An integer(1) vector a numeric(1) vector. What percentage or
#'   how many nucleotides should be removed from 5' extremities before searching
#'   for CP sites? It can be a decimal between 0, and 1, or an integer greater 
#'   than 1. 0.1 means 10 percent, 25 means cut first 25 bases
#' @param MINSIZE A integer(1) vector, specifying the minimal length in bp of a
#'   short/proximal 3' UTR. Default, 10
#' @param coverage_threshold An integer(1) vector, specifying the cutoff 
#'   threshold of coverage for first 100 nucleotides. If the coverage of first 
#'   100 nucleotides is lower than coverage_threshold, that transcript will be
#'   not considered for further analysis. Default, 5.
#'
#' @return A list as described below:
#'   \describe{\item{background}{The type of methods for background
#'                    coverage calculation}
#'             \item{z2s}{Z-score cutoff thresholds for each 3' UTRs}
#'             \item{depth.weight}{A named vector containing depth weight}}
#'
#' @import S4Vectors Biobase GenomicRanges GenomicFeatures methods
#' @importFrom BSgenome getSeq matchPWM
#' @importFrom utils object.size
#' @importFrom plyranges as_granges complement_ranges disjoin_ranges filter
#'   group_by mutate reduce_ranges reduce_ranges_directed remove_names select
#'   set_genome_info shift_downstream summarise
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle mapSeqlevels seqlevels
#' @importFrom parallel detectCores
#' @export
#' @author Jianhong Ou, Haibo Liu
#' 
#' @examples
#'  if (interactive()) {
#'     library(BSgenome.Mmusculus.UCSC.mm10)
#'     library("TxDb.Mmusculus.UCSC.mm10.knownGene")
#'     genome <- BSgenome.Mmusculus.UCSC.mm10
#'     TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'
#'     ## load UTR3 annotation and convert it into a GRangesList
#'     data(utr3.mm10)
#'     utr3 <- split(utr3.mm10, seqnames(utr3.mm10), drop = TRUE)
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
#'                              chr2exclude = "chrM")
#'     }
#'     data4CPsitesSearch <- setup_CPsSearch(sqlite_db,
#'                                           genome,
#'                                           chr.utr3 = utr3[["chr6"]],
#'                                           seqname = "chr6",
#'                                           background = "10K",
#'                                           TxDb = TxDb,
#'                                           chr2exclude = "chrM",
#'                                           hugeData = TRUE,
#'                                           outdir = outdir)
#'  }

setup_CPsSearch <- function(sqlite_db,
                            genome = getInPASGenome(), 
                            chr.utr3,
                            seqname,
                            background = c(
                                "same_as_long_coverage_threshold",
                                "1K", "5K", "10K", "50K"),
                            TxDb = getInPASTxDb(),
                            chr2exclude = getChr2Exclude(),
                            hugeData = TRUE,
                            outdir = getInPASOutputDirectory(),
                            silence = FALSE,
                            minZ = 2,
                            cutStart = 10,
                            MINSIZE = 10,
                            coverage_threshold = 20) {
    gcCompensation = NA
    mappabilityCompensation = NA
    FFT = FALSE
    fft.sm.power = 20
    
    if (!is.null(chr2exclude) && !is.character(chr2exclude))
    {
        stop("chr2Exclude must be NULL or a character vector")
    }
    if (!is.character(outdir) || length(outdir) != 1){
        stop("An explicit output directory is required")
    } else {
        outdir_total <- file.path(outdir, "004.TotalCov.out")
        if (!dir.exists(outdir_total)){
            dir.create(outdir_total, recursive = TRUE, 
                       showWarnings = FALSE)
        }
        outdir_total <- normalizePath(outdir_total)
        
        outdir_utr3 <- file.path(outdir, "005.UTR3TotalCov.out")
        if (!dir.exists(outdir_utr3)){
            dir.create(outdir_utr3, recursive = TRUE, 
                       showWarnings = FALSE)
        }
        outdir_utr3 <- normalizePath(outdir_utr3)
    }
    
    if (missing(chr.utr3)) {
        stop("chr.utr3 are required.")
    }
    if (!is(chr.utr3, "GRanges") ||
        !all(chr.utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
        stop("utr3 must be one element of the output of function of ",
             "extract_UTR3Anno()")
    }
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }

    if (seqlevelsStyle(chr.utr3) != seqlevelsStyle(genome)) {
        stop("The seqlevelsStyle of utr3 must be same as genome")
    }
    if (missing(sqlite_db) || length(sqlite_db) != 1 ||
        !file.exists(sqlite_db)) {
        stop("sqlite_db, a path to the SQLite database is required!")
    }
    lock_filename <- getLockName()
    if (!file.exists(lock_filename)) 
    {
        stop("lock_filename must be an existing file.",
             "Please call addLockName() first!")
    }

    background <- match.arg(arg = background, 
                            choices = c("same_as_long_coverage_threshold",
                                      "1K", "5K", "10K", "50K"))
    introns <- GRanges()
    if (background != "same_as_long_coverage_threshold") {
        if (missing(TxDb) || !is(TxDb, "TxDb")) {
            stop("TxDb is missing when you want local background")
        }
        if (!seqname %in% seqlevels(TxDb)){
            stop("seqname", seqname, "is not included in TxDb")
        }
        seqlevels(TxDb) <- seqname
        introns <- unlist(intronsByTranscript(TxDb, use.names = TRUE)) %>%
            plyranges::disjoin_ranges() ## strand-unaware
        exons <- exons(TxDb) %>% plyranges::disjoin_ranges()
        
        intron_exons <- disjoin(c(exons, introns))
        ol.exons <- findOverlaps(exons, intron_exons)
        
        ## pure intronic regions
        introns <- intron_exons[-unique(subjectHits(ol.exons))]
    }
    
    chr.utr3 <- chr.utr3[chr.utr3$feature != "CDS"]
    if (length(chr.utr3) == 0){return(NULL)}
    
    if (!silence) message("coverage per sample per chromosome start at ", 
                          date(), ".\n")
    ## assemble coverage data into chromosome-oriented list of Rle objects
    chr.cov <- assemble_allCov(sqlite_db = sqlite_db,
                               seqname = seqname,
                               outdir = outdir,
                               genome = genome, 
                               chr2exclude = getChr2Exclude())
    if (!silence) message("coverage per sample per chromosome done at ", 
                          date(), ".\n")
    
    ## calculate coverage weight
    file_lock <- flock::lock(lock_filename)
    db_conn <- dbConnect(drv = RSQLite::SQLite(), 
                                            dbname = sqlite_db)
    metadata <- dbReadTable(db_conn, "metadata") 
    dbDisconnect(db_conn)
    flock::unlock(file_lock)
    
    depth.weight <- get_depthWeight(metadata = metadata, 
                                    hugeData = hugeData)
    
    if (!silence) message("total coverage start at ", date(), ".\n")
    ## get pooled coverage per condition per chromosome
    chr.totalCov <- get_totalCov(outdir = outdir_total,
                                 sqlite_db = sqlite_db,
                                 chr.cov = chr.cov, 
                                 seqname = seqname, 
                                 metadata = metadata, 
                                 hugeData = hugeData)
    if (!silence) message("total coverage done at ", date(), ".\n")
    
    ## calculate background
    z2s <- get_zScoreCutoff(background = background, 
                            chr.introns = introns, 
                            chr.totalCov = chr.totalCov, 
                            chr.utr3 = chr.utr3, 
                            seqname = seqname, 
                            z = minZ)
    if (!silence) message("backgroud around 3utr done at ", date(), ".\n")
    
    ## get UTR3 total coverage
    chr.cov <- get_UTR3TotalCov(chr.utr3 = chr.utr3,
                                chr.totalCov = chr.totalCov,
                                gcCompensation = gcCompensation,
                                mappabilityCompensation =
                                    mappabilityCompensation,
                                FFT = FFT,
                                fft.sm.power = fft.sm.power)
    if (!silence) message("utr3 TotalCov done at ", date(), ".\n")
    
    chr.cov <- chr.cov[sapply(chr.cov, mean) > 0]
    if (length(chr.cov) == 0) {
        return(NULL)
    }
    utr3.utr <- chr.utr3[chr.utr3$feature == "utr3"]
    utr3.gap <- chr.utr3[chr.utr3$feature == "next.exon.gap"]
    co <- countOverlaps(utr3.gap, utr3.utr, maxgap = 1, 
                        ignore.strand = TRUE)
    utr3.gap <- utr3.gap[co > 1]
    chr.utr3$conn_next_utr3 <-
        chr.utr3$transcript %in% utr3.gap$transcript
    
    ## remove utr3 Granges without coverage
    chr.utr3 <- chr.utr3[names(chr.utr3) %in% names(chr.cov)]
    if (length(chr.utr3) == 0){
        return(NULL)
    }
    
    chr.utr3 <- split(chr.utr3, chr.utr3$transcript)
    conn_next_utr3 <- sapply(chr.utr3, function(.UTR) {
        .UTR$conn_next_utr3[1]
    })
    
    merge_chrCov <- function(.UTR) {
        .UTR <- .UTR[order(start(.UTR))]
        chr.utr3TotalCov <- chr.cov[names(.UTR)]
        chr.utr3TotalCov <-
            mapply(function(.covList, .start, .end, .property) {
                # set names for each position
                .posList <- .start:.end
                
                # if not a matrix
                if (length(dim(.covList)) == 0 && !is.null(.covList)) {
                    .covList <- t(.covList)
                } 
                rownames(.covList) <- paste(.property, .posList, sep = "_SEP_")
                .covList
            }, chr.utr3TotalCov, start(.UTR), end(.UTR), 
            .UTR$feature, SIMPLIFY = FALSE)
        
        chr.utr3TotalCov <- do.call(rbind, chr.utr3TotalCov)
        ## reverse the negative strand
        if (as.character(strand(.UTR))[1] == "-") { 
            chr.utr3TotalCov <-
                chr.utr3TotalCov[rev(rownames(chr.utr3TotalCov)), ,
                                 drop = FALSE]
        }
        
        # if the range of "cutstart" is (0, 1), percentage; 
        # otherwise, absolute bases
        if (!is.na(cutStart)) {
            if (cutStart < 1) {
                cutStart <- floor(length(chr.utr3TotalCov) * cutStart)
            }
            if (cutStart > 0) {
                chr.utr3TotalCov <- chr.utr3TotalCov[-(1:cutStart), ,
                                                     drop = FALSE]
            }
        }
        chr.utr3TotalCov
    }
    
    chr.cov.merge <- lapply(chr.utr3, merge_chrCov)
    if (!silence) message("chromsome ", seqname, " coverage merged.\n")
    
    ## filter UTR3 based on coverage of the first 100 bases
    coverage_quality <- sapply(chr.cov.merge, function(.ele) {
        if (nrow(.ele[grepl("utr3_SEP_", rownames(.ele)), ,
                      drop = FALSE]) > MINSIZE) {
            any(colMeans(.ele[1:min(nrow(.ele), 100), ,
                              drop = FALSE]) > coverage_threshold)
        } else {
            FALSE
        }})
    
    chr.cov.merge <- chr.cov.merge[coverage_quality]
    conn_next_utr3 <- conn_next_utr3[coverage_quality]
    chr.utr3 <- chr.utr3[coverage_quality]
    
    if (!silence) message("Preparation for CPsite search done at ", 
                          date(), ".\n")
    
    if (length(chr.cov.merge) > 0) {
        CPsSearch_data <- list(background = background,
                               z2s = z2s, 
                               depth.weight = depth.weight,
                               chr.cov.merge = chr.cov.merge,
                               conn_next_utr3 = conn_next_utr3,
                               chr.utr3 = chr.utr3)
        tmpfile <- file.path(outdir_utr3, 
                             paste(seqname, "data4CPsSearch.RDS",
                                   sep = "_"))
        saveRDS(CPsSearch_data, file = tmpfile)
        
        tryCatch({
            file_lock <- flock::lock(lock_filename)
            db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
            res <- dbSendStatement(db_conn,
                                   paste0("DELETE FROM utr3_total_coverage 
                                      WHERE chr = '", seqname, "';"))
            dbClearResult(res)
            ## insert into database the per-chromosome utr3_total_coverage
            res <- dbSendStatement(db_conn, 
                                   paste0("INSERT INTO 
                           utr3_total_coverage (chr, coverage_file) 
                           VALUES ('", seqname, "',", "'", tmpfile, "');"))
            dbClearResult(res)
        }, error = function(e) {
            print(paste(conditionMessage(e)))
        },finally = {dbDisconnect(db_conn)
            flock::unlock(file_lock)})
        return(CPsSearch_data)
    } else {
        return(NULL)
    }
}