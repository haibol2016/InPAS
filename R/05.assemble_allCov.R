#' Assemble coverage files for all samples
#' 
#' Process individual sample-chromosome-specific coverage files in an
#'   experiment into a file containing a list of chromosome-specific Rle
#'   coverage of all samples
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. 
#'   the output of setup_sqlitedb()
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   the coverage data. If it doesn't exist, it will be created.
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param removeScaffolds A logical(1) vector, whether the scaffolds should be
#'   removed from the genome If you use a TxDb containing alternative
#'   scaffolds, you'd better to remove the scaffolds.

#' @return A list of paths to per-chromosome coverage files of all samples.
#' \itemize{\item seqname, chromosome/scaffold name
#'             \itemize{
#'                         \item  tag1, name tag for sample1
#'                         \item  tag2, name tag for sample2
#'                         \item  tagN, name tag for sampleN
#'                     }
#'         }
#' @export
#' @import RSQLite
#' @author Haibo Liu
#' @examples 
#' if (interactive()) {
#'    library(BSgenome.Mmusculus.UCSC.mm10)
#'    genome <- BSgenome.Mmusculus.UCSC.mm10
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
#'                                                     "metadata.txt"),
#'                                outdir)
#'    coverage <- list()
#'    for (i in seq_along(bedgraphs)){
#'    coverage[[tags[i]]] <- get_ssRleCov(bedgraph = bedgraphs[i],
#'                             tag = tags[i],
#'                             genome = genome,
#'                             sqlite_db = sqlite_db,
#'                             outdir = outdir,
#'                             removeScaffolds = TRUE,
#'                             BPPARAM = NULL)
#'    }
#'    coverage_files <- assemble_allCov(sqlite_db, 
#'                                     outdir, 
#'                                     genome, 
#'                                     removeScaffolds = FALSE)
#' }

assemble_allCov <- function(sqlite_db, 
                          outdir,
                          genome,
                          removeScaffolds = FALSE){
    if (missing(genome)) {
        stop("genome is required.")
    }
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }
    if (missing(sqlite_db) || !file.exists(sqlite_db) || 
        length(sqlite_db) != 1){
        stop("The sqlite_db length is not 1 or it doesn't exist!")
    }
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
    res <- dbGetQuery(db_conn, "SELECT * FROM sample_coverage;")
    dbDisconnect(db_conn)
    
    if (nrow(res) < 1){
        stop("The sample_coverage table in the sqlite_db is empty")
    }
    if (is.null(lock_filename) || !file.exists(lock_filename)) 
    {
        stop("lock_filename must be an existing file.")
    }
    if (missing(outdir)){
        stop("A directory of write permission is required!")
    }
    od <- file.path(outdir, "chromosmewise_RleCov")
    if (!dir.exists(od)){
        dir.create(od, recursive = TRUE)
    }
    od <- normalizePath(od, mustWork = TRUE)
    
    seqLen <- get_seqLen(genome, removeScaffolds)
    
    filenames <- list()
    for (seqname in unique(res$chr)){
        chr_cov <- list()
        for (tag in unique(res$tag)){
            file_name <- res$coverage_file[res$tag == tag & 
                                               res$chr == seqname]
            if (length(file_name) == 1) {
                chr_cov[[seqname]][[tag]] <- readRDS(file_name)
            } else {
                chr_cov[[seqname]][[tag]] <- 
                    Rle(values = 0, lengths = seqLen[seqname])
            }
        }
        filename <- file.path(od, 
                              paste(seqname, 
                                    "RleCov.RDS", sep = "_"))
        
        ## chr_cov: chr_cov[[seqname]][[tag]]
        saveRDS(chr_cov, file = filename)
        filenames[[seqname]] <- filename
    }

    filename_df <- data.frame(chr = names(filenames), 
                              coverage_file = unlist(filenames))
    
    tryCatch({
        db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
        res <- dbSendStatement(db_conn, 
                               paste0("DELETE FROM chromosome_coverage;"))
        dbClearResult(res)
        res <- dbSendStatement(db_conn, 
                               "INSERT INTO 
                        chromosome_coverage (chr, coverage_file) 
                        VALUES (:chr, :coverage_file);", 
                               filename_df)
        dbClearResult(res)
    }, error = function(e) {
        print(paste(conditionMessage(e)))
    },
    finally = {dbDisconnect(db_conn)})
    
    ## filenames[[seqname]]
    filenames
}