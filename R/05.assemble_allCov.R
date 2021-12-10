#' Assemble coverage files for a given chromosome for all samples
#' 
#' Process individual sample-chromosome-specific coverage files in an
#'   experiment into a file containing a list of chromosome-specific Rle
#'   coverage of all samples
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. 
#'   the output of setup_sqlitedb()
#' @param seqname A character(1) vector, the name of a chromosome/scaffold
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or 
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#'   
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
#'                             chr2exclude = "chrM",
#'                             BPPARAM = NULL)
#'    }
#'    chr_coverage <- assemble_allCov(sqlite_db,
#'                                    seqname = "chr6",
#'                                    outdir = outdir, 
#'                                    genome = genome, 
#'                                    chr2exclude = "chrM")
#' }

assemble_allCov <- function(sqlite_db,
                            seqname,
                            outdir = getInPASOutputDirectory(),
                            genome = getInPASGenome(), 
                            chr2exclude = getChr2Exclude()){
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }
    if (missing(sqlite_db) || !file.exists(sqlite_db) || 
        length(sqlite_db) != 1){
        stop("The sqlite_db length is not 1 or it doesn't exist!")
    }
    if (!is.null(chr2exclude) && !is.character(chr2exclude))
    {
        stop("chr2Exclude must be NULL or a character vector")
    }
    lock_filename <- getLockName()
    if (!file.exists(lock_filename)) 
    {
        stop("lock_filename must be an existing file.",
             "Please call addLockName() first!")
    }
    
    file_lock <- lock(lock_filename)
    db_conn <- dbConnect(drv = RSQLite::SQLite(), 
                         dbname = sqlite_db)
    res <- dbGetQuery(db_conn, "SELECT * FROM sample_coverage;") 
    dbDisconnect(db_conn)
    unlock(file_lock)
    
    if (nrow(res) < 1){
        stop("The sample_coverage table in the sqlite_db is empty")
    }
    if (!seqname %in% res$chr){
        stop("seqname", seqname, "is not in the sample_coverage table")
    }
    if (!is.character(outdir) || length(outdir) != 1){
        stop("A directory of write permission is required!")
    }
    outdir <- file.path(outdir, "003.chromosomewise.RleCov")
    if (!dir.exists(outdir)){
        dir.create(outdir, recursive = TRUE)
    }
    outdir <- normalizePath(outdir, mustWork = TRUE)

    seqLen <- get_seqLen(genome, chr2exclude)
    
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
    filename_chr <- file.path(outdir, 
                          paste(seqname, 
                                "RleCov.RDS", sep = "_"))
    ## chr_cov: chr_cov[[seqname]][[tag]]
    saveRDS(chr_cov, file = filename_chr)

        tryCatch({
        file_lock <- lock(lock_filename)
        db_conn <- dbConnect(drv = RSQLite::SQLite(), 
                             dbname= sqlite_db)
        res <- dbSendStatement(db_conn, 
                               paste0("DELETE FROM chromosome_coverage 
                                      WHERE chr = '", seqname, "';"))
        dbClearResult(res)
        res <- dbSendStatement(db_conn, 
                               paste0("INSERT INTO 
                        chromosome_coverage (chr, coverage_file) 
                        VALUES ('", seqname,"',", "'", filename_chr, "');"))
        dbClearResult(res)
    }, error = function(e) {
        print(paste(conditionMessage(e)))
    },
    finally = {dbDisconnect(db_conn)
        unlock(file_lock)})
    chr_cov
}