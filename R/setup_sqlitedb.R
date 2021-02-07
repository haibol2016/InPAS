#' Create an SQLite database for storing metadata and paths to coverage files
#' 
#' Create an SQLite database with three tables, metadata, sample_coverage, 
#'   and chromosome_coverage, for storing metadata (sample tag, condition, 
#'   paths to begraph files, and sample total read coverage), paths 
#'   to sample-oriented coverage files, and paths to chromosome-oriented 
#'   coverage file, respectively
#'   
#' @param metadata A path to a tab-delimited file, with columns "tag", "condition",
#'   and "bedgraph_file", storing a unique name tag for each sample, a condition 
#'   name for each sample, such as "treatment" and "control", and a path to the
#'   bedgraph file for each sample
#' @param outdir A character(1), a path with write permission for storing the
#'   SQLite database
#'   
#' @import RSQLite
#'
#' @export
#' 
#' @return A character(1), the path to the SQLite database
#'

setup_sqlitedb <- function(metadata, outdir){
    if (missing(metadata) || !file.exists(metadata)){
        stop("A path to a tab-delimited file, with columns: ",
             "tag, condition, and bedgraph_file, is required")
    }
    if(missing(outdir)){
        stop("A path with write permission for storing the SQLite database ",
             "is required")
    }
    
    outdir <- normalizePath(outdir)
    inpas <- file.path(outdir, "InPAS_hugeData.sqlite")
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= inpas)
    
    ## Drop existing tables
    tables <- c("metadata", "chromosome_coverage", "sample_coverage")
    for (table in tables){
        dbSendQuery(db_conn, paste("DROP TABLE IF EXISTS ", table,";"))
    }
    
    ## create table: metadata
    dbSendQuery(db_conn, "CREATE TABLE 
                metadata (tag TEXT PRIMARY KEY,
                condition TEXT,
                bedgraph_file TEXT,
                depth INTEGER);")
    
    metadata <- read.delim(metadata, header = TRUE)
    required_cols <- c("tag", "condition", "bedgraph_file")
    if (any(!required_cols %in% colnames(metadata))) {
        stop("Your metadata file does NOT contain the essential columns: ",
        "tag, condition, bedgraph_file")
    }
    metadata$depth <- 0
    dbWriteTable(db_conn, "metadata", metadata, 
                 overwrite = TRUE, row.names = FALSE)
    
    ## create table: sample_coverage
    dbSendQuery(db_conn, "CREATE TABLE 
                sample_coverage (
                tag TEX,
                chr TEXT,
                coverage_file TEXT,
                PRIMARY KEY(tag, chr));")
    
    ## create table: chromosome_coverage
    dbSendQuery(db_conn, "CREATE TABLE 
                chromosome_coverage (
                chr TEXT,
                coverage_file TEXT,
                PRIMARY KEY(chr));")

    if (!any(dbListTables(db_conn) %in% tables)) {
        stop("The SQLite database has NOT been sep up correctly!")
    }
    dbDisconnect(db_conn)
    inpas
}