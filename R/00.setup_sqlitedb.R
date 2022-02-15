#' Create an SQLite database for storing metadata and paths to coverage files
#'
#' Create an SQLite database with five tables, "metadata","sample_coverage",
#'   "chromosome_coverage", "CPsites", and "utr3_coverage", for storing metadata
#'   (sample tag, condition, paths to bedgraph files, and sample total read
#'   coverage), sample-then-chromosome-oriented coverage files (sample tag,
#'   chromosome, paths to bedgraph files for each chromosome), and paths to
#'   chromosome-then-sample-oriented coverage files (chromosome, paths to
#'   bedgraph files for each chromosome), CP sites on each chromosome (chromosome,
#'   paths to cpsite files), read coverage for 3' UTR and last CDS regions on each
#'   chromosome (chromosome, paths to utr3 coverage file), respectively
#'
#' @param metadata A path to a tab-delimited file, with columns "tag", "condition",
#'   and "bedgraph_file", storing a unique name tag for each sample, a condition
#'   name for each sample, such as "treatment" and "control", and a path to the
#'   bedgraph file for each sample
#' @param outdir A character(1) vector, a path with write permission for storing
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @import RSQLite
#' @importFrom utils read.delim
#' @export
#' @author Haibo Liu
#' @return A character(1) vector, the path to the SQLite database
#' @examples
#' if (interactive()) {
#'   bedgraphs <- system.file("extdata", c(
#'     "Baf3.extract.bedgraph",
#'     "UM15.extract.bedgraph"
#'   ),
#'   package = "InPAS"
#'   )
#'   tags <- c("Baf3", "UM15")
#'   metadata <- data.frame(
#'     tag = tags,
#'     condition = c("Baf3", "UM15"),
#'     bedgraph_file = bedgraphs
#'   )
#'   outdir <- tempdir()
#'   write.table(metadata,
#'     file = file.path(outdir, "metadata.txt"),
#'     sep = "\t", quote = FALSE, row.names = FALSE
#'   )
#'   sqlite_db <- setup_sqlitedb(
#'     metadata =
#'       file.path(outdir, "metadata.txt"),
#'     outdir
#'   )
#' }
setup_sqlitedb <- function(metadata,
                           outdir = getInPASOutputDirectory()) {
  if (missing(metadata) || !file.exists(metadata)) {
    stop(
      "A path to a tab-delimited file, with columns:",
      "tag, condition, and bedgraph_file, is required"
    )
  }
  if (!is.character(outdir) || length(outdir) != 1) {
    stop(
      "A path with write permission for storing",
      "the SQLite database is required"
    )
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  outdir <- normalizePath(outdir)
  sqlite_db <- file.path(outdir, "InPAS_hugeData.sqlite")
  db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)

  ## Drop existing tables
  tables <- c(
    "metadata", "genome_anno",
    "chromosome_coverage",
    "sample_coverage", "total_coverage",
    "utr3_total_coverage", "CPsites",
    "utr3cds_coverage"
  )
  for (table in tables) {
    res <- dbSendStatement(db_conn, paste("DROP TABLE IF EXISTS ", table, ";"))
    dbClearResult(res)
  }

  ## create table: metadata
  res <- dbSendStatement(db_conn, "CREATE TABLE
                       metadata (tag TEXT,
                       condition TEXT,
                       bedgraph_file TEXT,
                       depth INTEGER,
                       PRIMARY KEY (tag));")
  dbClearResult(res)

  metadata <- read.delim(metadata, header = TRUE)
  required_cols <- c("tag", "condition", "bedgraph_file")
  if (any(!required_cols %in% colnames(metadata))) {
    stop(
      "Your metadata file does NOT contain the essential columns: ",
      "tag, condition, bedgraph_file"
    )
  }
  metadata$depth <- 0
  dbWriteTable(db_conn, "metadata", metadata,
    overwrite = TRUE, row.names = FALSE
  )
  ## create table: genome_anno
  res <- dbSendStatement(db_conn, "CREATE TABLE
                genome_anno (
                type TEXT,
                anno_file TEXT,
                PRIMARY KEY(type));")
  dbClearResult(res)

  ## create table: sample_coverage
  res <- dbSendStatement(db_conn, "CREATE TABLE
                sample_coverage (
                tag TEXT,
                chr TEXT,
                coverage_file TEXT,
                PRIMARY KEY(tag, chr),
                FOREIGN KEY(tag) REFERENCES metadata(tag));")
  dbClearResult(res)

  ## create table: chromosome_coverage
  res <- dbSendStatement(db_conn, "CREATE TABLE
                chromosome_coverage (
                chr TEXT,
                coverage_file TEXT,
                PRIMARY KEY(chr),
                FOREIGN KEY(chr) REFERENCES sample_coverage(chr));")
  dbClearResult(res)

  ## create table: total_coverage for CP site search
  res <- dbSendStatement(db_conn, "CREATE TABLE
                total_coverage (
                chr TEXT,
                coverage_file TEXT,
                PRIMARY KEY (chr),
                FOREIGN KEY(chr) REFERENCES sample_coverage(chr));")
  dbClearResult(res)

  ## create table: utr3_total_coverage for CP site search
  res <- dbSendStatement(db_conn, "CREATE TABLE
                utr3_total_coverage (
                chr TEXT,
                coverage_file TEXT,
                PRIMARY KEY(chr),
                FOREIGN KEY(chr) REFERENCES sample_coverage(chr));")
  dbClearResult(res)

  ## create table: CPsites
  res <- dbSendStatement(db_conn, "CREATE TABLE
                CPsites (
                chr TEXT,
                cpsites_file TEXT,
                PRIMARY KEY(chr),
                FOREIGN KEY(chr) REFERENCES sample_coverage(chr));")
  dbClearResult(res)

  ## create table utr3cds_coverage
  res <- dbSendStatement(db_conn, "CREATE TABLE
                utr3cds_coverage (
                chr TEXT,
                coverage_file TEXT,
                PRIMARY KEY (chr),
                FOREIGN KEY(chr) REFERENCES sample_coverage(chr));")
  dbClearResult(res)

  if (!any(dbListTables(db_conn) %in% tables)) {
    stop("The SQLite database has NOT been sep up correctly!")
  }
  dbDisconnect(db_conn)
  options(InPAS.sqliteDb = sqlite_db)
  sqlite_db
}

#' Get the path to an SQLite database
#'
#' @return A path to an SQLite database
#' @export
getInPASSQLiteDb <- function() {
  sqlite_db <- options()$InPAS.sqliteDb
  if (!file.exists(sqlite_db)) {
    stop(
      "SQLite database doesn't exist!\n",
      "Please call setup_sqlitedb() first"
    )
  }
  sqlite_db
}
