#' Get coverage for 3' UTR and last CDS regions on a single chromosome
#'
#' @param chr.utr3 one element of an output of [extract_UTR3Anno()]
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param outdir A path to a folder for storing coverage data of 3' UTRs and 
#'   last CDSs on a given chromosome/scaffold. If it doesn't exist, it will be
#'   created.
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam()
#' @param phmm A logical(1) vector, indicating whether data should be 
#'   prepared for singleSample analysis? By default, FALSE
#'
#' @return coverage view in GRanges
#' @export
#' @author Jianhong Ou, Haibo Liu

get_regionCov <- function(chr.utr3,
                          sqlite_db,
                          outdir,
                          BPPARAM = NULL,
                          phmm = FALSE){
  if (!is(chr.utr3, "GRanges") ||
      !all(chr.utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
    stop("chr.utr3 must be one element of an output of function of",
         "extract_UTR3Anno()")
  }
  if (missing(sqlite_db)){
      stop("sqlite_db, a path to the SQLite database is required")
  }
 if (missing(outdir)){
    stop("An explicit output directory is required")
  } else {
    outdir <- file.path(outdir, "UTR3_lastCDS.RleCov")
    if (!dir.exists(outdir)){
      dir.create(outdir, recursive = TRUE, 
                 showWarnings = FALSE)
    }
    outdir <- normalizePath(outdir)
  }
  
  ## get 3'UTRs and their last CDS regions
  chr.regions <- get_UTR3CDS(sqlite_db,
                             chr.utr3,
                             BPPARAM = BPPARAM)
  if (length(chr.regions) == 0) {return(NULL)}
  
  chr <- unique(as.character(seqnames(chr.regions)))
  if (length(chr) != 1){
    stop("chr.regions should only be for a single chromosome/scaffold")
  }
  
  end <- end(chr.regions)
  maxEnd <- max(end)
  
  db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
  chr.coverage <- dbReadTable(db_conn, "chromosome_coverage")
  dbDisconnect(db_conn)
  chr.coverage <- 
    chr.coverage$coverage_file[chr.coverage$chr == chr]
  if (length(chr.coverage) != 1){
    stop(paste0("No coverage data for ", chr))
  }
  ## load coverage of each sample for chr
  chr.coverage <- readRDS(chr.coverage)[[1]]

  .cov <- list()
  if (phmm) {all.tx <- list()}
  for (i in seq_along(chr.coverage)) {
      .ele <- chr.coverage[[i]]
      
      ## This should not happen for my implementation
      if (maxEnd > length(.ele)) {
          .ele <- append(.ele, rep(0, maxEnd - length(.ele) + 1))
      }
      if (is(.ele, "Rle")) {
          .cvg <- Views(.ele, start(chr.regions), end(chr.regions))
          if (phmm) {
              all.tx[[names(chr.coverage)[i]]] <- 
                  viewApply(.cvg, as.integer)
          }
          .cvg <- viewMeans(.cvg, na.rm = TRUE)
      } else {
          if (sum(.ele) != 0) {
              stop("Sum of current chromosome, which is not covered,",
              "is not zero.")
          }
          .cvg <- rep(0, length(chr.regions))
      }
      .cov[[names(chr.coverage)[i]]] <- .cvg
  }
      
  this.trans <- list()
  for (i in 1:length(.cov[[1]])) {
      this.trans[[i]] <- list()
      for (j in names(.cov)) {
          this.trans[[i]][[j]] <- .cov[[j]][[i]]
      }
      this.trans[[i]] <- do.call(cbind, this.trans[[i]])
  }
  chr.regions$data <- this.trans
  if (phmm) {chr.regions$data2 <- all.tx[[1]]}
  
  filename <- file.path(outdir, paste(chr, 
                                      "utr3_CDS.RleCov.RDS",
                                      sep = "_"))
  saveRDS(chr.regions, file = filename)
  
  insist_execute_sqlite <- fix_dbLockError()
  tryCatch({
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
    
    ## remove existing records for this chr
    insist_execute_sqlite(db_conn, "BEGIN IMMEDIATE TRANSACTION")
    # Rollback on failure
    on.exit(try(dbExecute(db_conn, "ROLLBACK TRANSACTION")))
    res <- dbSendStatement(db_conn, 
                           paste0("DELETE FROM   
                                utr3cds_coverage WHERE chr = '", chr, "';"))
    dbClearResult(res)
    dbExecute(db_conn, "COMMIT TRANSACTION")
    
    ## insert new data
    insist_execute_sqlite(db_conn, "BEGIN IMMEDIATE TRANSACTION")
    # Rollback on failure
    on.exit(try(dbExecute(db_conn, "ROLLBACK TRANSACTION")))
    
    res <- dbSendStatement(db_conn, 
                           paste0("INSERT INTO 
                         utr3cds_coverage (chr, coverage_file) 
                         VALUES ('", chr, "','", filename,"');"))
    dbClearResult(res)
    dbExecute(db_conn, "COMMIT TRANSACTION")
    
    # Don't Roll back on success
    on.exit(NULL)
  }, error = function(e) {
    print(paste(conditionMessage(e)))
  },
  finally = dbDisconnect(db_conn))
  chr.regions
}
