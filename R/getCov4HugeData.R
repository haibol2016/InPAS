#' Process coverage files for hugData
#' 
#' Process per-sample coverage files of an experiment into a list
#'   per chromosome coverage file and depth per sample
#'
#' @param coverageFiles A vector of character(n), paths to R objects storing
#'   coverage filenames and depth for each sample, which are output of
#'   bedgraph2RleCov() of Bam2RleCov() with hugeData set to "TRUE"
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   setup_sqlitedb()
#' @param outdir A path to a directory/folder where output is written

#' @return A list of coverage with tags as names
#' @export
#' @import RSQLite


getCov4HugeData <- function(coverageFiles,
                            sqlite_db, 
                            outdir){
    if (missing(sqlite_db) || !file.exists(sqlite_db) || 
        length(sqlite_db) != 1){
        stop("The sqlite_db length is not 1 or it doesn't exist!")
    }
    
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
    res <- dbGetQuery(db_conn, "SELECT * FROM sample_coverage;")
    if (nrow(res) < 1 && missing(coverageFiles)){
        stop("coverageFiles is required because the ",
        "sample_coverage table in the sqlite_db is empty")
    } 
    
    if (missing(outdir)){
        stop("A directory of write permission is required!")
    }
    
    od <- file.path(outdir, "chromosmewise_RleCov")
    if (!dir.exists(od)){
        dir.create(od, recursive = TRUE)
    }
    
    filenames <- list()
    
    ## extract sample read depth and rearrange sample-wise coverage filenames 
    ## into chromosome-wise coverage filenames
    if (nrow(res) < 1){
        depth <- vector("integer", length(coverageFiles))
        cov_files <- list(list())
        
        for (i in seq_along(coverageFiles)) {
            cov_depth <- readRDS(coverageFiles[i])
            depth[i] <- cov_depth$depth
            names(depth)[i] <- names(cov_depth$depth)
            
            tag <- names(cov_depth)[1]
            chr <- names(cov_depth[[1]])
            
            for (seqname in chr){
                cov_files[[seqname]][[tag]] <- cov_depth[[1]][[seqname]]
            }
        }
        
        ## rearrange sample-wise coverage Rle objects into chromosome-wise coverage
        ## Rle objects
        for (seqname in names(cov_files)){
            chr_cov <- list()
            for (tag in names(cov_files[[seqname]])){
                chr_cov[[tag]] <- readRDS(cov_files[[seqname]][[tag]])
            }
            filename <- file.path(od, paste(seqname, 
                                            Sys.Date(), "_RleCov.RDS"))
            
            ## chr_cov: a list with each of its element containing a Rle of a given
            ## chromosome of each sample
            saveRDS(chr_cov, file = filename)
            filenames[[seqname]] <- filename
        }
        cov_depth <- list(cov = filenames, depth = depth)
    } else {
        for (seqname in unique(res$chr)){
            chr_cov <- list()
            for (tag in unique(res$tag)){
                chr_cov[[tag]] <- readRDS(res$coverage_file[res$tag == tag & 
                                                  res$chr == seqname])
            }
            filename <- file.path(od, paste0(seqname, 
                                             Sys.Date(), "_RleCov.RDS"))
            
            ## chr_cov: a list with each of its element containing a Rle of a given
            ## chromosome of each sample
            saveRDS(chr_cov, file = filename)
            filenames[[seqname]] <- filename
        }
        cov_depth <- list(cov = filenames, depth = depth)
    }
    
    filename_df <- data.frame(chr = names(filenames), 
                              coverage_file = unlist(filenames))
    dbBegin(db_conn)
    res <- dbSendQuery(db_conn, 
                       "INSERT INTO 
                           chromosome_coverage (chr, coverage_file) 
                           VALUES (:chr, :coverage_file);", 
                       filename_df)
    dbClearResult(res)
    dbCommit(db_conn)
    dbDisconnect(db_conn)
    
    ## return a list of filenames pointing to a list of chromosome-wise
    ## coverage of each sample and total read depth of each sample
    cov_depth
}