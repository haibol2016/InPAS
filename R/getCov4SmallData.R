#' Get coverage for a small-scale experiment
#'
#' @param bedgraphs a vector of character, containing paths to bedGraph files,
#'   absolute paths or path relative to the current working directory 
#' @param tags unique sample name tags matching each bedGraph file
#' @param genome an object [BSgenome::BSgenome-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be removed
#'   from the genome
#' @param hugeData a vector of logical(1), indicating whether the data is huge.
#'   For experiments with more than 3 samples, hugeData = TRUE should be 
#'   considered.
#' @param outdir a directory with write permission for storing intermediate
#'   coverage data. This is required for hughe data, i.e., when hugeData is set
#'   to TRUE
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam()
#'
#' @return a list of coverage files or objects of [S4Vectors::Rle-class] 
#' @export

getCov4SmallData <- function(sqlite_db,
                             genome, 
                             removeScaffolds = FALSE,
                             BPPARAM = NULL){
    if (missing(sqlite_db) || length(sqlite_db) != 1) {
        stop("A path to a sqlite database is required!")
    }
    if (missing(genome)) {
        stop("genome is required.")
    }
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }
    
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
    metadata <- dbGetQuery(db_conn, "SELECT * FROM metadata;")
    bedgraphs <- metadata$bedgraph_file
    tags <- metadata$tag
    
    if (nrow(metadata) > 6){
        stop("For an experiment with more than 6 samples, please use the ",
             "functions bedgraph2RleCov() and getCov4HugeData subsequentially")
    }

    if (is.null(BPPARAM)){
        cov <- mapply(function(.x, .y){
            bedgraph2RleCov(bedgraph = .x,
                            tag = .y,
                            genome = genome, 
                            removeScaffolds = removeScaffolds,
                            hugeData = FALSE)
        }, bedgraphs, tags, SIMPLIFY = FALSE)
    } else {
        cov <- bptry(bpmapply(bedgraph2RleCov, 
                              bedgraphs, tags, 
                              MoreArgs = list(genome = genome, 
                                              removeScaffolds = removeScaffolds,
                                              hugeData = FALSE), 
                              SIMPLIFY = FALSE, 
                              USE.NAMES = FALSE, 
                              BPPARAM = BPPARAM))
        while (!all(bpok(cov))) {
            cov <- bptry(bpmapply(bedgraph2RleCov, 
                                  bedgraphs, tags, 
                                  MoreArgs = list(genome = genome, 
                                                  removeScaffolds = removeScaffolds,
                                                  hugeData = FALSE), 
                                  SIMPLIFY = FALSE, 
                                  USE.NAMES = FALSE,
                                  PREDO = cov,
                                  BPPARAM = BPPARAM))
        }
    }

    # cov[[tag]][[tag]] to get a list of Rle coverage per chromosome
    # cov[[tag]][[depth]] to get depth
    names(cov) <- tags
    
    ## convert sample-oriented coverage to chromosome-oriented coverage
    cov_depth <- list(cov = list(), 
                      depth = vector("integer", length(bedgraphs)))
    
    for (i in seq_along(tags)){
        for (seqname in names(cov[[1]][[1]])) {
            cov_depth$cov[[seqname]][[tags[i]]] <- cov[[tags[i]]][[tags[i]]][[seqname]]
        }
        cov_depth$depth[i] <- cov[[tags[i]]]$depth
    }
    names(cov_depth$depth) <- tags
    
    for (d in names(cov_depth$depth)){
        dbSendQuery(db_conn,
                    paste("UPDATE metadata SET depth = ",
                          cov_depth$depth[d], "WHERE tag = ",
                          d, ";"))
    }
    
    
    ## return a list with the first element containing a list of chromosome-wise
    ## Rle coverage of each sample and the first element containing total read
    ## depth of each sample
    cov_depth
}