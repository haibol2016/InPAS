#' Get Rle coverage from a bedgraph file for a sample
#'
#' Get RLe coverage from a bedgraph file for a sample
#'
#' @param bedgraph A path to a bedGraph file
#' @param tag A character(1) vector, a name tag used to label the bedgraph file.
#'   It must match the tag specified in the metadata file used to setup the
#'   SQLite database
#' @param genome an object [BSgenome::BSgenome-class]. To make things easy, we 
#'   suggest users creating a [BSgenome::BSgenome-class] instance from the 
#'   reference genome used for read alignment. For details, see the 
#'   documentation of [BSgenome::forgeBSgenomeDataPkg()].
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param future.chunk.size The average number of elements per future 
#'   ("chunk"). If Inf, then all elements are processed in a single future.
#'   If NULL, then argument future.scheduling = 1 is used by default. Users can
#'   set future.chunk.size = total number of elements/number of cores set for 
#'   the backend. See the future.apply package for details. You may adjust 
#'   this number based based on the available computing resource: CPUs and RAM. 
#'   This parameter affects the time for converting coverage from bedgraph to Rle.
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or 
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @return A data frame, as described below.
#'   \describe{
#'            \item{tag}{the sample tag}
#'            \item{chr}{chromosome name}
#'            \item{coverage_file}{path to Rle coverage files for each chromosome
#'             per sample tag}
#'           }
#' @importFrom dplyr as_tibble mutate filter arrange bind_rows group_by
#'   left_join summarise n bind_cols syms desc pull
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle mapSeqlevels
#' @import readr RSQLite S4Vectors GenomicRanges flock
#' @importFrom future.apply future_lapply
#' @importFrom magrittr %>%
#' 
#' @export
#' @author Jianhong Ou, Haibo Liu
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
#'    addLockName()
#'    coverage_info <- get_ssRleCov(bedgraph = bedgraphs[1], 
#'                           tag = tags[1],
#'                           genome = genome,
#'                           sqlite_db = sqlite_db,
#'                           outdir = outdir,
#'                           chr2exclude = "chrM")
#'    # check read coverage depth
#'    db_connect <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
#'    dbReadTable(db_connect, "metadata")
#'    dbDisconnect(db_connect)
#' }

get_ssRleCov <- function(bedgraph,
                         tag,
                         genome = getInPASGenome(),
                         sqlite_db,
                         future.chunk.size = NULL,
                         outdir = getInPASOutputDirectory(),
                         chr2exclude = getChr2Exclude()) {
    if (!is(genome, "BSgenome")) {
        stop("genome,an object of BSgenome, required.")
    }
    if (missing(tag) || missing(bedgraph)) {
        stop("tags and bedgraphs are required.")
    }
    if (length(bedgraph) != 1 || length(tag) != 1) {
        stop("length of bedgraph must be 1")
    }
    if (!file.exists(bedgraph)) {
        stop("bedgraph file does not exist")
    }
    if (!is.character(outdir) || length(outdir) != 1)
    {
        stop("A writable output directory is required!")
    }
    if (missing(sqlite_db) || !file.exists(sqlite_db)){
        stop("The SQLite database for InPAS is required!")
    }
    if (!is.null(chr2exclude) && !is.character(chr2exclude))
    {
        stop("chr2Exclude must be NULL, or a character vector")
    }
    lock_filename <- getLockName()
    if (!file.exists(lock_filename)) 
    {
        stop("lock_filename must be an existing file.",
             "Please call addLockName() first!")
    }
    seqnames.bedfile <-
        as.data.frame(read_tsv(
            bedgraph, comment = "#",
            col_names = FALSE, skip = 0,
            col_types = cols(X1 = col_factor(), 
                             X2 = "-", X3 = "-", X4 = "-")
        ))[, 1]
    
    seqnames <- trim_seqnames(genome, chr2exclude)
    
    seqStyle <- seqlevelsStyle(genome)
    seqStyle.bed <- seqlevelsStyle(levels(seqnames.bedfile))
    
    if (any(seqStyle.bed != seqStyle)) {
        stop("seqlevelsStyle of genome is different from bedgraph file.")
     }
    seqnames <- sort(intersect(levels(seqnames.bedfile), seqnames))
    if (length(seqnames) < 1) {
        stop(paste("there is no intersect seqname in", 
                   bedgraph, "and the genome"))
    }
    
    ## save coverage file of each chr of a sample into a directory
    ## called "samplewise_RleCov"
    outdir <- file.path(outdir, "002.samplewise.RleCov", tag)
    if (!dir.exists(outdir))
    {
        dir.create(outdir, recursive = TRUE)
    }
    outdir <- normalizePath(outdir, mustWork = TRUE)
    
    seqLen <- get_seqLen(genome, chr2exclude)
    
    summaryFunction <- function(seqname) {
        seqL <- seqLen[seqname]
        
        ## apply some trick here by using "FALSE"
        lines2read <- Rle(c(FALSE, seqnames.bedfile == seqname))
        true <- which(runValue(lines2read))
        skip <- runLength(lines2read)[true - 1]
        skip[1] <- skip[1] - 1
        nrow <- runLength(lines2read)[true]
        
        dat <- read_tsv(bedgraph,
                        comment = "#",
                        col_names = FALSE,
                        skip = skip[1],
                        n_max = nrow[1],
                        col_types = cols(X1 = "-", X2 = "i",
                                         X3 = "i", X4 = "d"))
        
        if (length(true) > 1) { ## this should never happen for bedgraph files
            cul_skip <- skip[1] + nrow[1]
            for (i in 2:length(true)) {
                cul_skip <- cul_skip + skip[i]
                lines <-
                    read_tsv(bedgraph,
                             comment = "#",
                             col_names = FALSE,
                             skip = cul_skip,
                             n_max = nrow[i],
                             col_types = cols(
                                 X1 = "-", X2 = "i",
                                 X3 = "i", X4 = "d"))
                dat <- dplyr::bind_rows(dat, lines)
                cul_skip <- cul_skip + nrow[i]
            }
        }
        ## convert bedgraph 0-based index to 1-based index for GRanges
        dat <- dat %>% dplyr::mutate(X2 = X2 + 1)
        
        ## convert uncovered regions as gaps
        ## This step can be avoided when generate bedgraph with bedtools:
        ## bedtools genomecov -bga -split -ibam  <coordinate-sorted BAM file>
        gaps <-
            as_tibble(as.data.frame(gaps(IRanges(
                dplyr::pull(dat, 1),
                dplyr::pull(dat, 2)),
            start = 1, end = seqL)))
        if (nrow(gaps) > 0) {
            gaps <- dplyr::bind_cols(gaps[, 1:2],
                                     V4 = rep(0, nrow(gaps)))
            colnames(gaps) <- colnames(dat)
            dat <- dplyr::bind_rows(dat, gaps)
        }
        rm(gaps)
        gc()
        
        dat <- dat %>%
            dplyr::arrange(X2) %>%
            dplyr::mutate(width = X3 - X2 + 1)
        seqname_cov <- Rle(dplyr::pull(dat, 3), dplyr::pull(dat, width))
        rm(dat)
        gc()
        
        filename <- file.path(outdir, 
                              paste(tag, seqname, "RleCov.RDS",sep = "_"))
        saveRDS(seqname_cov, file = filename)
        cov_file <- list(seqname = seqname_cov, 
                         filename = filename,
                         depth = sum(as.double(runLength(seqname_cov) *
                                               runValue(seqname_cov))))
        names(cov_file)[1] <- seqname
        cov_file
    }
    
    ## get coverage
    cvg <- future_lapply(seqnames, summaryFunction, 
                         future.stdout = NA,
                         future.chunk.size = future.chunk.size)
    depth <- sum(sapply(cvg, function(.cov){
             .cov[[3]]
         }))
    coverage_file <- sapply(cvg, function(.cov){
        .cov[[2]]
    })

    rm(cvg)
    gc()
    
    filename_df <- data.frame(tag = tag, 
                              chr = seqnames, 
                              coverage_file = coverage_file)
    tryCatch(
        {
            file_lock <- flock::lock(lock_filename)
            db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
            res <- dbSendStatement(db_conn,
                                   paste0("UPDATE metadata SET depth = ",
                                          depth, " WHERE tag = '",
                                          tag, "';"))
            dbClearResult(res) 
            ## remove existing record for the same "tag"
            res <- dbSendStatement(db_conn,
                                   paste0("DELETE FROM sample_coverage WHERE tag = '",
                                          tag, "';"))
            dbClearResult(res)
            res <- dbSendStatement(db_conn, 
                                   "INSERT INTO 
                       sample_coverage (tag, chr, coverage_file) 
                       VALUES (:tag, :chr, :coverage_file);", 
                                   filename_df)
            dbClearResult(res)
        },
        error = function(e) {
            print(paste(conditionMessage(e)))
        }, finally = {
            dbDisconnect(db_conn)
            unlock(file_lock)})
    filename_df
}
