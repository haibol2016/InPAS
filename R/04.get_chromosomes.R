#' Identify chromosomes/scaffolds for CP site discovery
#' 
#' Identify chromosomes/scaffolds which have both coverage and annotated
#' 3' utr3 for CP site discovery
#'
#' @param utr3 An object of [GenomicRanges::GRangesList-class]. An output of
#'   [extract_UTR3Anno()].
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#'
#' @return A vector of characters, containing names of chromosomes/scaffolds
#'   for CP site discovery
#' @export
#' @examples 
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' data(utr3.mm10)
#' utr3 <- split(utr3.mm10, seqnames(utr3.mm10), drop = TRUE)

#' bedgraphs <- system.file("extdata",c("Baf3.extract.bedgraph",
#'                                      "UM15.extract.bedgraph"), 
#'                          package = "InPAS")
#' tags <- c("Baf3", "UM15")
#' metadata <- data.frame(tag = tags, 
#'                        condition = c("Baf3", "UM15"),
#'                        bedgraph_file = bedgraphs)
#' outdir = tempdir()
#' write.table(metadata, file =file.path(outdir, "metadata.txt"), 
#'             sep = "\t", quote = FALSE, row.names = FALSE)
#' 
#' sqlite_db <- setup_sqlitedb(metadata = file.path(outdir, 
#'                                                  "metadata.txt"),
#'                             outdir)
#' addLockName(filename = tempfile())
#' coverage <- get_ssRleCov(bedgraph = bedgraphs[1], 
#'                          tag = tags[1],
#'                          genome = genome,
#'                          sqlite_db = sqlite_db,
#'                          outdir = outdir,
#'                          chr2exclude = "chrM",
#'                          BPPARAM = NULL)
#' get_chromosomes(utr3, sqlite_db)  
#'

get_chromosomes <- function(utr3, sqlite_db){
    if (missing(utr3) || !is(utr3, "GRangesList")){
        stop("utr3 must be a GRangesList object")
    }
    if (missing(sqlite_db) || !file.exists(sqlite_db)){
        stop("sqlite_db must be an exisiting path to SQLite database")
    }
    
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
    covered_chrs <- unique(dbReadTable(db_conn, "sample_coverage")$chr)
    dbDisconnect(db_conn)
    
    ## only consider chromosomes with utr3 annotation and coverage
    chromosomes <- intersect(names(utr3), covered_chrs)
    chromosomes
}