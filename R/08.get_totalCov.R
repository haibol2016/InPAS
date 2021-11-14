#' Calculate the total coverage
#'
#' For hugeData, coverage of samples in each condition is merged chromosome
#' by chromosome. For non-hugeData, per-chromosome coverage of all samples
#' are returned.
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param seqname A character(1), the chromosome/scaffold name
#' @param hugeData A logical(1), indicating whether it is huge data
#' @return A list containing coverage data. For hugeData, coverage of samples
#'   in each condition is merged chromosome by chromosome. For non-hugeData, 
#'   per-chromosome coverage of all samples are returned.
#'   \describe{
#'           \item{seqname}{chromosome/scaffold name}
#'             \describe{
#'                      \item{tag1}{name tag for sample1}
#'                      \item{tag2}{name tag for sample2}
#'                      \item{tagN}{name tag for sampleN}
#'                      }
#'   }
#' @keywords internal
#' @import S4Vectors RSQLite
#' @author  Haibo Liu, Jianhong Ou

get_totalCov <- function(sqlite_db, seqname, hugeData) {
  if (length(seqname)!= 1 || !is.character(seqname)) {
    stop("seqname must be a character vector of length 1")
  }
  
  db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
  metadata <- dbReadTable(db_conn, "metadata")
  chromsome_coverage <- dbGetQuery(db_conn, 
                      paste0("SELECT * FROM chromosome_coverage WHERE chr = '",
                            seqname, "';"))
  dbDisconnect(db_conn)
  
  if (nrow(chromsome_coverage) == 0){
    stop("The seqname is invalid, no coverage data for this seqname!")
  }
  seqname_cov <- readRDS(chromsome_coverage$coverage_file)
  
  # num_conditions <- length(unique(metadata$condition))
  if (hugeData) {
    cov <- list()
    condition_tag <- split(metadata$tag, metadata$condition)
    for (condition in names(condition_tag)){
      cov[[seqname]][[condition]] <- Rle(values = 0, lengths = 1)
      for (tag in condition_tag[[condition]]){
        cov[[seqname]][[condition]] <- cov[[seqname]][[condition]] + 
          seqname_cov[[seqname]][[tag]]
      }
    }
  } else {
    cov <- seqname_cov
  }
  
  ## cov[[seqname]][[tag]]
  cov
}
