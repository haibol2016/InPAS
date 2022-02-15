#' Calculate the total coverage
#'
#' For hugeData, coverage of samples in each condition is merged chromosome
#' by chromosome. For non-hugeData, per-chromosome coverage of all samples
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param chr.cov A list of Rle objects storing coverage per sample for a given
#'   chromosome/scaffold
#' @param seqname A character(1), the chromosome/scaffold name
#' @param metadata A data frame containing the metadata for a RNA-seq experiment,
#'  which can be extract from the SQLite database set up by [setup_sqlitedb()]
#' @param outdir A character(1) vector, a path with write permission for storing
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param hugeData A logical(1), indicating whether it is huge data
#' @return A list containing pooled coverage data. For hugeData, coverage of
#'   samples under each condition is merged chromosome by chromosome.
#'   For non-hugeData, per-chromosome coverage of all samples are returned.
#'   \describe{
#'           \item{seqname}{chromosome/scaffold name}
#'             \describe{
#'                      \item{condition1}{condition name 1}
#'                      \item{condition1}{condition name 2}
#'                      }
#'   }
#' @keywords internal
#' @import S4Vectors RSQLite
#' @author Haibo Liu, Jianhong Ou

get_totalCov <- function(sqlite_db,
                         chr.cov,
                         seqname,
                         metadata,
                         outdir,
                         hugeData) {
  lock_filename <- getLockName()
  if (!file.exists(lock_filename)) {
    stop(
      "lock_filename must be an existing file.",
      "Please call addLockName() first!"
    )
  }
  if (hugeData) {
    cov <- list()
    condition_tag <- split(metadata$tag, metadata$condition)
    for (condition in names(condition_tag)) {
      cov[[seqname]][[condition]] <- Rle(values = 0, lengths = 1)
      # for (tag in condition_tag[[condition]]){
      #   cov[[seqname]][[condition]] <- cov[[seqname]][[condition]] +
      #     seqname_cov[[seqname]][[tag]]
      # }
      # Reduce() is faster
      cov[[seqname]][[condition]] <-
        Reduce(function(.x, .y) {
          .x + .y
        },
        chr.cov[[seqname]][condition_tag[[condition]]],
        init = cov[[seqname]][[condition]]
        )
    }
  } else {
    cov <- chr.cov
  }
  filename_total <- file.path(
    outdir,
    paste(seqname,
      "totalCov.RDS",
      sep = "_"
    )
  )

  ## chr_cov: chr_cov[[seqname]][[tag]]
  saveRDS(cov, file = filename_total)

  tryCatch(
    {
      file_lock <- flock::lock(lock_filename)
      db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
      res <- dbSendStatement(
        db_conn,
        paste0("DELETE FROM total_coverage
                                  WHERE chr = '", seqname, "';")
      )
      dbClearResult(res)

      ## insert into database the per-chromosome utr3_total_coverage
      res <- dbSendStatement(
        db_conn,
        paste0(
          "INSERT INTO
                           total_coverage (chr, coverage_file)
                           VALUES ('", seqname, "',", "'",
          filename_total, "');"
        )
      )
      dbClearResult(res)
    },
    error = function(e) {
      print(paste(conditionMessage(e)))
    },
    finally = {
      dbDisconnect(db_conn)
      unlock(file_lock)
    }
  )
  cov
}
