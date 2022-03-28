#' Get 3' UTRs and their last CDS regions based on CP sites
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   setup_sqlitedb().
#' @param chr.utr3 An object of [GenomicRanges::GRanges-class], specifying UTR3
#'   GRanges for a chromosome. It must be one element of an output of
#'   [extract_UTR3Anno()].
#' @param outdir A character(1) vector, a path with write permission for storing
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param min.length.diff An integer(1) vector, specifying minimal length 
#' difference between proximal and distal APA sites which should be met to be 
#' considered for differential APA analysis. Default is 200 bp.
#'
#' @return An object of [GenomicRanges::GRanges-class] containing GRanges for
#'   UTRs with alternative CP sites and the corresponding last CDSs.
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu

get_UTR3CDS <- function(sqlite_db,
                        chr.utr3,
                        outdir = getInPASOutputDirectory(),
                        min.length.diff = 200) {
  if (missing(sqlite_db) || missing(chr.utr3)) {
    stop("sqlite_db and chr.utr3 are required")
  }
  if (length(sqlite_db) != 1 || !file.exists(sqlite_db)) {
    stop("The path to the SQLite database is invalid!")
  }
  if (!is(chr.utr3, "GRanges") ||
    !all(chr.utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
    stop(
      "chr.utr3 must be one element of an output of function of",
      "extract_UTR3Anno()"
    )
  }
  if (!is.character(outdir) || length(outdir) != 1) {
    stop("An explicit output directory is required")
  } else {
    outdir <- file.path(outdir, "007.UTR3CDS.out")
    if (!dir.exists(outdir)) {
      dir.create(outdir,
        recursive = TRUE,
        showWarnings = TRUE
      )
    }
    outdir <- normalizePath(outdir)
  }
  lock_filename <- getLockName()
  if (is.null(lock_filename) || !file.exists(lock_filename)) {
    stop(
      "lock_filename must be an existing file.",
      "Please call addLockName() first!"
    )
  }

  seqname <- unique(as.character(seqnames(chr.utr3)))
  if (length(seqname) != 1) {
    stop("chr.utr3 should contains only one chromosome/scaffold")
  }
  file_lock <- flock::lock(lock_filename)
  db_conn <- dbConnect(
    drv = RSQLite::SQLite(),
    dbname = sqlite_db
  )
  CPsites <- dbReadTable(db_conn, "CPsites")
  dbDisconnect(db_conn)
  unlock(file_lock)
  chr.CPsites <- CPsites$cpsites_file[CPsites$chr == seqname]
  if (length(chr.CPsites) != 1) {
    stop(paste("No CP sites for", seqname))
  }

  chr.CPsites <- readRDS(chr.CPsites)
  
  ## filter invalid APA entries
  chr.CPsites <- chr.CPsites[!is.na(chr.CPsites$Predicted_Proximal_APA) & 
                 abs(chr.CPsites$Predicted_Proximal_APA - 
                       chr.CPsites$Predicted_Distal_APA) >= min.length.diff]
  if (length(chr.CPsites) == 0) {
    return(NULL)
  }
  chr.CPsites <- split(chr.CPsites, names(chr.CPsites))
  utr3.regions <- lapply(chr.CPsites, get_UTR3region)
  utr3.regions <- unlist(GRangesList(utr3.regions))

  ## merge utr3.regions with CDS of the same set of transcripts whose
  ## utr3.regions are considered
  CDS <- chr.utr3 %>%
    plyranges::filter(feature == "CDS" &
      transcript %in% unique(utr3.regions$transcript))
  utr3.cds.regions <- c(utr3.regions, CDS)
  saveRDS(utr3.cds.regions, 
          file = file.path(outdir, paste0(seqname, "_UTR3CDS.RDS")))
  utr3.cds.regions
}
