#' prepare coverage data and fitting data for plot
#'
#' @param gr An object of [GenomicRanges::GRanges-class]
#' @param proximalSites An integer(n) vector, specifying the coordinates of
#'   proximal CP sites. Each of the proximal sites must match one entry in the
#'   GRanges object, \code{gr}.
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param hugeData A logical(1), indicating whether it is huge data
#'
#' @return An object of [GenomicRanges::GRanges-class] with metadata:
#'   \item{dat}{A data.frame, first column is the position, the other columns are
#'   Coverage and value}
#'   \item{offset}{offset from the start of 3' UTR}
#' @import GenomicRanges
#' @importFrom reshape2 melt
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'
#' ## load UTR3 annotation and convert it into a GRangesList
#' data(utr3.mm10)
#' utr3 <- split(utr3.mm10, seqnames(utr3.mm10), drop = TRUE)
#'
#' bedgraphs <- system.file("extdata", c(
#'   "Baf3.extract.bedgraph",
#'   "UM15.extract.bedgraph"
#' ),
#' package = "InPAS"
#' )
#' tags <- c("Baf3", "UM15")
#' metadata <- data.frame(
#'   tag = tags,
#'   condition = c("baf", "UM15"),
#'   bedgraph_file = bedgraphs
#' )
#' outdir <- tempdir()
#' write.table(metadata,
#'   file = file.path(outdir, "metadata.txt"),
#'   sep = "\t", quote = FALSE, row.names = FALSE
#' )
#'
#' sqlite_db <- setup_sqlitedb(
#'   metadata = file.path(
#'     outdir,
#'     "metadata.txt"
#'   ),
#'   outdir
#' )
#' addLockName(filename = tempfile())
#' coverage <- list()
#' for (i in seq_along(bedgraphs)) {
#'   coverage[[tags[i]]] <- get_ssRleCov(
#'     bedgraph = bedgraphs[i],
#'     tag = tags[i],
#'     genome = genome,
#'     sqlite_db = sqlite_db,
#'     outdir = outdir,
#'     chr2exclude = "chrM"
#'   )
#' }
#' data4CPsSearch <- setup_CPsSearch(sqlite_db,
#'   genome,
#'   chr.utr3 = utr3[["chr6"]],
#'   seqname = "chr6",
#'   background = "10K",
#'   TxDb = TxDb,
#'   hugeData = TRUE,
#'   outdir = outdir
#' )
#'
#' gr <- GRanges("chr6", IRanges(128846245, 128850081), strand = "-")
#' names(gr) <- "chr6:128846245-128850081"
#' data4plot <- get_usage4plot(gr,
#'   proximalSites = 128849148,
#'   sqlite_db,
#'   hugeData = TRUE
#' )
#' plot_utr3Usage(
#'   usage_data = data4plot,
#'   vline_color = "purple",
#'   vline_type = "dashed"
#' )
get_usage4plot <- function(gr,
                           proximalSites,
                           sqlite_db,
                           hugeData) {
  gcCompensation <- NA
  mappabilityCompensation <- NA
  FFT <- FALSE
  fft.sm.power <- 20

  if (!is(gr, "GRanges")) {
    stop("gr must be an object of GRanges")
  }
  if (is.null(names(gr))) {
    names(gr) <-
      paste("X", formatC(1:length(gr),
        width = nchar(length(gr)),
        flag = "0"
      ), sep = "")
  }
  if (is.null(names(proximalSites))) {
    names(proximalSites) <- names(gr)
  } else {
    if (!all(names(gr) %in% names(proximalSites))) {
      stop("not all GRanges has proximalSites")
    }
  }

  grs <- split(gr, as.character(seqnames(gr)))
  seqnames <- names(grs)

  db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
  metadata <- dbReadTable(db_conn, "metadata")
  total_cov <- dbReadTable(db_conn, "total_coverage")
  on.exit(dbDisconnect(db_conn))
  ## only need to do this for the specified gr??
  depth.weight <- get_depthWeight(
    metadata = metadata,
    hugeData
  )

  total_cov <- total_cov[total_cov$chr %in% seqnames, ]
  if (nrow(total_cov) < 1) {
    stop("chromosomes in gr are not covered")
  }
  totalCov <- lapply(total_cov$coverage_file, function(.x) {
    readRDS(.x)[[1]]
  })
  names(totalCov) <- total_cov$chr

  ## get coverage for each region
  utr3TotalCov <- list()
  for (seqname in seqnames) {
    utr3TotalCov[[seqname]] <-
      get_UTR3TotalCov(grs[[seqname]], totalCov[seqname],
        gcCompensation = gcCompensation,
        mappabilityCompensation =
          mappabilityCompensation,
        FFT = FFT,
        fft.sm.power = fft.sm.power
      )
  }

  ## get least square error
  datInfo <- lapply(utr3TotalCov, function(.cov) {
    .gr <- gr[names(.cov)]
    data <- mapply(function(.ele, .str) {
      # reverse coverage data if strand is "-"
      if (.str) .ele <- .ele[nrow(.ele):1, , drop = FALSE]
      se <- nrow(.ele) - 1
      ## normalized coverage
      .ele <- .ele / depth.weight
      cov_diff <- apply( .ele, 2, calculate_mse,
        search_point_START = 1,
        search_point_END = se
      )
      cov_diff <- rowMeans(cov_diff)
      cvg <- .ele[1:se, , drop = FALSE]

      # reformat data from wide to long
      dat <- as.data.frame(cbind(
        Position = 1:nrow(cvg),
        MSE = cov_diff, cvg
      ))
      dat <- reshape2::melt(dat,
        variable.name = "Coverage",
        measure.vars = colnames(dat)[2:4]
      )
      dat$Coverage <- as.character(dat$Coverage)
      dat$Coverage <- ifelse(dat$Coverage != "MSE",
        paste(dat$Coverage, "coverage"),
        dat$Coverage
      )
      dat$Coverage <- factor(dat$Coverage,
        levels = c(
          "MSE",
          unique(dat$Coverage[dat$Coverage != "MSE"])
        )
      )
      dat
    }, .cov, as.character(strand(.gr)) == "-", SIMPLIFY = FALSE)
    .gr$dat <- data
    .gr$offset <- ifelse(strand(.gr) == "+",
      proximalSites[names(.cov)] - start(.gr),
      end(.gr) - proximalSites[names(.cov)] + 1
    )
    .gr
  })
  datInfo <- unlist(GRangesList(datInfo))
}
