#' prepare coverage data and fitting data for plot
#'
#'
#' @param gr an object of [GenomicRanges::GRanges-class]
#' @param coverage an object of [S4Vectors::Rle-class]
#' @param proximalSites proximal CP sites
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param groupList a list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#'
#' @return an object of [GenomicRanges::GRanges-class] with metadata:
#'   \item{dat}{matrix, first column is the fit data, the other columns are
#'   coverage data for each sample}
#'   \item{offset}{offset from the start of 3UTR}
#' @import GenomicRanges
#' @export
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' path <- system.file("extdata", package = "InPAS")
#' bedgraphs <- c(
#'   file.path(path, "Baf3.extract.bedgraph"),
#'   file.path(path, "UM15.extract.bedgraph")
#' )
#' coverage <- coverageFromBedGraph(
#'   bedgraphs,
#'   tags = c("Baf3", "UM15"),
#'   genome = Mmusculus,
#'   removeScaffolds = TRUE,
#'   hugeData = FALSE
#' )
#' gr <- GRanges("chr6", IRanges(128846245, 128850081),
#'   strand = "-"
#' )
#' dat <- usage4plot(gr, coverage,
#'   proximalSites = 128849148, Mmusculus
#' )
#' data <- dat$dat[[1]]
#' op <- par(mfrow = c(3, 1))
#' plot(data[, 1],
#'   type = "l",
#'   xlab = "",
#'   ylab = "The fitted value"
#' )
#' abline(v = dat$offset)
#' plot(data[, 2],
#'   type = "l",
#'   xlab = "",
#'   ylab = "Baf3"
#' )
#' plot(data[, 3],
#'   type = "l",
#'   xlab = "",
#'   ylab = "UM15"
#' )
#' par(op)
usage4plot <- function(gr, coverage, proximalSites, genome,
                       groupList = NULL) {
  gcCompensation <- NA
  mappabilityCompensation <- NA
  FFT <- FALSE
  fft.sm.power <- 20
  if (!is(gr, "GRanges")) stop("gr must be an object of GRanges")
  if (is.null(names(gr))) {
    names(gr) <-
      paste("X", formatC(1:length(gr),
        width = nchar(length(gr)), flag = "0"
      ), sep = "")
  }
  if (is.null(names(proximalSites))) {
    names(proximalSites) <- names(gr)
  } else {
    if (!all(names(gr) %in% names(proximalSites))) {
      stop("not all GRanges has proximalSites")
    }
  }

  ## only need to do this for the specified gr??

  hugeData <- is.character(coverage[[1]])
  depth.weight <- depthWeight(coverage, hugeData)
  totalCov <- totalCoverage(coverage, genome, hugeData, groupList)
  ## get coverage for each region
  utr3TotalCov <-
    UTR3TotalCoverage(gr, totalCov,
      gcCompensation = gcCompensation,
      mappabilityCompensation = mappabilityCompensation,
      FFT = FFT, fft.sm.power = fft.sm.power
    )
  ## get least square error
  grs <- split(gr, as.character(seqnames(gr)))
  seqnames <- names(grs)
  datInfo <- lapply(utr3TotalCov[seqnames], function(.cov) {
    .gr <- gr[names(.cov)]
    data <- mapply(function(.ele, .str) {
      if (hugeData) {
        if (.str) .ele <- rev(.ele)
        se <- length(.ele) - 1
        os <- optimalSegmentation(.ele,
          search_point_START = 1,
          search_point_END = se
        )
        cov_diff <- os$cov_diff
        cvg <- .ele[1:se]
      } else {
        if (.str) .ele <- .ele[nrow(.ele):1, , drop = FALSE]
        se <- nrow(.ele) - 1
        os <- apply(.ele, 2, optimalSegmentation,
          search_point_START = 1, search_point_END = se
        )
        cov_diff <- sapply(os, "[[", "cov_diff")
        cov_diff <- rowMeans(cov_diff)
        cvg <- .ele[1:se, , drop = FALSE]
      }
      cbind(cov_diff, cvg)
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
