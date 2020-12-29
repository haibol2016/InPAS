#' sort GRanges
#'
#' sort GRanges on the positive/negative strand by start/end position in
#' ascending/descending order, respectively
#'
#' @param .ele an object of [GenomicRanges::GRanges-class]
#'
#' @return a sorted [GenomicRanges::GRanges-class] object
#' @import GenomicRanges
#' @keywords internal
#'

sortGR <- function(.ele) {
  str <- as.character(strand(.ele))[1] == "+"
  if (str) {
    .ele <- .ele[order(start(.ele))]
  } else {
    .ele <- .ele[order(-end(.ele))]
  }
  .ele
}
