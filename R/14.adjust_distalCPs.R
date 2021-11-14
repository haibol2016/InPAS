#' adjust distal CP sites by cleanUpdTSeq
#'
#' adjust distal CP sites by cleanUpdTSeq
#'
#' @param distalCPs the output of [search_distalCPs()]
#' @param classifier An R object for Naive Bayes classifier model, like the one
#'   in the cleanUpdTSeq package.
#' @param classifier_cutoff A numeric(1) vector. A cutoff of probability that a
#'   site is classified as true CP sites. The value should be between 0.5 and 1.
#'   Default, 0.8.
#' @param shift_range An integer(1) vector, specifying a shift range for 
#'   adjusting the proximal and distal CP sites. Default, 100. It determines the
#'   range flanking the candidate CP sites to search the most likely real 
#'   CP sites.
#' @param genome a [BSgenome::BSgenome-class] object
#' @param step An integer (1) vector, specifying the step size used for adjusting
#'   the proximal or distal CP sites using the Naive Bayes classifier from the
#'   cleanUpdTSeq package. Default 1. It can be in the range of 1 to 10.
#' @seealso [search_proximalCPs()], [get_PAscore2()]
#' @keywords internal
#' @author Jianhong Ou

adjust_distalCPs <- function(distalCPs, 
                      classifier, 
                      classifier_cutoff,
                      shift_range, 
                      genome, 
                      step = 1) {
  dCPs <- distalCPs$dCPs
  next.exon.gap <- distalCPs$next.exon.gap
  gap.cov <- mapply(function(gap, cp, ID, strand) {
    if (cp > 0) {
      coor <- as.integer(gsub("^.*_SEP_", "", names(gap[1:cp])))
      start <- coor[length(coor)]
      end <- ifelse(length(coor) > 2 * shift_range,
        coor[length(coor) - 2 * shift_range],
        coor[1]
      )
      pos <- seq(start, end, by = ifelse(strand == "+",
                                         -1 * step, 1 * step))
      idx <- match(pos, coor)
      cbind(pos, idx, ID)
    } else {
      NULL
    }
  }, next.exon.gap, dCPs$distalCP,
  1:length(next.exon.gap), dCPs$strand,
  SIMPLIFY = FALSE
  )

  gap.cov <- do.call(rbind, gap.cov)
  if (length(gap.cov) > 0) {
    idx <- get_PAscore2(
      dCPs$seqnames[gap.cov[, "ID"]],
      gap.cov[, "pos"],
      dCPs$strand[gap.cov[, "ID"]],
      gap.cov[, "idx"],
      gap.cov[, "ID"],
      genome, classifier, classifier_cutoff
    )
    distalCPs$dCPs[idx$idx.gp, "distalCP"] <- idx$idx
  }
  distalCPs
}
