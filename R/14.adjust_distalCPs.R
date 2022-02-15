#' Adjust distal CP sites by the cleanUpdTSeq algorithm
#'
#' Adjust distal CP sites by the cleanUpdTSeq algorithm
#'
#' @param distalCPs the output of [search_distalCPs()]
#' @param classifier An R object for Naive Bayes classifier model, like the one
#'   in the cleanUpdTSeq package.
#' @param classifier_cutoff A numeric(1) vector. A cutoff of probability that a
#'   site is classified as true CP sites. The value should be between 0.5 and 1.
#'   Default, 0.8.
#' @param shift_range An integer(1) vector, specifying a shift range for
#'   adjusting the proximal and distal CP sites. Default, 50. It determines the
#'   range flanking the candidate CP sites to search the most likely CP sites.
#' @param genome a [BSgenome::BSgenome-class] object
#' @param step An integer (1) vector, specifying the step size used for adjusting
#'   the proximal or distal CP sites using the Naive Bayes classifier from the
#'   cleanUpdTSeq package. Default 1. It can be in the range of 1 to 5.
#' @seealso [search_proximalCPs()], [get_PAscore2()]
#' @keywords internal
#' @author Jianhong Ou

adjust_distalCPs <- function(distalCPs,
                             classifier,
                             classifier_cutoff,
                             shift_range,
                             genome,
                             step = 1) {
  if (is(classifier, "PASclassifier")) {
    dCPs <- distalCPs$dCPs
    final.utr3 <- distalCPs$final.utr3
    generate_gapCov <- function(gap, cp, ID, strand) {
      if (cp > 0) {
        ## base coordinates of refined gap
        coor <- as.integer(gsub("^.*_SEP_", "", names(gap)))
        ## adjusting start from the end of gap
        start <- coor[length(coor)]
        ## adjust ending backward only?
        end <- ifelse(length(coor) > 2 * shift_range,
                      coor[length(coor) - 2 * shift_range],
                      coor[1]
        )
        pos <- seq(start, end, by = ifelse(strand == "+",
                                           -1 * step, 1 * step
        ))
        idx <- match(pos, coor) ## not needed
        
        # pos: relative position of search start to end
        # idx: absolute coordinates for pos
        # ID: final.utr3 unique id (from 1 to Nth)
        return(cbind(pos, idx, ID))
      } else {
        return(NULL)
      }
    }
    
    gap.cov <- mapply(generate_gapCov,
                      final.utr3,
                      dCPs$distalCP, 
                      seq_along(final.utr3),
                      dCPs$strand,
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
        genome, classifier,
        classifier_cutoff
      )
      ## add adjusted distal CP sites to dCPs data.frame
      if (!is.null(idx)) {
        dCPs[idx$idx.gp, "adjustedDistalCP"] <-
          as.numeric(gsub(
            ".*_(\\d+)_.+",
            "\\1", idx$peak_name
          ))
      }
      dCPs$Predicted_Distal_APA[!is.na(dCPs$adjustedDistalCP)] <-
        dCPs$adjustedDistalCP[!is.na(dCPs$adjustedDistalCP)]
    }
    distalCPs$dCPs <- dCPs
  }
  distalCPs
}
