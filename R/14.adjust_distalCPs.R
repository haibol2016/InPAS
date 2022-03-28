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
#' @param seqname A character(1) vector, specifying a chromososome/scaffold name
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
                             seqname,
                             step = 1) {
  if (is(classifier, "PASclassifier")) {
    dCPs <- distalCPs$dCPs
    final.utr3 <- distalCPs$final.utr3
    truncated <- dCPs$truncated
    chr.length <- seqlengths(genome)[seqname]
    ## only adjust distal APA if there is a proximal APA
    ## to speed up
    proximal.APA <- lapply(distalCPs$Predicted_Proximal_APA, 
                          function(x) {!all(is.na(x))})
    
    generate_gapCov <- function(gap, .proximal.APA, truncated,
                                cp, ID, strand, chr.length) {
      ## longer 3'UTR must be at least 200 nt long
      if (cp > 0 && .proximal.APA && !truncated) {
        ## base coordinates of refined gap
        coor <- as.integer(gsub("^.*_SEP_", "", names(gap)))
        ## adjusting start from the end of gap
        start <- coor[length(coor)]
        ## adjust ending backward only? No, search for downstream is preferred
        end <- ifelse(strand == "+",
                      coor[length(coor)] + shift_range,
                      coor[length(coor)] - shift_range)
        pos <- seq(start, end, by = ifelse(strand == "+",
                                            1 * step, -1 * step))
        idx <- match(pos, coor) ## not needed
        # pos: relative position of search start to end
        # idx: absolute coordinates for pos
        # ID: final.utr3 unique id (from 1 to Nth)
        pos.mat <- cbind(pos, idx, ID)
        
        # remove out of range position
        pos.mat <- pos.mat[pos.mat[, "pos"] > 0 & 
                             pos.mat[, "pos"] <= chr.length, ]
        return(pos.mat)
      } else {
        return(NULL)
      }
    }
    
    gap.cov <- mapply(generate_gapCov,
                      final.utr3,
                      proximal.APA,
                      truncated,
                      dCPs$preliminary_distal_APA, 
                      seq_along(final.utr3),
                      dCPs$strand,
                      chr.length,
                      SIMPLIFY = FALSE
    )
    gap.cov <- do.call(rbind, gap.cov)
    if (length(gap.cov) > 0) {
      idx <- InPAS:::get_PAscore2(
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
        dCPs[idx$idx.gp, "NBC_adjusted_distal_APA"] <-
          as.numeric(gsub(
            ".*_(\\d+)_.+",
            "\\1", idx$peak_name
          ))
      }
      dCPs$Predicted_Distal_APA[!is.na(dCPs$NBC_adjusted_distal_APA)] <-
        dCPs$NBC_adjusted_distal_APA[!is.na(dCPs$NBC_adjusted_distal_APA)]
    }
    distalCPs$dCPs <- dCPs
  }
  distalCPs
}
