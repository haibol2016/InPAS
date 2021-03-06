#' adjust the proximal CP sites
#'
#' adjust the proximal CP sites by PolyA PWM and cleanUpdTSeq
#'
#' @param CPs the outputs of [search_proximalCPs()]
#' @param MINSIZE min size for short from
#' @param PolyA_PWM PolyA position weight matrix
#' @param genome a [BSgenome::BSgenome-class] object
#' @param classifier cleanUpdTSeq classifier
#' @param classifier_cutoff cutoff value of the classifier
#' @param shift_range the searching range for the better CP sites
#' @param search_point_START just in case there is no better CP sites
#' @param step adjust step, default 1, means adjust by each base by
#'   cleanUpdTSeq.
#'
#' @return keep same as [search_proximalCPs()], which can be handled by
#'   [polish_CPs()].
#' @seealso [search_proximalCPs()], [polish_CPs()], [adjust_proximalCPsByPWM()],
#'   [adjust_proximalCPsByNBC()], [get_PAscore()],[get_PAscore2()]
#' @keywords internal
#' @author Jianhong Ou

adjust_proximalCPs <- function(CPs, 
                        MINSIZE, 
                        PolyA_PWM,
                        genome, 
                        classifier, 
                        classifier_cutoff,
                        shift_range, 
                        search_point_START, 
                        step = 1) {
  dCPs <- CPs$dCPs
  flag <- CPs$flag
  seqnames <- as.character(dCPs$seqnames)
  strands <- as.character(dCPs$strand)
  starts <- coors <-
    lapply(CPs$chr.cov.merge, function(.ele) as.numeric(rownames(.ele)))
  starts[strands == "-"] <- lapply(starts[strands == "-"], rev)
  starts <- sapply(starts, `[`, 1)
  idx.list <- CPs$Predicted_Proximal_APA
  if (is(PolyA_PWM, "matrix")) {
    idx.list <- adjust_proximalCPsByPWM(
      idx.list, PolyA_PWM, seqnames, starts,
      strands, genome, shift_range,
      search_point_START
    )
  }
  cov_diff.list <- CPs$fit_value
  if (is(classifier, "PASclassifier")) {
    idx.list <- adjust_proximalCPsByNBC(
      idx.list, cov_diff.list,
      seqnames, starts, strands,
      genome,
      classifier, classifier_cutoff,
      shift_range, search_point_START,
      step
    )
  }
  CPs$Predicted_Proximal_APA[flag] <- idx.list[flag]

  CPs
}
