#' remove the converging candidates 3' UTRs LIKE UTR3___UTR3
#'
#' some of the results is from connected two 3' UTRs. We want to remove them.
#'
#' The algorithm need to be improved.
#'
#' @param x the collapsed next.exon.gap coverage
#'
#' @return the collapsed next.exon.gap after removing the next 3UTR
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu

remove_convergentUTR3s <- function(x) {
  if (length(x) > 100) {
    id <- InPAS:::find_valleyBySpline(x,
      ss = 1,
      se = length(x), 
      n = 1,
      filter.last = FALSE
    )
    if (!is.na(id)) {
      x <- x[seq_len(id)]
    }
  }
  x
}
