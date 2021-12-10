#' remove the converging candidates 3' UTRs LIKE UTR3___UTR3
#'
#' some of the results is from connected two 3' UTRs. We want to remove them.
#'
#' The algorithm need to be improved.
#'
#' @param x the distal 3UTR coverage
#'
#' @return the 3UTR coverage after removing the next 3UTR
#' @keywords internal
#' @author Jianhong Ou

remove_convergentUTR3s <- function(x) {
  ## smooth by window_size=10
  ws <- 10
  len <- length(x)
  if (len > 100) {
    y <- rowsum(x,
      group = rep(1:ceiling(len / ws), each = ws)[1:len],
      reorder = FALSE
    )
    y <- y[-length(y)]
    id <- InPAS:::find_valley(y, 1, length(y), 1, filterByPval = FALSE)
    if (length(id) > 0) {
      if (id < 3) id <- 3
      pval <- try(t.test(y[1:(id - 2)], y[-(1:(id + 2))])$p.value,
        silent = TRUE
      )
      if (is.numeric(pval) && length(pval) == 1 && !is.na(pval)) {
        if (pval < 0.001) {
          x <- x[1:((id - 1) * ws)]
        }
      }
    }
  }
  x
}
