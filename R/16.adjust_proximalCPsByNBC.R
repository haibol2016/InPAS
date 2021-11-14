#' adjust the proximal CP sites by using Naive Bayes classifier from cleanUpdTSeq
#'
#' adjust the proximal CP sites by using Naive Bayes classifier from cleanUpdTSeq
#'
#' @param idx.list the offset of positions of CP sites
#' @param cov_diff.list  the MSE values
#' @param seqnames a character(n) vector, the chromosome/scaffolds' names
#' @param starts starts
#' @param strands strands
#' @param genome a [BSgenome::BSgenome-class] object
#' @param classifier cleanUpdTSeq classifier
#' @param classifier_cutoff  cutoff value of the classifier
#' @param shift_range the searching range for the better CP sites
#' @param search_point_START just in case there is no better CP sites
#' @param step adjusting step, default 1, means adjust by each base by
#'   cleanUpdTSeq.
#'
#' @return the offset of positions of CP sites after filter
#' @details the step for calculating is 10, can not do every base base it is
#'   really very slow.
#' @keywords internal
#' @author Jianhong Ou
#' @seealso [adjust_proximalCPsByPWM()], [get_PAscore2()]
#'

adjust_proximalCPsByNBC <- function(idx.list, 
                                      cov_diff.list,
                                      seqnames, 
                                      starts, 
                                      strands,
                                      genome, 
                                      classifier, 
                                      classifier_cutoff,
                                      shift_range, 
                                      search_point_START, 
                                      step = 1) {
  idx.len <- sapply(idx.list, length)
  offsite <- 10^nchar(as.character(max(idx.len) * ceiling(shift_range / step)))
  pos.matrix <- mapply(function(idx, start, strand, cov_diff, ID) {
    if (length(idx) == 0) {
      return(NULL)
    }
    if (is.na(idx[1])) {
      return(NULL)
    }
    idx.gp <- 1:length(idx)
    if (shift_range > step) {
      idx_lo <- idx - shift_range
      idx_up <- idx + shift_range
      idx <- mapply(function(a, b, c) {
        unique(sort(c(seq(a, b, by = step), c)))
      },
      idx_lo, idx_up, idx,
      SIMPLIFY = FALSE
      )
      idx.gp <- rep(1:length(idx_lo), sapply(idx, length))
      idx <- as.integer(unlist(idx))
      idx <- cbind(idx, idx.gp)
      idx <- idx[idx[, 1] >= search_point_START &
        idx[, 1] < length(cov_diff), , drop = FALSE]
      idx <- idx[!duplicated(idx[, 1]), , drop = FALSE]
      idx.gp <- idx[, 2]
      idx <- idx[, 1]
    }
    if (length(idx) > 0) {
      pos <- if (strand == "+") start + idx - 1 else start - idx + 1
      cbind(pos, idx, idx.gp, ID)
    } else {
      NULL
    }
  }, idx.list, starts, strands, cov_diff.list,
  1:length(idx.list),
  SIMPLIFY = FALSE
  )
  pos.matrix <- do.call(rbind, pos.matrix)

  if (length(pos.matrix) > 0) {
    idx <- get_PAscore2(
      seqnames[pos.matrix[, "ID"]],
      pos.matrix[, "pos"],
      strands[pos.matrix[, "ID"]],
      pos.matrix[, "idx"],
      pos.matrix[, "ID"] * offsite + pos.matrix[, "idx.gp"],
      genome, classifier, classifier_cutoff
    )
    idx$ID <- floor(idx$idx.gp / offsite)
    idx <- idx[!duplicated(idx$ID), ]
    IDs <- unique(pos.matrix[, "ID"])
    idx.idx <- idx$idx[match(IDs, idx$ID)]
    idx.idx[is.na(idx.idx)] <- sapply(idx.list[IDs[is.na(idx.idx)]], `[`, 1)
    idx.list[IDs] <- idx.idx
  }
  idx.list
}
