#' Find major valleys after spline smoothing
#'
#' Find major valleys after spline smoothing
#'
#' @param x A vector of numeric(n), containing MSEs for a given range
#' @param ss An positive integer, search start site relative to the leftmost base
#' @param se An positive integer, search end site relative to the leftmost base
#' @param nknots An positive integer, the number of knots for smoothing using
#'   spline[stats::smooth.spline()]. By default, set to 10 knots per kb.
#' @param n An integer, specifying the number of location where MSE are local
#'   minima (candidate CP sites). If set to -1, return all candidate CP sites.
#' @DIST2END An integer, specifying a cutoff of the distance between last valley
#'   and the end of the 3' UTR (where MSE of the last base is calculated). If 
#'   the last valley is closer to the end than the specified distance, it will be
#'   not be considered because it is very likely due to RNA coverage decay at the
#'   end of mRNA. Default is 1200. User can consider a value between 1000 and 
#'   1500, depending on the library preparation procedures: RNA fragmentation and
#'   size selection.
#' @param min.dist An integer, minimal distance allowed between two adjacent
#'   candidate CP sites otherwise collapsed by selecting the one with lower MSE.
#' @param filter.last A logical(1), whether to filter out the last valley, which
#'   is likely the 3' end of the longer 3' UTR if no novel distal CP site is
#'   detected and the 3' end excluded by setting cutEnd/search_point_END is small.
#' @param plot A logical(1), whether to plot the MSE profile and the candidate
#'   valleys.
#' @importFrom graphics abline lines
#' @importFrom stats smooth.spline
#' @return A vector of integer.
#' @keywords internal

find_valleyBySpline <- function(x,
                                ss,
                                se = length(x),
                                nknots = ceiling((se - ss + 1) / 1000 * 10),
                                n = -1,
                                min.dist = 200,
                                filter.last = TRUE,
                                DIST2END = 1200,
                                plot = FALSE) {
  df <- data.frame(index = seq.int(length(x)), mse = x)

  ## remove the pos < ss, which are 0s
  df <- df[-c(1:(ss - 1)), ]
  ## smoothing MSE profile
  pos <- try({
    ss_pred <-suppressMessages(smooth.spline(df[, 2], 
                                             nknots = nknots, cv = TRUE))
    ## find local minima
    pos <- which(diff(c(sign(diff(c(ss_pred$y[1], ss_pred$y))), 0)) == 2) +
      ss - 1
    pos
  }, silent  = TRUE)
  while (is(pos, "try-error")) {
    nknots <- nknots * 2
    pos <- try({
      ss_pred <- suppressMessages(smooth.spline(df[, 2], 
                                                nknots = nknots, cv = TRUE))
      ## find local minima
      pos <- which(diff(c(sign(diff(c(
        ss_pred$y[1],
        ss_pred$y
      ))), 0)) == 2) + ss - 1
      pos
    }, silent  = TRUE)
  }

  ## add ss or se if it is smaller than mse[pos]
  tobeadd <- ifelse(x[ss] < x[se], ss, se)
  if (length(pos) > 0) {
    w <- x[tobeadd] < min(x[pos])
    if (length(w) > 0 && !is.na(w[1]) && is.logical(w[1])) {
      if (any(w)) pos <- unique(c(pos, tobeadd))
    }
  } else {
    pos <- tobeadd
  }

  if (length(pos) > 0) {
    pos <- pos[!is.na(x[pos])]
  }

  ## remove minor minima
  if (length(pos) > 0) {
    pos <- pos[max(x[-c(1:ss)]) - x[pos] > (max(x[-c(1:ss)]) - 
                                              min(x[-c(1:ss)])) * 0.1]
  }

  if (length(pos) >= 1) {
    ## get original position with local minimal unsmoothed MSE
    pos <- unlist(lapply(pos, function(.pos) {
      left <- .pos - 100
      right <- .pos + 100
      if (left < ss) {
        left <- ss
      }
      if (right > se) {
        right <- se
      }
      raw_mse <- x[left:right]
      .pos <- which(raw_mse == min(raw_mse)) + left - 1
      .pos
    }))
    
    ## merge adjacent local minima with min.dist by picking the smallest
    ## from neighboring positions
    if (length(pos) > 1) {
      interval_dist <- diff(pos)
      while (any(interval_dist <= min.dist)) {
        index_rmv <- numeric()
        for (i in seq_along(interval_dist))
        {
          if (interval_dist[i] <= min.dist) {
            index_rmv <- ifelse(x[pos[i]] > x[pos[i + 1]], i, i + 1)
          }
        }
        index_rmv <- unique(index_rmv)
        pos <- pos[-(index_rmv)]
        interval_dist <- diff(pos)
      }
    }
    
    ## filter last
    if (length(pos) > 0 && filter.last) {
      pos <- sort(pos)
      if (length(x) - pos[length(pos)] <= DIST2END) {
         if (x[pos[length(pos)]] == min(x[pos]))
         {
             return(NA)
         } else {
             pos <- pos[-length(pos)]
         }
      }
    }
    if (length(pos) > 0) {
      pos <- pos[x[pos] < 0.95*mean(x[ss:(ss+49)])]
    }
    
    if (length(pos) > 0) {
      if (plot) {
        plot(x = df[, 1], y = df[, 2], xlab = "relative position", ylab = "MSE")
        lines(ss_pred$x + ss - 1, ss_pred$y, col = "red")
        abline(v = pos, col = "purple")
      }
      if (n == -1 || n > length(pos)) {
        n <- length(pos)
      }
      pos <- pos[order(x[pos])][1:n]
    } else {
      pos <- NA
    }
  } else {
    pos <- NA
  }
  return(pos)
}
