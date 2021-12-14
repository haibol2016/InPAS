#' find the local minimal mean squared error (MSE)
#'
#' for a giving numeric vectors, calculate the top N local minimal mean squared
#' error. The sites of local minimal MSE are saved if they are in the range of
#' searching start site to searching end site
#'
#' @param x  a numeric(n) vector containing MSEs for a given range.
#' @param ss position to start searching positions of local minimal MSE
#' @param se position to end searching positions of local minimal MSE
#' @param n  the length of output. If n=-1, output all the local minimal MSE
#'   positions.
#' @param savedID saved positions of local minimal MSE
#' @param filterByPval A logical(1) vector, indicating whether tp filter the
#'   candidate positions by significant levels (p values) or not.
#'
#' @return a numeric vector containing a number of candidate CP sites.
#' @keywords internal
#' @importFrom stats quantile
#' @author Jianhong Ou

find_valley <- function(x, 
                        ss, 
                        se,
                        n = 1,
                        savedID = NA,
                        filterByPval = TRUE) {
  ## local minimal: like sign switching of second-order derivative,
  ## assuming there is no flat valley
  pos <- which(diff(c(sign(diff(c(x[1], x[ss:se]))), 0)) == 2) + ss - 1
  if (length(pos) > 0 && filterByPval) {
    ## peaks
    pos.pos <-
      which(diff(c(sign(diff(c(x[1], x[ss:se]))), 0)) == -2) + ss - 1
    pos.pos <- sort(c(pos, pos.pos, ss, se))
    w <- which(pos.pos %in% pos)
    space <- 10
    y <- mapply(function(a, b, p) {
      c <- unique(sort(seq(a, b)))
      c <- c[c != p & c > ss & c < se]
      d <- x[c]
    }, pos.pos[w - 1] - space,
    pos.pos[w - 1] + space,
    pos,
    SIMPLIFY = FALSE)

    mu <- sapply(y, mean, na.rm = TRUE)
    sigma <- sapply(y, sd, na.rm = TRUE)
    z <- (x[pos] - mu) / sigma
    p <- 2 * pnorm(-abs(z))
    pos <- pos[p < 0.001]
  }
  # pos <- unique(c(ss, pos, se))
  pos <- pos[!is.na(x[pos])]

  tobeadd <- ifelse(x[ss] < x[se], ss, se)
  if (length(pos) > 0) {
    w <- x[tobeadd] < min(x[pos])
    if (length(w) > 0 && !is.na(w[1]) && is.logical(w[1])) {
      if (any(w)) pos <- unique(c(pos, tobeadd))
    }
  } else {
    pos <- tobeadd
  }
  if (!is.na(savedID) && savedID > ss && savedID < se) {
    pos <- unique(c(pos, savedID))
  }

  pos <- pos[x[pos] <= quantile(x[ss:se], probs = .5, na.rm = TRUE)]
  if (n == -1 || n > length(pos)) n <- length(pos)
  if (length(pos) < 1) {
    return(pos)
  } else {
    pos <- pos[order(x[pos])][1:n]
  }
}
