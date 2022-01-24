#' Calculate mean squared errors (MSE)
#'
#' Calculate mean squared errors (MSE) for each searched site which is assumed
#' bisection site (i.e. potential CP site).
#'
#' @param .ele A numeric vector, storing 3' UTR coverage for a give sample or 
#'   collapsed 3' UTR coverage for a given condition
#' @param search_point_START An integer, specifying the start position to 
#'   calculate MSE
#' @param search_point_END An integer, specifying end position to calculate MSE
#'
#' @return a vector of numeric, containing mean squared errors for each searched
#'   site when which is assumed as a bisection site (i.e. potential CP site).
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu

calculate_mse <- function(.ele,
                          search_point_START,
                          search_point_END) {
  .l <- length(.ele)
  short_UTR_abun <- cumsum(.ele) / (1:.l) ## averaged coverage as abundance

  ## cumulative coverage from end to start
  long_UTR_abun <- cumsum(rev(.ele)) / (1:.l)

  ## reverse the cumulative coverage from end to start so it
  ## the matches the indexing order of "short_UTR_abun"
  long_UTR_abun <- rev(long_UTR_abun)

  short_UTR_abun <- short_UTR_abun[-length(short_UTR_abun)]
  long_UTR_abun <- long_UTR_abun[-1]
  short_UTR_abun <- short_UTR_abun - long_UTR_abun
  short_UTR_abun <- ifelse(short_UTR_abun < 0, 0, short_UTR_abun)

  ## initiate a vector of length set
  cov_diff <- numeric(length(short_UTR_abun)) ## initiate with 0s
  ss <- max(search_point_START, 1)
  se <- min(search_point_END, .l)

  for (i in ss:se) {
    cov_diff_tmp <- .ele
    cov_diff_tmp <- cov_diff_tmp - long_UTR_abun[i]
    cov_diff_tmp[1:i] <- cov_diff_tmp[1:i] - short_UTR_abun[i]
    cov_diff[i] <- mean(cov_diff_tmp^2)
  }
  cov_diff
}
