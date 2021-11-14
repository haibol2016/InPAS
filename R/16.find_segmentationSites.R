#' Find optimal segmentation sites by minimizing MSE
#'
#' Find optimal segmentation sites by minimizing MSE
#'
#' @param .ele 3' UTR coverage
#' @param search_point_START start position to calculate
#' @param search_point_END end position to calculate
#' @param n the length of output. If set to -1, return all candidate CP sites
#' @param savedID the proximal CPsites
#'
#' @return a list of MSE and candidate APA sites
#' @keywords internal
#' @author Jianhong Ou

find_segmentationSites <- function(.ele, 
                                search_point_START,
                                search_point_END,
                                n = 1, 
                                savedID = NA) {
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
  cov_diff <- numeric(length(short_UTR_abun))
  ss <- max(search_point_START, 1)
  se <- min(search_point_END, .l)

  for (i in ss:se) {
    cov_diff_tmp <- .ele
    cov_diff_tmp <- cov_diff_tmp - long_UTR_abun[i]
    cov_diff_tmp[1:i] <- cov_diff_tmp[1:i] - short_UTR_abun[i]
    cov_diff[i] <- mean(cov_diff_tmp^2)
  }
  idx <- find_valley(cov_diff, ss, se, n, savedID)
  return(list(cov_diff = cov_diff, idx = idx))
}
