#' search proximal CPsites
#'
#' @param CPs output from [search_distalCPs()] or [adjust_distalCPs()]
#' @param curr_UTR GRanges for current 3' UTR
#' @param window_size window size
#' @param MINSIZE MINSIZE for short form
#' @param cutEnd a numeric(1) between 0 and 1. The percentage of nucleotides
#'   should be removed from the end before search, 0.1 means 10 percent.
#' @param search_point_START An integer, specifying the start position to 
#'   calculate MSE
#' @param search_point_END An integer, specifying end position to calculate MSE
#'
#' @return a list
#' @seealso [adjust_proximalCPs()], [polish_CPs()],[adjust_proximalCPsByPWM()],
#'   [adjust_proximalCPsByNBC()], [get_PAscore()], [get_PAscore2()]
#' @keywords internal
#' @author Jianhong Ou

search_proximalCPs <- function(CPs, 
                              curr_UTR,
                              window_size, 
                              MINSIZE,
                              cutEnd,
                              search_point_START,
                              search_point_END) {
  dCPs <- CPs$dCPs
  chr.cov.merge <- CPs$chr.cov.merge
  next.exon.gap <- CPs$next.exon.gap
  annotated.utr3 <- CPs$annotated.utr3
  saved.id <- dCPs$length <- sapply(annotated.utr3, length)
  saved.proximal.apa <- mapply(function(.ele, .len) {
    as.numeric(gsub("^.*_SEP_", "", names(.ele)[.len]))
  }, annotated.utr3, saved.id)
  flag <- dCPs$distalCP > window_size
  dCPs$length[flag] <- dCPs$length[flag] + dCPs$distalCP[flag]
  dCPs$type <- ifelse(flag, "novel distal", "novel proximal")
  dCPs$type[grepl("proximalCP", dCPs$annotatedProximalCP) &
    !flag] <- "annotated proximal"
  dist_apa <- function(d, id) {
    ifelse(id > 0, as.numeric(rownames(d)[id]), 0)
  }
  chr.cov.merge <- lapply(chr.cov.merge, function(.ele) {
    rownames(.ele) <- gsub("^.*_SEP_", "", rownames(.ele))
    .ele
  })
  dCPs$Predicted_Distal_APA <- mapply(dist_apa, chr.cov.merge, dCPs$length)
  chr.cov.merge <- mapply(function(d, len) {
    if (len > 0) {
      d[1:len, , drop = FALSE]
    } else {
      d[FALSE, , drop = FALSE]
    }
  }, chr.cov.merge, dCPs$length, SIMPLIFY = FALSE)
  proximalCP <- sapply(curr_UTR, function(.ele) {
    grepl("proximalCP", .ele[.ele$feature == "utr3"]$annotatedProximalCP[1])
  })
  Predicted_Proximal_APA <- vector("list", length = nrow(dCPs))
  fit_value <- vector("list", length = nrow(dCPs))
  Predicted_Proximal_APA_rev <- vector("list", length = nrow(dCPs))
  fit_value_rev <- vector("list", length = nrow(dCPs))
  dCPs$fit_value <- NA
  if (sum(proximalCP) > 0) {
    Predicted_Proximal_APA[proximalCP] <-
      lapply(curr_UTR[proximalCP], function(.ele) {
        as.integer(unlist(strsplit(
          .ele[.ele$feature == "utr3"]$annotatedProximalCP[1], "_"
        )[2]))
      })
  }
  if (!is.na(cutEnd)) {
    if (cutEnd < 1) {
      chr.cov.merge <- lapply(chr.cov.merge, function(.ele) {
        .ele[1:floor((nrow(.ele) - 1) * (1 - cutEnd)), ,
          drop = FALSE]})
    } else {
      chr.cov.merge <- lapply(chr.cov.merge, function(.ele) {
        .ele[1:max(nrow(.ele) - 1 - floor(cutEnd), 1), ,
          drop = FALSE]})
    }
  }

  minStartPos <- dCPs$length >= max(c(search_point_START, MINSIZE))
  len <- sapply(chr.cov.merge, nrow)
  search_point_END <- rep(abs(search_point_END), nrow(dCPs))
  search_point_end <- ifelse(is.na(search_point_END),
    len - 1,
    ifelse(search_point_END < 1,
      ceiling(len * (1 - search_point_END)),
      floor(len - search_point_END)
    )
  )
  flag <- minStartPos & (search_point_end > search_point_START) & (!proximalCP)

  ## if no annotated proximal CPsites, find them
  fit_value[flag] <- mapply(function(.ele, search_point_END) {
   cov_diff <- apply(.ele, 2, calculate_mse,
                 search_point_START = search_point_START,
                 search_point_END = search_point_END
    )
    cov_diff <- rowMeans(cov_diff)
  }, chr.cov.merge[flag], search_point_end[flag], SIMPLIFY = FALSE)
  
  Predicted_Proximal_APA[flag] <-
    mapply(function(cov_diff, search_point_END, savedID) {
      idx <- find_valley(cov_diff, search_point_START,
                         search_point_END,
                         n = 5, savedID = savedID
      )
      if (search_point_START < MINSIZE) {
        idx <- idx[idx != search_point_START]
      }
      idx
    }, fit_value[flag], search_point_end[flag], saved.id[flag],
    SIMPLIFY = FALSE)
  
  idx1 <- lapply(Predicted_Proximal_APA, `[`, 1)
  idx1[sapply(idx1, length) == 0] <- NA
  idx1 <- unlist(idx1)

  list(dCPs = dCPs, chr.cov.merge = chr.cov.merge,
       next.exon.gap = next.exon.gap,
       annotated.utr3 = annotated.utr3,
       flag = flag, fit_value = fit_value,
       Predicted_Proximal_APA = Predicted_Proximal_APA,
       saved.id = saved.id, idx1 = idx1)
}
