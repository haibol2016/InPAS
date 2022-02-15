#' search proximal CPsites
#'
#' @param CPs output from [search_distalCPs()] or [adjust_distalCPs()]
#' @param curr_UTR GRanges for current 3' UTR
#' @param window_size window size
#' @param MINSIZE MINSIZE for short form
#' @param cutEnd A numeric(1) between 0 and 1 or an integer(1) greater than 1,
#'   specifying the percentage of or the number of nucleotides should be removed
#'   from the end before search for proximal CP sites, 0.1 means 10 percent.
#'   It is recommended to use an integer great than 1, such as 200, 400 or 600,
#'   because read coverage at 3' extremities is determined by fragment size due
#'   to RNA fragmentation and size selection during library construction.
#' @param search_point_START An integer, specifying the start position to
#'   calculate MSE
#' @param search_point_END A numeric(1) between 0 and 1 or an integer(1) greater
#'  than 1, specifying the percentage of or the number of nucleotides should not
#'  be excluded from the end to calculate MSE.
#' @param filter.last A logical(1), whether to filter out the last valley, which
#'   is likely the 3' end of the longer 3' UTR if no novel distal CP site is
#'   detected and the 3' end excluded by setting cutEnd/search_point_END is small.
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
                               cutEnd = 0,
                               search_point_START,
                               search_point_END = NA,
                               filter.last = TRUE) {
  dCPs <- CPs$dCPs
  chr.cov.merge <- CPs$chr.cov.merge
  utr_len <- sapply(CPs$final.utr3, length)

  ## annotated proximal CP sites per 3' UTR
  ## known candidate proximal CP sites recorded in UTR3 annotation
  utr_end <- ifelse(dCPs$strand == "-", dCPs$start, dCPs$end)
  saved.id <-
    mapply(function(anno.cp, truncated, end) {
      if (!truncated) {
        if (anno.cp == "unknown") {
          return(end)
        } else {
          annotated.proximal <- as.numeric(unlist(strsplit(gsub(
            "proximalCP_", "",
            anno.cp
          ), "_")))
          annotated.proximal <- c(annotated.proximal, end)
          return(annotated.proximal)
        }
      } else {
        if (anno.cp == "unknown") {
          return(NA)
        } else {
          annotated.proximal <- as.numeric(unlist(strsplit(gsub(
            "proximalCP_", "",
            anno.cp
          ), "_")))
          return(annotated.proximal)
        }
      }
    }, dCPs$annotatedProximalCP, dCPs$truncated,
    utr_end,
    SIMPLIFY = FALSE
    )
  names(saved.id) <- rownames(dCPs)
  flag <- dCPs$distalCP > 0

  ## trimming chr.cov.merge according to current 3' utr length
  chr.cov.merge <- mapply(function(d, len) {
    if (len > 0) {
      d[1:len, , drop = FALSE]
    } else {
      d[FALSE, , drop = FALSE]
    }
  }, chr.cov.merge, utr_len, SIMPLIFY = FALSE)

  if (!is.na(cutEnd)) {
    if (cutEnd < 1) { ## By proportion
      chr.cov.merge <- lapply(chr.cov.merge, function(.ele) {
        .ele[1:floor((nrow(.ele) - 1) * (1 - cutEnd)), ,
          drop = FALSE
        ]
      })
    } else { ## actual length
      chr.cov.merge <- lapply(chr.cov.merge, function(.ele) {
        .ele[1:max(nrow(.ele) - 1 - floor(cutEnd), 1), ,
          drop = FALSE
        ]
      })
    }
  }

  minStartPos <- utr_len >= max(c(search_point_START, MINSIZE))
  len <- sapply(chr.cov.merge, nrow)
  search_point_END <- rep(abs(search_point_END), nrow(dCPs))
  search_point_end <- ifelse(is.na(search_point_END),
    len - 1,
    ifelse(search_point_END < 1,
      ceiling(len * (1 - search_point_END)),
      floor(len - search_point_END)
    )
  )
  flag <- minStartPos & (search_point_end > search_point_START)

  Predicted_Proximal_APA <- vector("list", length = nrow(dCPs))
  fit_value <- vector("list", length = nrow(dCPs))
  dCPs$fit_value <- NA
  dCPs$Predicted_Proximal_APA <- NA
  dCPs$Predicted_Proximal_APA_Type <- NA

  ## find proximal CP sites if the utr3 is valid in term of length
  fit_value[flag] <- mapply(function(.ele, search_point_END) {
    cov_diff <- apply(.ele, 2, calculate_mse,
      search_point_START = search_point_START,
      search_point_END = search_point_END
    )
    cov_diff <- rowMeans(cov_diff)
  }, chr.cov.merge[flag],
  search_point_end[flag],
  SIMPLIFY = FALSE
  )

  Predicted_Proximal_APA[flag] <-
    mapply(function(cov_diff, search_point_END) {
      idx <- find_valleyBySpline(cov_diff,
        search_point_START,
        search_point_END,
        n = 1,  #be more accurate, only return the one with minimal MSE
        filter.last = filter.last
      )
      if (search_point_START < MINSIZE) {
        idx <- idx[idx != search_point_START]
      }
      idx
    }, fit_value[flag],
    search_point_end[flag],
    SIMPLIFY = FALSE
    )

  list(
    dCPs = dCPs,
    chr.cov.merge = chr.cov.merge,
    saved.id = saved.id,
    flag = flag,
    fit_value = fit_value,
    Predicted_Proximal_APA = Predicted_Proximal_APA
  )
}
