#' search distal CP sites
#'
#' search distal CP sites
#'
#' @param chr.cov.merge merged coverage data for a given chromosome
#' @param conn_next_utr3 A logical(1) vector, indicating whether joint to next
#'   3UTR or not (used by [remove_convergentUTR3s()])
#' @param curr_UTR GRanges of 3' UTR for a given chromosome
#' @param window_size An integer(1) vector, the window size for novel distal or 
#'   proximal CP site searching. default: 100.
#' @param depth.weight A named vector. One element of an output of
#'   [setup_CPsSearch()] for coverage depth weight, which is the output of
#'   [get_depthWeight()]
#' @param long_coverage_threshold An integer(1) vector, specifying the cutoff 
#'   threshold of coverage for the terminal of long form 3' UTRs. If the coverage
#'   of first 100 nucleotides is lower than coverage_threshold, that transcript
#'   will be not considered for further analysis. Default, 2.
#' @param background A character(1) vector, the range for calculating cutoff
#'   threshold of local background. It can be "same_as_long_coverage_threshold",
#'   "1K", "5K","10K", or "50K".
#' @param z2s one element of an output of [setup_CPsSearch()] for Z-score cutoff
#'   values, which is the output of [get_zScoreCutoff()]
#'   
#' @return a list
#' #' \itemize{\item dCPs, a data frame converted from GRanges
#'             \item chr.cov.merge, depth-normalized sample/condition specific 
#'             coverage
#'             \item next.exon.gap, all-in-one collapsed, refined next.exon.gap
#'              coverage
#'             \item annotated.utr3,all-in-one collapsed coverage for annotated 
#'             proximal UTRs
#'         }
#' @seealso [get_PAscore2()]
#' @import S4Vectors
#' @keywords internal
#' @author Jianhong Ou

search_distalCPs <- function(chr.cov.merge,
                             conn_next_utr3,
                             curr_UTR,
                             window_size,
                             depth.weight,
                             long_coverage_threshold,
                             background,
                             z2s) {
  find_distalCPs <- function(chr.cov.merge.ele,
                        conn_next_utr,
                        curr_UTR.ele) {
    ## sequencing depth-normalized sample/condition-specific utr3_total_coverage
    ## next.exon.gap total coverage
    chr.cov.merge.ele <-
      t(t(chr.cov.merge.ele) / depth.weight[colnames(chr.cov.merge.ele)])
    ## all-in-one collapsed utr3_total_coverage
    .ele <- rowSums(chr.cov.merge.ele)
    # By Jianhong, if there are reads covered to the next.exon.gap the proximal 
    # CP site should be the known-utr3-end and distal site should be the end
    # of gap.
    # But per Haibo, this might be NOT correct in many cases.
    next.exon.gap <- .ele[grepl("next.exon.gap", names(.ele))]
    # remove the gaps with 0 coverage and their width > window_size
    next.exon.gap.rle <- Rle(next.exon.gap)
    id <- which(runLength(next.exon.gap.rle) > window_size &
      runValue(next.exon.gap.rle) == 0)
    if (length(id) >= 1) {
      conn_next_utr <- FALSE
      id <- id[1] - 1
      if (id > 0) {
        id <- sum(runLength(next.exon.gap.rle)[1:id])
        next.exon.gap <- next.exon.gap[1:id]
      } else {
        next.exon.gap <- numeric(0)
      }
    }
    ## stop if any connected two windows different from 30 times,
    ## to separate connected two UTRs or strange peaks.
    if (length(next.exon.gap) > 1) {
      next.exon.gap.ids <- rep(1:ceiling(length(next.exon.gap) / window_size),
                               each = window_size)[1:length(next.exon.gap)]
      next.exon.gap.split <- split(next.exon.gap, next.exon.gap.ids)
      next.exon.gap.split <- sapply(next.exon.gap.split, mean)
      next.exon.gap.split.diff <-
        log2(next.exon.gap.split[-1] + 0.01) -
        log2(next.exon.gap.split[-length(next.exon.gap.split)] + 0.01)
      
      #coverage ratio greater than 30 between two adjacent windows?
      id <- which(abs(next.exon.gap.split.diff) > log2(30))
      ### drop by background by window
      ## the gap == 2X window size with low coverage
      cont.num <- function(x) {
        if (length(x) == 1) {
          return(x)
        }
        ## the index with two neighboring windows with coverage < specified 
        ## background
        r <- which(diff(x) == 1)
        if (length(r) > 0) {
          return(x[r[1] + 1])
        } else {
          return(x[length(x)])
        }
      }
      if (background != "same_as_long_coverage_threshold") {
        z2_threshold <- z2s[names(curr_UTR.ele[curr_UTR.ele$feature == "utr3"])]
        if (!is.na(z2_threshold)) {
          z2 <- which(next.exon.gap.split < z2_threshold)
          if (length(z2) > 0) {
            z2 <- cont.num(z2)
            id <- min(c(z2, id)) ## end of two adjacent window with < background
                                 ## coverage or abs(coverage ratio) > 30.
          }
        }
      } else {
        ## background == long_coverage_threshold
        z2 <- which(next.exon.gap.split < long_coverage_threshold)
        if (length(z2) > 0) {
          z2 <- cont.num(z2)
          id <- min(c(z2, id))
        }
      }
      if (length(id) > 0) {
        conn_next_utr <- FALSE
        id <- id[1] - 1
        if (id > 0) {
          next.exon.gap <- next.exon.gap[1:(window_size * id)]
        } else {
          next.exon.gap <- numeric(0)
        }
      }
    }
    ## remove utr3---___---utr3, need to improve.
    if (conn_next_utr && length(next.exon.gap) > 50) {
      next.exon.gap <- remove_convergentUTR3s(next.exon.gap)
    }
    annotated.utr3 <- .ele[grepl("utr3", names(.ele))]
    utr3start <- as.numeric(gsub("^.*_SEP_", "", names(annotated.utr3)[1]))
    utr3end <- as.numeric(gsub("^.*_SEP_", "",
                          names(annotated.utr3)[length(annotated.utr3)]))
    # last.annotated.utr3 <- annotated.utr3[length(annotated.utr3)]
    
    ## drop low coverage from the 3' end if basewise coverage < 
    ## long_coverage_threshold
    if (length(next.exon.gap) > 5) {
      ## make it simple, smooth and cutoff
      fit.loess <- loess.smooth(1:length(next.exon.gap),
                                next.exon.gap,
                                evaluation = length(next.exon.gap),
                                family = "gaussian",
                                span = .75, degree = 1)

      next.exon.gap.abun <- fit.loess$y
      ## check smoothed coverage
      id <- which(next.exon.gap.abun >= long_coverage_threshold)
      if (length(id) > 0) {
        ## check real coverage
        next.exon.gap.abun <- next.exon.gap[1:id[length(id)]]
        id <- which(next.exon.gap.abun >= long_coverage_threshold)
      }
      if (length(id) > 0) {
        next.exon.gap.end.pos <- id[length(id)]
      } else {
        next.exon.gap.end.pos <- 0
      }
    } else {
      next.exon.gap.end.pos <- 0
    }

    info <- as.data.frame(curr_UTR.ele)[1, ]
    info$utr3start <- utr3start
    info$utr3end <- utr3end
    info$distalCP <- next.exon.gap.end.pos
    list(info = info,
         cov = chr.cov.merge.ele,
         gap = next.exon.gap,
         annotated.utr3 = annotated.utr3)
  }
  
  distalCPs <- mapply(find_distalCPs,
                      chr.cov.merge,
                      conn_next_utr3,
                      curr_UTR,
                      SIMPLIFY = FALSE)
  dCPs <- do.call(rbind, lapply(distalCPs, `[[`, "info"))
  ## depth-normalized sample/condition-specific coverage (utr3 + gap)
  chr.cov.merge <- lapply(distalCPs, `[[`, "cov")
  ##refined next.exon.gap coverage
  next.exon.gap <- lapply(distalCPs, `[[`, "gap") 
  ## annotated 3' UTR coverage
  annotated.utr3 <- lapply(distalCPs, `[[`, "annotated.utr3")
  
  list(dCPs = dCPs,   ## see info above
       chr.cov.merge = chr.cov.merge,
       next.exon.gap = next.exon.gap,
       annotated.utr3 = annotated.utr3)
}
