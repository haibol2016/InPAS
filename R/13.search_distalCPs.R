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
    annotated.utr3 <- .ele[grepl("utr3", names(.ele))]
    
    ## due to cutStart = 10, the first 10 bases are missing.
    utr3start <- as.numeric(gsub("^.*_SEP_", "", 
                                 names(annotated.utr3)[1]))
    utr3end <- as.numeric(gsub(
        "^.*_SEP_", "",
        names(annotated.utr3)[length(annotated.utr3)]
    ))
    
    ## prepare output data.frame
    info <- as.data.frame(curr_UTR.ele)[1, ]
    info$utr3start <- utr3start
    info$utr3end <- utr3end
    info$preliminary_distal_APA <- 0
    info$NBC_adjusted_distal_APA <- NA 
    info$Predicted_Distal_APA <- NA
    info$Predicted_Distal_APA_type <- NA

    # detect gaps > window_size with 0 coverage and their width 
    # and trimming the collapsed next.exon.gap coverage.
    refine_gap <- function(gap.cov, window_size,
                           conn_next_utr,
                           n_window = 1,
                           is_first = TRUE) {
      gap.cov.rle <- Rle(gap.cov)
      id <- which(runLength(gap.cov.rle) > n_window * window_size &
        runValue(gap.cov.rle) == 0)
      ## find first uncovered region
      if (length(id) >= 1) {
        conn_next_utr <- FALSE
        id <- {
          if (is_first) {
            id[1] - 1
          } else {
            id[length(id)] - 1
          }
        }
        if (id > 0) {
          tot_len <- sum(runLength(gap.cov.rle)[1:id])
          gap.cov <- gap.cov[1:tot_len]
        } else {
          gap.cov <- numeric(0)
        }
      }
      list(
        conn_next_utr = conn_next_utr,
        gap.cov = gap.cov
      )
    }

    if (length(next.exon.gap) > 0) {
        # If the next.exon.gap is shared by two convergent 3' UTRs, cut at the first
        # gap with 0 coverage,otherwise, cut at the first gap with 0 coverage
        refined_gap <- refine_gap(
            gap.cov = next.exon.gap,
            window_size = window_size,
            conn_next_utr = conn_next_utr,
            n_window = 1,
            is_first = TRUE
        )
        next.exon.gap <- refined_gap$gap.cov
        conn_next_utr <- refined_gap$conn_next_utr
    } 
        ## By Jianhong, cut if any connected two windows different from 30 times,
        ## to separate connected two UTRs or strange peaks. This cut is not
        ## implemented any more.
    if (length(next.exon.gap) > 1) {
      if (background != "same_as_long_coverage_threshold") {
        threshold <- z2s[names(curr_UTR.ele[curr_UTR.ele$feature == "utr3"])]
        if (is.na(threshold)) {
          threshold <- long_coverage_threshold
        }
      } else {
        threshold <- long_coverage_threshold
      }
      next.exon.gap.thresholded <- next.exon.gap
      
      # if coverage < threshold, set to 0
      next.exon.gap.thresholded[next.exon.gap.thresholded < threshold] <- 0
      
      ## depending on if this next.exon.gao is shared by two convergent
      ## transcripts
      refined_gap <- refine_gap(
        gap.cov = next.exon.gap.thresholded,
        window_size = window_size,
        conn_next_utr = conn_next_utr,
        n_window = 1,
        is_first = TRUE
      )
      next.exon.gap.thresholded <- refined_gap$gap.cov
      conn_next_utr <- refined_gap$conn_next_utr
      if (length(next.exon.gap.thresholded) > 0) {
        next.exon.gap <- next.exon.gap[1:length(next.exon.gap.thresholded)]
      } else {
        next.exon.gap <- numeric(0)
      }
    }
    ## remove utr3---___---utr3, need to improve.
    if (conn_next_utr && length(next.exon.gap) > 50) {
      ## find split position by the deepest valley in coverage profile
      next.exon.gap <- InPAS:::remove_convergentUTR3s(next.exon.gap)
      conn_next_utr <- FALSE
    }
    
    ## By Jianhong, trim low coverage from the 3' end if basewise coverage <
    ## long_coverage_threshold.
    next.exon.gap.end.pos <- 0
    if (length(next.exon.gap) > 5) {
      ## make it simple, smooth and cutoff
      fit.loess <- loess.smooth(1:length(next.exon.gap),
                                next.exon.gap,
                                evaluation = length(next.exon.gap),
                                family = "gaussian",
                                span = .75, degree = 1
      )
      next.exon.gap.abun <- fit.loess$y
      id <- which(next.exon.gap.abun >= long_coverage_threshold)
      if (length(id) > 0) {
        next.exon.gap.abun <- next.exon.gap[1:id[length(id)]]
        id <- which(next.exon.gap.abun >= long_coverage_threshold)
      }
      if (length(id) > 0) {
        next.exon.gap.end.pos <- id[length(id)]
      } 
    }
    
    if (next.exon.gap.end.pos > 0) # longer 3' utr ending in next.exon.gap
    {
      distal.utr3.len <- length(annotated.utr3) + next.exon.gap.end.pos
      info$Predicted_Distal_APA_type <- "extended novel distal"
      ## final utr3 extended by valid next exon gap
      final.utr3 <- .ele[1:distal.utr3.len]
      ## absolute coordinates of distal CP sites
      info$preliminary_distal_APA <- as.numeric(gsub(
        "^.*_SEP_", "",
        names(final.utr3)[length(final.utr3)]
      ))
    } else { ## check annotated 3' UTR region
      original_len <- length(annotated.utr3)
      annotated.utr3.rle <- Rle(annotated.utr3)

      ## trimming uncovered tails > 100 bp
      if (runLength(annotated.utr3.rle)[nrun(annotated.utr3.rle)] >= 100 &&
        runValue(annotated.utr3.rle)[nrun(annotated.utr3.rle)] == 0) {
        annotated.utr3.rle <-
          Rle(
            values = runValue(annotated.utr3.rle)[-nrun(annotated.utr3.rle)],
            lengths = runLength(annotated.utr3.rle)[-nrun(annotated.utr3.rle)]
          )
      }
      len <- sum(runLength(annotated.utr3.rle))
      annotated.utr3 <- annotated.utr3[seq_len(len)]
      
      final.utr3 <- numeric(0)
      ## no 3' UTR of valid length
      info$preliminary_distal_APA <- -1
      info$Predicted_Distal_APA_type <- "undetectable"
      
      ## trimming by background threshold from 3' end
      if (length(annotated.utr3) > 0) {
        if (background != "same_as_long_coverage_threshold") {
          threshold <- z2s[names(curr_UTR.ele[curr_UTR.ele$feature == "utr3"])]
          if (is.na(threshold)) {
            threshold <- long_coverage_threshold
          }
        } else {
          threshold <- long_coverage_threshold
        }
        annotated.utr3.threshold <- annotated.utr3

        # if coverage < threshold, set to 0
        annotated.utr3.threshold[annotated.utr3.threshold < threshold] <- 0
        annotated.utr3.threshold.rle <- Rle(annotated.utr3.threshold)

        ## trimming annotated 3'UTR from 3' end conservatively
        ## remove trailing zeroes if there are more than 100 zeros
        if (runLength(annotated.utr3.threshold.rle)[nrun(annotated.utr3.threshold.rle)] >= 100 &&
          runValue(annotated.utr3.threshold.rle)[nrun(annotated.utr3.threshold.rle)] == 0) {
          annotated.utr3.threshold.rle <-
            Rle(
              values = runValue(annotated.utr3.threshold.rle)[-nrun(annotated.utr3.threshold.rle)],
              lengths = runLength(annotated.utr3.threshold.rle)[-nrun(annotated.utr3.threshold.rle)]
            )
        }
        trimmed_len <- sum(runLength(annotated.utr3.threshold.rle))

        if (trimmed_len > 0) {
          final.utr3 <- annotated.utr3[1:trimmed_len]
          ## just a little shortened
          if (trimmed_len >= original_len - 200) {
            info$Predicted_Distal_APA_type <- "known distal"
            info$preliminary_distal_APA <- {if (info$strand == "-"){info$start} 
                else {info$end}}
          } else {
            info$Predicted_Distal_APA_type <- "shortened novel distal"
            info$preliminary_distal_APA <- as.numeric(gsub(
              "^.*_SEP_", "",
              names(final.utr3)[length(final.utr3)]
            ))
          }
        }
      }
    }
    info$Predicted_Distal_APA <- info$preliminary_distal_APA
    return(list(
      info = info,
      cov = chr.cov.merge.ele,
      final.utr3 = final.utr3
    ))
  }

  distalCPs <- mapply(find_distalCPs,
    chr.cov.merge,
    conn_next_utr3,
    curr_UTR,
    SIMPLIFY = FALSE
  )
  dCPs <- do.call(rbind, lapply(distalCPs, `[[`, "info"))
  
  ## depth-normalized sample/condition-specific coverage (utr3 + gap)
  chr.cov.merge <- lapply(distalCPs, `[[`, "cov")
  ## annotated 3' UTR coverage
  final.utr3 <- lapply(distalCPs, `[[`, "final.utr3")

  list(
    dCPs = dCPs, ## see info above
    chr.cov.merge = chr.cov.merge,
    final.utr3 = final.utr3
  )
}
