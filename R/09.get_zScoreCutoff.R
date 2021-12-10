#' Calculate local background cutoff value
#'
#' calculate local background z-score cutoff 
#'
#' @param background A character(1) vector, indicating how background coverage
#'   is defined.
#' @param chr.introns An object of [GenomicRanges::GRanges-class] for introns 
#'   of a give chromosome/scaffold
#' @param chr.totalCov total coverage for a given chromosome/scaffold, an output 
#'   from [get_totalCov()] for a given chromosome/scaffold
#' @param chr.utr3 An object of [GenomicRanges::GRanges-class], an element of 
#'   the output of [extract_UTR3Anno()] for a given chromosome/scaffold
#' @param seqname A character(1), the name of a chromosome/scaffold
#' @param z Z score cutoff value
#'
#' @return A named numeric vector containing local background Z-score cutoff
#'   values. The names are GRanges's name for 3' UTRs. 
#' @keywords internal
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @author Jianhong Ou, Haibo Liu

get_zScoreCutoff <- function(background,
                             chr.introns, 
                             chr.totalCov,
                             chr.utr3,
                             seqname,
                             z = 2) {
  ## output of get_totalCov() is a list of list:
  ## chr.totalCov[[seqname]][[condition]]
  chr.totalCov <- chr.totalCov[[1]]
  if (background == "same_as_long_coverage_threshold" || 
      length(chr.introns) < 1) {
    return(0)
  }
  background <- switch(background,
                       "1K" = 1000,
                       "5K" = 5000,
                       "10K" = 10000,
                       "50K" = 50000,
                       1000)
  
  curr_range <- chr.utr3[chr.utr3$feature == "utr3"]
  curr_str <- as.character(strand(curr_range)) == "+"
  start(curr_range) <- ifelse(curr_str,
                              start(curr_range) - background - 1,
                              end(curr_range))
  width(curr_range) <- background
  
  ol <- findOverlaps(curr_range, chr.introns)
  if (length(ol) < 1) {
    return(0)
  }
  chr.introns <- chr.introns[subjectHits(ol)]
  chr.introns$id <- names(curr_range)[queryHits(ol)]
  
  #no introns for some seqnames
  get_cvg <- function(.intron, .cov) {
    .cvg <-
      lapply(.cov, function(.cv) {
        .cv <- Views(.cv,
                     start = start(.intron),
                     end = end(.intron))
        .cv <- viewApply(.cv, as.integer, simplify = FALSE)
      })
    
    # cvg is a list of length(#number_sample) of list of length(#number_introns) 
    # correct Jianhong's bug here
    .cvg1 <- vector("list", length(.cvg[[1]]))
    for (k in seq_along(.cvg1))
    {
      .cvg1[[k]] <- list()
    }
    
    for (i in names(.cvg)) {
      for (j in seq_along(.cvg[[1]])) {
        .cvg1[[j]][[i]] <- .cvg[[i]][[j]]
      }
    }
    .cvg1 <- lapply(.cvg1, function(.e) {
      colSums(do.call("rbind", .e))
    })
    .intron$cvg <- .cvg1
    .intron
  }
  
  cvg <- get_cvg(.intron = chr.introns, .cov = chr.totalCov)

  cvg <- split(cvg$cvg, cvg$id)
  cvg <- sapply(cvg, function(.ele) {
    .ele <- unlist(.ele)
    mu <- mean(.ele, na.rm = TRUE)
    std <- sd(.ele, na.rm = TRUE)
    z * std + mu ## assuming normal distribution?
  })
  
  cvg <- cvg[match(names(curr_range), names(cvg))]
  names(cvg) <- names(curr_range)
  cvg[is.na(cvg)] <- 0
  cvg
}