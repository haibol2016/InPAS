#' Calculate local background cutoff value
#'
#' calculate local background z-score cutoff 
#'
#' @param background A character(1) vector, indicating how background coverage
#'   is defined.
#' @param introns An object of [GenomicRanges::GRanges-class] for introns
#' @param totalCov total coverage, a list of output from [get_totalCov()] for
#'   each chromosome
#' @param utr3 An object of [GenomicRanges::GRangesList-class], output of
#'   [extract_UTR3Anno()]
#' @param z Z score cutoff value
#'
#' @return A named numeric vector containing local background Z-score cutoff
#'   values. The names are GRanges's name for 3' UTRs. 
#' @keywords internal
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @author Jianhong Ou, Haibo Liu

get_zScoreCutoff <- function(background,
                            introns, 
                            totalCov,
                            utr3, 
                            z = 2) {
  if (background == "same_as_long_coverage_threshold" || 
      length(introns) < 1) {
    return(0)
  }
  
  ## convert GRangesList to GRanges
  utr3 <- unlist(utr3)
  curr_range <- utr3[utr3$feature == "utr3"]
  seqnames <-  sort(intersect(as.character(seqnames(curr_range)), 
                              names(totalCov)))
  
  background <- switch(background,
                       "1K" = 1000,
                       "5K" = 5000,
                       "10K" = 10000,
                       "50K" = 50000,
                       1000)
  
  curr_str <- as.character(strand(curr_range)) == "+"
  start(curr_range) <- ifelse(curr_str,
                              start(curr_range) - background - 1,
                              end(curr_range))
  width(curr_range) <- background
  
  introns <- introns[as.character(seqnames(introns)) %in% seqnames]
  ol <- findOverlaps(curr_range, introns)
  if (length(ol) < 1) {
    return(0)
  }
  introns <- introns[subjectHits(ol)]
  introns$id <- names(curr_range)[queryHits(ol)]
  
  intron_containing_seqnames <- as.character(seqnames(introns))
  introns.s <- split(introns, intron_containing_seqnames)
  intron_containing_seqnames <- unique(intron_containing_seqnames)
  
  #no introns for some seqnames
  cvg <- mapply(function(.intron, .cov) {
    .cvg <-
      lapply(.cov, function(.cv) {
        .cv <- Views(.cv,
                     start = start(.intron),
                     end = end(.intron)
        )
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
  }, introns.s[intron_containing_seqnames], 
  totalCov[intron_containing_seqnames], SIMPLIFY = FALSE)
  
  cvg <- unlist(GRangesList(cvg))
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