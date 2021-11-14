#' Compensate the coverage with GC-content or mappability
#'
#' @param view A list of view object
#' @param comp A numeric vector of weight for GC composition or mappability
#' @param start An integer vector, starting coordinates
#' @param end  An integer vector, end coordinates
#'
#' @return a list of GC composition or mappability corrected coverage
#' @keywords internal
#' @author Jianhong Ou

compensation <- function(view, comp, start, end) {
  mapply(function(.ele, .s, .e) {
    .ele * comp[.s:.e]
  }, view, start, end, SIMPLIFY = FALSE)
}


#' Smoothing using Fast Discrete Fourier Transform
#'
#' @param sn a real or complex array containing the values to be transformed.
#'   see [stats::fft()]
#' @param p An integer(1), fft smoothing power
#'
#' @return a numeric vector, the real part of inverse fft-transformed signal
#' @keywords internal
#' @author Jianhong Ou

fft.smooth <- function(sn, p) {
  if (length(sn) <= p) {
    return(sn)
  }
  sn.fft <- fft(sn)
  sn.fft[p:length(sn.fft)] <- 0 + 0i
  sn.ifft <- fft(sn.fft, inverse = TRUE) / length(sn.fft)
  Re(sn.ifft)
}

#' extract coverage of 3' UTR for CP sites prediction
#'
#' extract 3' UTR coverage from totalCov according to the 
#'   [GenomicRanges::GRanges-class] object utr3.
#'
#' @param utr3 An object of [GenomicRanges::GRangesList-class]. It must be the
#'   output of [extract_UTR3Anno()]
#' @param totalCov total coverage for each sample. It must be a list of the 
#'   output of [get_totalCov()]
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam()
#' @param gcCompensation GC compensation vector. Not support yet.
#' @param mappabilityCompensation mappability compensation vector. Not support
#'   yet.
#' @param FFT Use FFT smooth or not.
#' @param fft.sm.power the cut-off frequency of FFT smooth.
#'
#' @return A list containing total coverage for each 3' UTR. 
#'   \describe{
#'      \item{seqname}{chromosome/scaffold name}
#'      \describe{
#'      \item{transcript1}{data matrix containing summarized View
#'            for transcript1}
#'      \item{transcript2}{data matrix containing summarized View
#'            for transcript2}
#'      \item{transcriptN}{data matrix containing summarized View
#'            for transcriptN}
#'         }
#'      }
#' @keywords internal
#' @import GenomicRanges
#' @importFrom IRanges Views viewApply IRanges viewMeans
#' @importFrom BiocParallel bptry bpok bplapply bpparam bpmapply
#' @author Jianhong Ou

get_UTR3TotalCov <- function(utr3, totalCov,
                             BPPARAM = NULL,
                             gcCompensation = NA,
                             mappabilityCompensation = NA,
                             FFT = FALSE,
                             fft.sm.power = 20) {
  seqnames <- sort(intersect(names(utr3), names(totalCov)))
  
  ## memory consuming if gcComp, and mapComp are used for coverage correction
  utr3_totalCov <- function(.utr, .cvg, gcComp, mapComp) {
    strand <- as.character(strand(.utr))
    start <- start(.utr)
    end <- end(.utr)
    maxEnd <- max(end)
    for (i in 1:length(.cvg)) {
      ## why this can happen?
      if (maxEnd > length(.cvg[[i]])) {
        .cvg[[i]] <-
          append(
            .cvg[[i]],
            ## no  + 1 by Haibo
            rep(0, maxEnd - length(.cvg[[i]]))
          )
      }
    }
    
    view <- lapply(.cvg, Views,
                   start = start,
                   end = end, names = names(.utr)
    )
    view <- sapply(view, function(.view) {
      viewApply(.view, as.integer)
    })
    if (is.list(view[1, ])) {
      view <- apply(view, 1, function(.ele) {
        do.call(cbind, .ele)
      })
    } else {
      view <- list(view)
      names(view) <- names(.utr)
    }
    if (length(view) > 0) view <- view[sapply(view, mean) > 0]
    if (!is.na(gcComp[1]) && length(gcComp) == length(.cvg[[1]])) {
      view <- compensation(view, gcComp, start, end)
    }
    if (!is.na(mapComp[1]) && length(mapComp) == length(.cvg[[1]])) {
      view <- compensation(view, mapComp, start, end)
    }
    if (FFT) {
      lapply(view, function(.ele) {
        apply(.ele, 2, fft.smooth, p = fft.sm.power)
      })
    } else {
      view
    }
  }
  
  if(is.null(BPPARAM)){
    utr3TotalCov <- mapply(utr3_totalCov, utr3[seqnames], 
                           totalCov[seqnames],
                           gcCompensation[seqnames],
                           mappabilityCompensation[seqnames],
                           SIMPLIFY = FALSE)
  } else {
    utr3TotalCov <- bptry(bpmapply(utr3_totalCov, utr3[seqnames], 
                           totalCov[seqnames],
                           gcCompensation[seqnames],
                           mappabilityCompensation[seqnames],
                           SIMPLIFY = FALSE,
                           BPPARAM = BPPARAM))
    while (!all(bpok(utr3TotalCov))){
      utr3TotalCov <- bptry(bpmapply(utr3_totalCov, utr3[seqnames], 
                                     totalCov[seqnames],
                                     gcCompensation[seqnames],
                                     mappabilityCompensation[seqnames],
                                     SIMPLIFY = FALSE,
                                     BPREDO = utr3TotalCov,
                                     BPPARAM = BPPARAM))
    }
  }
  utr3TotalCov
}
