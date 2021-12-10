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
#' @param chr.utr3 An object of [GenomicRanges::GRanges-class]. It must be an 
#'   element of the output of [extract_UTR3Anno()] for a given chromosome.
#' @param chr.totalCov total coverage for each condition of a given chromosome. It
#'   must be an output of [get_totalCov()]
#' @param gcCompensationensation GC compensation vector. Not support yet.
#' @param mappabilityCompensation mappability compensation vector. Not support
#'   yet.
#' @param FFT Use FFT smooth or not.
#' @param fft.sm.power the cut-off frequency of FFT smooth.
#'
#' @return path to a file storing the UTR3 total coverage for a given 
#'   chromosome/scaffold
#' @keywords internal
#' @import GenomicRanges
#' @importFrom IRanges Views viewApply IRanges viewMeans
#' @author Jianhong Ou

get_UTR3TotalCov <- function(chr.utr3, 
                             chr.totalCov,
                             gcCompensation = NA,
                             mappabilityCompensation = NA,
                             FFT = FALSE,
                             fft.sm.power = 20) {
  ## output of get_totalCov() is a list of list:
  ## chr.totalCov[[seqname]][[condition]]
  chr.totalCov <- chr.totalCov[[1]]
  ## memory consuming if gcCompensation, and mappabilityCompensation are 
  ## used for coverage correction
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
      view <- lapply(view, function(.ele) {
        apply(.ele, 2, fft.smooth, p = fft.sm.power)
      })
    } else {
      view
    }
  }
  
  view <- utr3_totalCov(.utr = chr.utr3,
                        .cvg = chr.totalCov, 
                        gcComp = gcCompensation,
                        mapComp = mappabilityCompensation)
  view
}
