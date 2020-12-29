#' Compensate the coverage with GC-content or mappability
#'
#' @param view A list of view object
#' @param comp A numeric vector of weight for GC composition or mappability
#' @param start An integer vector, starting coordinates
#' @param end  An integer vector, end coordinates
#'
#' @return a list of GC composition or mappability corrected coverage
#' @keywords internal
#'

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
#'

fft.smooth <- function(sn, p) {
  if (length(sn) <= p) {
    return(sn)
  }
  sn.fft <- fft(sn)
  sn.fft[p:length(sn.fft)] <- 0 + 0i
  sn.ifft <- fft(sn.fft, inverse = TRUE) / length(sn.fft)
  Re(sn.ifft)
}

#' extract coverage of 3UTR for CP sites prediction
#'
#' extract 3UTR coverage from totalCov according to the
#' [GenomicRanges::GRanges-class] object utr3.
#'
#' @param utr3 an [GenomicRanges::GRanges-class] object. It must be the output
#'   of [utr3Annotation()]
#' @param totalCov total coverage for each sample. It must be the output of
#'   [totalCoverage()]
#' @param gcCompensation GC compensation vector. Not support yet.
#' @param mappabilityCompensation mappability compensation vector. Not support
#'   yet.
#' @param FFT Use FFT smooth or not.
#' @param fft.sm.power the cut-off frequency of FFT smooth.
#'
#' @return a list. level 1: chromosome; level 2: each transcripts; level3: data
#'   matrix
#' @keywords internal
#' @import GenomicRanges
#' @importFrom IRanges Views
#'
UTR3TotalCoverage <- function(utr3, totalCov,
                              gcCompensation = NA,
                              mappabilityCompensation = NA,
                              FFT = FALSE,
                              fft.sm.power = 20) {
  utr3.s <- split(utr3, as.character(seqnames(utr3)))
  seqnames <- unlist(sapply(totalCov, names))
  seqnames <- table(seqnames)
  seqnames <- names(seqnames)[seqnames == length(totalCov)]
  seqnames <- sort(intersect(names(utr3.s), seqnames))
  cov <- list()
  for (seq in seqnames) {
    cov[[seq]] <- list()
    for (n in names(totalCov)) {
      cov[[seq]][[n]] <- totalCov[[n]][[seq]]
    }
  }
  mapply(function(.utr, .cvg, gcComp, mapComp) { ## memory consuming
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
  }, utr3.s[seqnames], cov[seqnames],
  gcCompensation[seqnames],
  mappabilityCompensation[seqnames],
  SIMPLIFY = FALSE
  )
}
