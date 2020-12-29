#' calculate the usage of long and short forms of UTR3
#'
#' calculate the usage of long and short forms of UTR3 for the results of
#' [CPsites()]
#'
#' @param CPsites outputs of [CPsites()]
#' @param coverage coverage for each sample, outputs of [coverageFromBedGraph()]
#' @param hugeData is this dataset consume too much memory? if it is TRUE, the
#'   coverage will be saved into tempfiles.
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply.
#' @param phmm logical(1), should data be prepared for singleSample analysis? By
#'   default, FALSE
#'
#' @return A [GenomicRanges::GRanges-class] object
#' @seealso [CPsites()]
#' @keywords internal
#' @import GenomicRanges
#'

UTR3usage <- function(CPsites, coverage, hugeData,
                      BPPARAM = NULL, phmm = FALSE) {
  utr3.selected.shorten.UTR <-
    split(CPsites, as.character(seqnames(CPsites)))
  seqnames <- names(utr3.selected.shorten.UTR)

  ## convert coverage from list of Rle to list of matrix
  ## step1 prepare the short and long utr3 regions
  utr3.trans.shorten.UTR <- split(CPsites, CPsites$transcript)
  if (!is.null(BPPARAM)) {
    utr3.regions <- bplapply(utr3.trans.shorten.UTR,
      getUTR3region,
      BPPARAM = BPPARAM
    )
  } else {
    utr3.regions <- lapply(utr3.trans.shorten.UTR, getUTR3region)
  }
  utr3.regions <- unlist(GRangesList(utr3.regions))
  ## step2, calculate coverage and merge into a matrix for long and short
  utr3.regions.chr <- split(utr3.regions, as.character(seqnames(utr3.regions)))
  utr3.regions.chr <- utr3.regions.chr[seqnames]
  if (!is.null(BPPARAM)) {
    utr3.regions.chr <- bplapply(seqnames, getRegionCoverage,
      BPPARAM = BPPARAM,
      utr3.regions.chr = utr3.regions.chr,
      hugeData = hugeData,
      coverage = coverage,
      phmm = phmm
    )
  } else {
    utr3.regions.chr <- lapply(seqnames, getRegionCoverage,
      utr3.regions.chr = utr3.regions.chr,
      hugeData = hugeData,
      coverage = coverage,
      phmm = phmm
    )
  }
  cov.ranges <- unlist(GRangesList(utr3.regions.chr))
}
