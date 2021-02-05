#' extract coverage of last CDS exon region
#'
#' extract coverage of last CDS exon region
#'
#' @param CDS a [GenomicRanges::GRanges-class] object for CDS regions
#' @param coverage output of [coverageFromBedGraph()]
#' @param hugeData logical(1), indicating whether the input is a huge dataset
#' @param BPPARAM An optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply.
#' @param phmm logical(1), indicating whether output data is prepared for
#'   singleSample analysis
#'
#' @return the average coverage of last CDS for each transcript
#' @import GenomicRanges BiocParallel
#' @keywords internal

lastCDSusage <- function(CDS, coverage, 
                         hugeData,
                         BPPARAM = NULL, 
                         phmm = FALSE) {
  
  CDS.regions.chr <- split(CDS, as.character(seqnames(CDS)))
  seqnames <- names(CDS.regions.chr)
  if (!is.null(BPPARAM)) {
    CDS.regions.chr <- bplapply(seqnames, getRegionCoverage,
      BPPARAM = BPPARAM,
      utr3.regions.chr = CDS.regions.chr,
      hugeData = hugeData,
      coverage = coverage,
      phmm = phmm
    )
  } else {
    CDS.regions.chr <- lapply(seqnames, getRegionCoverage,
      utr3.regions.chr = CDS.regions.chr,
      hugeData = hugeData,
      coverage = coverage,
      phmm = phmm
    )
  }
  cov.ranges <- unlist(GRangesList(CDS.regions.chr))
}
