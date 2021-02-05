#' Get coverage for 3' UTR and last CDS regions on a single chromosome
#'
#' @param utr3.cds.regions.chr An object of [GenomicRanges::GRangesList-class],
#'   with a GRanges for a single seqname (chromosome/scaffold). It is an  
#'   output from the [get3UTRCDSRegionsPerSeqname()]
#' @param hugeData logical(1), indicating whether it is huge data
#' @param coverage coverage for each sample, output of [coverageFromBedGraph()],
#'   [getCov4SmallExperiment()], or [getCovFileList()] 
#' @param tmpfolder A path to a folder for storing coverage data of 3' UTRs on a
#'   given chromosome/scaffold
#' @param phmm logical(1), indicating whether data should be prepared for 
#'   singleSample analysis? By default, FALSE
#'
#' @return coverage view
#' @export
#' 
get3UTRCDSRegionCoveragePerSeqname <- function(utr3.cds.regions.chr,
                                  hugeData = TRUE,
                                  coverage,
                                  tmpfolder = NULL,
                                  phmm = FALSE){
  stopifnot(is(utr3.cds.regions.chr, "GRangesList"))
  if (hugeData && is.null(tmpfolder)){
      stop("A folder musted be specified to store coverage", 
      "data for 3' UTR when it is huge data")
  }
  stopifnot(length(coverage) > 0)
  
  seqname <- names(utr3.cds.regions.chr)
  
  if (length(seqname) != 1){
    stop("utr3.cds.regions.chr should only be for a single chromosome/scaffold")
  }
  
  utr3.regions.cov <- 
        getRegionCoverage(chr = seqname,
                          utr3.regions.chr = utr3.cds.regions.chr,
                          hugeData = hugeData,
                          coverage = coverage,
                          phmm = phmm)
  if (hugeData){
      saveRDS(utr3.regions.cov, 
              file = file.path(tmpfolder, 
                            paste(seqname, "utr3.coverage.RDS")))
      return(NULL)
  } else {
      return(utr3.regions.cov)
  }
}