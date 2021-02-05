
#' Assemble chromosome-wise 3'UTR coverage data into a GRanges object
#'
#' Assemble chromosome-wise 3'UTR coverage data for each sample into a GRanges
#'   object
#' 
#' @param UTR3CoverageFiles paths to files storing chromosome-wise 3'UTR 
#'   coverage data for each sample
#'
#' @return An object of [GenomicRanges::GRanges-class]
#' @export

integrate3UTRUsage <- function(UTR3CoverageFiles){
    
    if (!any(file.exists(UTR3CoverageFiles))){
        stop("Some files for UTR3Coverage don't exist")
    }
    utr.cvg <- lapply(UTR3CoverageFiles, function(.x){
        chr.utr3.cvg <- readRDS(.x)
    })
    utr.cvg <- unlist(GRangesList(utr.cvg))
    utr.cvg <- sort(utr.cvg)
}
    