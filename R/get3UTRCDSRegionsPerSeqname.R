#' Get 3' UTRs and their last CDS regions based on CP sites
#'
#' @param CPsites CP sites on a single chromosome/scaffold, as detected by
#'   [getCPsitesPerSeqname()]
#' @param utr3 output of [utr3Annotation()]
#' @param seqname character(1), chromosome/scaffold name
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply.
#'
#' @return An object of [GenomicRanges::GRangesList-class], with a GRanges for
#'   a seqname (chromosome/scaffold)
#' @export

get3UTRCDSRegionsPerSeqname <- function(CPsites,
                           utr3,
                           BPPARAM = NULL){
    stopifnot(length(CPsites) > 0 && is(CPsites, "GRanges"))
    
    if (!is(utr3, "GRanges") ||
        !all(utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
        stop("utr3 must be output of function of utr3Annotation")
    }
    
    utr3.trans.shorten.UTR <- split(CPsites, CPsites$transcript)
    
    if (!is.null(BPPARAM)) {
        utr3.regions <- bptry(bplapply(utr3.trans.shorten.UTR,
                                       getUTR3region,
                                       BPPARAM = BPPARAM
        ))
        while (!all(bpok(utr3.regions))) {
            utr3.regions <- bptry(bplapply(utr3.trans.shorten.UTR,
                                           getUTR3region,
                                           BPREDO = utr3.regions,
                                           BPPARAM = BPPARAM
            ))
        }
    } else {
        utr3.regions <- lapply(utr3.trans.shorten.UTR, getUTR3region)
    }
    utr3.regions <- unlist(GRangesList(utr3.regions))
    
    seqname <- unique(as.character(seqnames(utr3.regions)))
    
    ## merge utr3.regions with CDS of the same set of transcripts whose 
    ## utr3.regions are considered
    CDS <- utr3 %>% plyranges::filter(seqnames == seqname & 
                                      feature == "CDS" &
                        transcript %in% unique(utr3.regions$transcript))
    utr3.cds.regions <- c(utr3.regions, CDS)

    utr3.cds.regions.chr <- split(utr3.cds.regions, 
                                  as.character(seqnames(utr3.cds.regions)))
    return(utr3.cds.regions.chr)
}