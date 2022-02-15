#' Annotation of 3' UTRs for mouse (mm10)
#'
#' A dataset containing the annotation of the 3' UTRs of the mouse
#'
#' @format An object of [GenomicRanges::GRanges-class] with 7 metadata columns
#' \describe{
#'   \item{feature}{feature type, utr3, CDS, next.exon.gap}
#'   \item{annotatedProximalCP}{candidate proximal CPsites}
#'   \item{exon}{exon ID}
#'   \item{transcript}{transcript ID}
#'   \item{gene}{gene ID}
#'   \item{symbol}{gene symbol}
#'   \item{truncated}{whether the 3' UTR is trucated}
#' }
"utr3.mm10"
