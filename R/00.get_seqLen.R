#' Get sequence lengths for chromosomes/scaffolds
#'
#' Get sequence lengths for chromosomes/scaffolds from a 
#'   [BSgenome::BSgenome-class] object
#' 
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param removeScaffolds A logical(1) vector, whether the scaffolds should be
#'   removed from the genome. If you use a TxDb containing alternative
#'   scaffolds, you'd better to remove the scaffolds. To make things easy, we 
#'   suggest users creating a [BSgenome::BSgenome-class] instance from the 
#'   reference genome used for read alignment. For details, see the 
#'   documentation of [BSgenome::forgeBSgenomeDataPkg()].
#'
#' @return A named numeric vector containing lengths per seqname, with the
#'   seqnames as the names
#' @seealso [GenomeInfoDb::Seqinfo-class]
#' @importFrom BSgenome getSeq matchPWM
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' InPAS:::get_seqLen(genome, removeScaffolds = FALSE)

get_seqLen <- function(genome, removeScaffolds = FALSE) {
  if (missing(genome) || !is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome class")
  }
  len <- seqlengths(genome)
  if (removeScaffolds) {
    seqnames <- trim_seqnames(genome, removeScaffolds = removeScaffolds)
    len <- len[seqnames]
  }
  len
}
