#' Filter sequence names from a BSgenome object
#'
#' Filter sequence names for scaffolds from a BSgenome object so that only 
#'   chromosome-level seqnames are kept.
#'
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param removeScaffolds A logical(1) vector, whether the scaffolds should be
#'   removed from the genome. If you use a TxDb containing alternative
#'   scaffolds, you'd better to remove the scaffolds. To make things easy, we 
#'   suggest users creating a [BSgenome::BSgenome-class] instance from the 
#'   reference genome used for read alignment. For details, see the 
#'   documentation of [BSgenome::forgeBSgenomeDataPkg()].
#'   
#' @return An character vector containing filtered seqnames
#' @author Jianhong Ou, Haibo Liu
#' @importFrom BSgenome getSeq matchPWM
#' @keywords internal

trim_seqnames <- function(genome, removeScaffolds = FALSE) 
{
  if (!is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome")
  } 
  seqnames <- seqnames(genome)
  if (removeScaffolds) {
    seqnames <- seqnames[grepl("^(chr)?(\\d+|[XY])$", seqnames)]
  }
  seqnames
}
