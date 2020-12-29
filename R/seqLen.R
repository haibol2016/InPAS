#' Get sequence lengths
#'
#' Get sequence lengths from a [BSgenome::BSgenome-class] object
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be removed
#'   from the genome
#'
#' @return A numeric vector containing lengths per seqnames
#' @seealso [GenomeInfoDb::Seqinfo-class]
#' @importFrom BSgenome getSeq matchPWM
#' @keywords internal

seqLen <- function(genome, removeScaffolds = FALSE) {
  if (!is(genome, "BSgenome")) stop("genome must be an object of BSgenome")
  len <- seqlengths(genome)
  if (removeScaffolds) {
    seqnames <- trimSeqnames(genome, removeScaffolds = removeScaffolds)
    len <- len[seqnames]
  }
  len
}
