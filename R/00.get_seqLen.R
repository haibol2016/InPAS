#' Get sequence lengths for chromosomes/scaffolds
#'
#' Get sequence lengths for chromosomes/scaffolds from a 
#'   [BSgenome::BSgenome-class] object
#' 
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or 
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#'
#' @return A named numeric vector containing lengths per seqname, with the
#'   seqnames as the names
#' @seealso [GenomeInfoDb::Seqinfo-class]
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' InPAS:::get_seqLen(genome = genome, 
#'                    chr2exclude = "chrM")

get_seqLen <- function(genome = getInPASGenome(), 
                       chr2exclude = getChr2Exclude()) {
  if (!is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome class")
  }
  if (!is.null(chr2exclude) && !is.character(chr2exclude))
  {
    stop("chr2Exclude must be NULL or a character vector")
  }
  len <- seqlengths(genome)
  seqnames <- trim_seqnames(genome, chr2exclude)
  len <- len[seqnames]
  len
}
