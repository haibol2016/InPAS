#' Filter sequence names from a BSgenome object
#'
#' Filter sequence names for scaffolds from a BSgenome object so that only 
#'   chromosome-level seqnames are kept.
#'
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or 
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#'   
#' @return An character vector containing filtered seqnames
#' @author Jianhong Ou, Haibo Liu
#' @keywords internal

trim_seqnames <- function(genome = getInPASGenome(), 
                          chr2exclude = getChr2Exclude()) 
{
  if (!is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome")
  }
  if (!is.null(chr2exclude) && !is.character(chr2exclude))
  {
    stop("chr2exclude must be NULL or a character vector")
  }
  seqnames <- seqnames(genome)
  seqnames <- seqnames[!seqnames %in% chr2exclude]
  seqnames
}
