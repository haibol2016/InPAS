#' helper function to calculate chromosome/scaffold level GC content
#'
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param seqname a character(1) vector, the chromosome/scaffold's name
#' @param nonATCGExclude a logical(1) vector, whether nucleotides other than A,
#'   T, C, and G should be excluded when GC content is calculated
#'
#' @return a numeric(1) vector, containing the chromosome/scaffold -specific GC
#'   content in the range of 0 to 1
#' @importFrom Biostrings alphabetFrequency
#' @keywords internal
#' @author Haibo Liu
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' InPAS:::gcContents(genome, "chr1")
#' }

gcContents <- function(genome, seqname, nonATCGExclude = TRUE) {
  if (!is(genome, "BSgenome")) {
    stop("Genome must be a BSgenome object!")
  }
  if (!seqname %in% seqnames(genome)) {
    stop("seqname doesn't appear in the genome!")
  }
  freq <- alphabetFrequency(genome[[seqname]], as.prob = FALSE)

  if (nonATCGExclude) {
    gc <- sum(freq[c("G", "C")]) / sum(freq[c("A", "T", "G", "C")])
  } else {
    gc <- sum(freq[c("G", "C")]) / sum(freq)
  }
  gc
}

#' Calculate weights for GC composition
#'
#' Calculate read weights for GC composition-based coverage correction
#'
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param seqnames a character(n) vector, the chromosome/scaffolds' names
#'  in the same forms of seqnames in the BSgenome
#' @param window size of a sliding window, which optimally is set to the read
#'   length
#' @param future.chunk.size The average number of elements per future 
#'   ("chunk"). If Inf, then all elements are processed in a single future.
#'   If NULL, then argument future.scheduling = 1 is used by default. Users can
#'   set future.chunk.size = total number of elements/number of cores set for 
#'   the backend. See the future.apply package for details.
#' @importFrom Biostrings alphabetFrequency toString
#'   letterFrequencyInSlidingView DNAString
#' @return A list of numeric vectors containing the weight (scaffold-level GC\%
#'   / GC\% per sliding window) for GC composition-based correction for each 
#'   chromosome/scaffold.
#' @importFrom future.apply future_lapply
#' @author Jianhong Ou, Haibo Liu
#' @keywords internal
#' @references Cheung et al. Systematic bias in high-throughput sequencing data
#'   and its correction by BEADS. Nucleic Acids Res. 2011 Aug;39(15):e103.
#' @examples
#' \dontrun{
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' InPAS:::gcComp(genome, "chr1")
#' }

gcComp <- function(genome, seqnames, window = 50, 
                   future.chunk.size = NULL) {
  fref <- sapply(
    seqnames,
    function(seqname) {
      gcContents(genome, seqname)
    }
  )

  win_size <- floor(window / 2)

  ## helper function
  window_gc <- function(seqname){
    dna <- toString(genome[[seqname]])
    dna <- paste0(
      dna,
      paste(rep("N", win_size - 1), collapse = ""))
    window_gc <- 
      letterFrequencyInSlidingView(DNAString(dna),
                                   view.width = win_size,
                                   letters = "CG",
                                   OR = "|",
                                   as.prob = TRUE)
    
    ## seqname-level GC% / per-window GC%
    ## what if window_gc = 0?
    wk <- fref[seqname] / window_gc[, 1]
  }
    fdat <- future_lapply(seqnames, window_gc,
                          future.chunk.size = future.chunk.size,
                          future.stdout = NA)
    fdat
}

#' Calculate weights for mappability-base coverage correction
#'
#' Calculate weights for mappability-base coverage correction
#' 
#' @description mappability is calculated by using 
#' \href{http://algorithms.cnag.cat/wiki/Man:gem-mappability}{GEM} 
#'  with the following command lines: 
#' PATH=$PATH:~/bin/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
#' ./gem-indexer -i genome.fa -o mm10.index.gem
#' ./gem-mappability -I mm10.index.gem.gem -l 100 -o mm10.mappability
#' ./gem-2-wig -I mm10.index.gem.gem -i mm10.mappability -o mm10.mappability.wig
#' 
#' @param mi A numeric vector of mappability along per
#' chromosome/scaffold
#'
#' @return A numeric vector of weights for mappability-based
#' coverage correction
#' @keywords internal
#' @author Jianhong Ou
#' @references Derrien et al. Fast computation and applications of genome
#'   mappability. PLoS One. 2012;7(1):e30377. doi: 10.1371/journal.pone.0030377.
#' 

mapComp <- function(mi) {
  max(mi) / mi
}
