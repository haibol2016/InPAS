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
#' @importFrom BSgenome getSeq matchPWM
#' @keywords internal
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' InPAS:::gcContents(genome, "chr1")
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
#' @param seqnames chromosome/scaffold names in the same forms of seqnames in
#'   the genome
#' @param window size of a sliding window, which optimally is set to the read
#'   length
#' @importFrom Biostrings alphabetFrequency toString
#'   letterFrequencyInSlidingView DNAString
#' @return A list of numeric vectors containing the weight (scaffold-level GC\%
#'   / GC\% per sliding window) for GC composition-based correction (Cheung et
#'   al. 2011).
#'
#' @keywords internal
#' @references Cheung MS, Down TA, Latorre I, Ahringer J. Systematic bias in
#'   high-throughput sequencing data and its correction by BEADS. Nucleic Acids
#'   Res. 2011 Aug;39(15):e103. doi: 10.1093/nar/gkr425. Epub 2011 Jun 6. PubMed
#'   PMID: 21646344;
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' genome <- BSgenome.Mmusculus.UCSC.mm10
#' InPAS:::gcComp(genome, "chr1")
gcComp <- function(genome, seqnames, window = 50) {
  fref <- sapply(
    seqnames,
    function(seqname) {
      gcContents(genome, seqname)
    }
  )

  win_size <- floor(window / 2)

  fdat <- lapply(seqnames, function(seqname) {
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
  })
}

#' Calculate weights for mappability-base coverage correction
#'
#' Calculate weights for mappability-base coverage correction
#'
#' @param mi A numeric vector of mappability along per
#' chromosome/scaffold
#'
#' @return A numeric vector of weights for mappability-based
#' coverage correction
#' @keywords internal
#'

# mappability is calculated by
##      [GEM](http://algorithms.cnag.cat/wiki/Man:gem-mappability)
## ref: Derrien T, Estellé J, Marco Sola S, Knowles DG, Raineri E,
##      Guigó R, Ribeca P.
##      Fast computation and applications of genome mappability. PLoS One.
##      2012;7(1):e30377. doi: 10.1371/journal.pone.0030377. Epub 2012 Jan 19.
##      PubMed PMID: 22276185; PubMed Central PMCID: PMC3261895.
## keywords internal PATH=$PATH:~/bin/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
## ./gem-indexer -i \
## genome.fa \
## -o mm10.index.gem
##
## ./gem-mappability -I mm10.index.gem.gem -l 100 -o mm10.mappability
## ./gem-2-wig -I mm10.index.gem.gem -i mm10.mappability -o mm10.mappability.wig

mapComp <- function(mi) {
  max(mi) / mi
}
