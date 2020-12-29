#' calculate the CP score
#'
#' calculate the CP score by using PWM of polyadenylation signal
#'
#' @param seqname sequence names
#' @param pos genomic positions
#' @param str strands
#' @param idx offset position
#' @param PWM polyA position weight matrix
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param ups upstream base
#' @param dws downstream base
#'
#' @return A list containing offset positions after filtering
#' @import GenomicRanges
#' @importFrom BSgenome getSeq matchPWM
#' @importClassesFrom IRanges IRanges
#' @seealso [PAscore2()]
#' @keywords internal
#'

PAscore <- function(seqname, pos, str, idx, PWM, genome, ups = 50, dws = 50) {
  pos <- pos[!is.na(pos)]
  if (length(pos) < 1) {
    return(NULL)
  }
  start <- pos - ups
  start[start < 1] <- 1
  end <- pos + dws
  gr <- GRanges(seqname, IRanges(start, end, names = as.character(pos)),
    strand = str
  )
  seq <- getSeq(genome, gr)
  mT <- lapply(seq, matchPWM, pwm = PWM, min.score = "70%", with.score = TRUE)
  hits <- sapply(mT, function(.ele) {
    if (!is(.ele, "XStringViews")) {
      return(FALSE)
    }
    if (length(.ele) == 0) {
      return(FALSE)
    }
    TRUE
  })
  idx[hits]
}
