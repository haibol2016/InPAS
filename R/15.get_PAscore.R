#' calculate the CP score
#'
#' calculate the CP score by using PWM of polyadenylation signal
#'
#' @param seqname a character(n) vector, the chromosome/scaffold' name
#' @param pos genomic positions
#' @param str DNA strand
#' @param idx offset position
#' @param PWM An R object for a position weight matrix (PWM) for a hexamer
#'   polyadenylation signal (PAS), such as AAUAAA.
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param ups the number of upstream bases for PAS search.
#' @param dws the number of downstream bases for PAS search.
#' @param mc.cores integer(1), number of cores for the mc*apply function of the 
#'   parallel package
#' @return A list containing offset positions after PA score-based filtering
#' @import GenomicRanges
#' @importFrom BSgenome getSeq matchPWM
#' @seealso [get_PAscore2()]
#' @keywords internal
#' @author Jianhong Ou

get_PAscore <- function(seqname,
                        pos, 
                        str, 
                        idx, 
                        PWM, 
                        genome, 
                        ups = 50, 
                        dws = 50,
                        mc.cores = 1){
  pos <- pos[!is.na(pos)]
  if (length(pos) < 1) {
    return(NULL)
  }
  start <- pos - ups
  start[start < 1] <- 1
  end <- pos + dws
  gr <- GRanges(seqname, IRanges(start, end, 
                                 names = as.character(pos)),
                strand = str)
  seq <- getSeq(genome, gr)
  
  if (.Platform$OS.type == "windows" || mc.cores == 1){
    mT <- lapply(seq, matchPWM, pwm = PWM, 
                 min.score = "70%", with.score = TRUE)
  } else {
    mT  <- mclapply(seq, matchPWM, pwm = PWM, 
                    min.score = "70%", with.score = TRUE,
                    mc.cores = mc.cores)
  }
  
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
