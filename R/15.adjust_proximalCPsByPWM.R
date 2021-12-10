#' adjust the proximal CP sites by matching PWM
#'
#' adjust the proximal CP sites by polyA Position Weight Matrix. It only need
#' the PWM to get match in upstream or downstream shift_range nr.
#'
#' @param idx  the offset of positions of CP sites
#' @param PolyA_PWM polyA PWM
#' @param seqnames a character(n) vector, the chromosome/scaffolds' names
#' @param starts start position in the genome
#' @param strands strands
#' @param genome an [BSgenome::BSgenome-class] object
#' @param shift_range  the shift range of PWM hits
#' @param search_point_START  Not use
#' @param mc.cores An integer(1) vector, number of cores for the mc*apply 
#'   function of the parallel package
#' @return the offset of positions of CP sites after filter
#' @details the hits is searched by [Biostrings::matchPWM()] and the cutoff is
#'   70\%
#' @seealso [adjust_proximalCPsByNBC()], [get_PAscore()]
#' @keywords internal
#' @author Jianhong Ou

adjust_proximalCPsByPWM <- function(idx, 
                             PolyA_PWM,
                             seqnames, 
                             starts,
                             strands,
                             genome,
                             shift_range,
                             search_point_START,
                             mc.cores = 1) {
  mapply(function(id, seqname, start, strand) {
    if (length(id) == 1) {
      if (is.na(id)) {
        return(NULL)
      }
    }
    if (length(id) > 0) {
      pos <- if (strand == "+") start + id - 1 else start - id + 1
      id <- InPAS:::get_PAscore(seqname, pos, strand,
                        id,
                        PWM = PolyA_PWM, 
                        genome = genome,
                        ups = shift_range + 25,
                        dws = shift_range + 25,
                        mc.cores)
    }
    id
  }, idx, seqnames, starts, strands, SIMPLIFY = FALSE)
}
