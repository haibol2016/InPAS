#' polish the searching results of CP sites
#'
#' remove the multiple positions of CP sites for the same 3' UTRs and only keep
#' the best CP sites for proximal and distal.
#'
#' @param CPs output of [search_proximalCPs()] or [adjust_proximalCPs()]
#'
#' @return a matrix with columns: "fit_value", "Predicted_Proximal_APA",
#'   "Predicted_Distal_APA", "utr3start", "utr3end", "distalCPtype"
#' @keywords internal
#' @seealso [adjust_proximalCPs()], [adjust_proximalCPsByPWM()],
#'   [adjust_proximalCPsByNBC()], [get_PAscore2()]
#' @author Jianhong Ou

polish_CPs <- function(CPs) {
  coors <- lapply(
    CPs$chr.cov.merge,
    function(.ele) {
      as.numeric(gsub(".*_SEP_", "", rownames(.ele)))
    }
  )
  flag <- CPs$flag

  ## make a copy of fit_value and predicted_proximal_APA for visualization if
  ## necessary
  CPs$Predicted_Proximal_APA <- mapply(function(.ele, .id) {
    ifelse(is.na(.id[1]) | is.null(.id[1]), NA, .ele[.id[1]])
  }, coors, CPs$Predicted_Proximal_APA, SIMPLIFY = FALSE)
  
  CPs$dCPs$fit_value <- unlist(CPs$fit_value_min)
  CPs$Predicted_Proximal_APA[sapply(
    CPs$Predicted_Proximal_APA,
    length
  ) == 0] <- NA
  CPs$dCPs$Predicted_Proximal_APA[flag] <- 
    unlist(CPs$Predicted_Proximal_APA[flag])
  CPs$dCPs
}
