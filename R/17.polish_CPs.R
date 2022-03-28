#' polish the searching results of CP sites
#'
#' remove the multiple positions of CP sites for the same 3' UTRs and only keep
#' the best CP sites for proximal and distal.
#'
#' @param CPs output of [search_proximalCPs()] or [adjust_proximalCPs()]
#' @param output.all A logical(1), indicating whether to output entries with only 
#'   single CP site for a 3' UTR.
#' @param DIST2END An integer(1) vector, specifying minimal length 
#' difference between proximal and distal APA sites which should be met to be 
#' considered for outputted if *output.all* is set to TRUE. Default is 200 bp.
#'
#' @return a data.frame with columns: "fit_value", "Predicted_Proximal_APA",
#'   "Predicted_Distal_APA", "utr3start", "utr3end", "Predicted_Distal_APA_type"
#' @keywords internal
#' @seealso [adjust_proximalCPs()], [adjust_proximalCPsByPWM()],
#'   [adjust_proximalCPsByNBC()], [get_PAscore2()]
#' @author Jianhong Ou

polish_CPs <- function(CPs, output.all, DIST2END = 200) {
  dCPs <- CPs$dCPs
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
  
  dCPs$min_fit_value <- unlist(CPs$fit_value_min)
  CPs$Predicted_Proximal_APA[sapply(
    CPs$Predicted_Proximal_APA,
    length
  ) == 0] <- NA
  dCPs$Predicted_Proximal_APA <- 
      unlist(ifelse(is.na(dCPs$Predicted_Proximal_APA),
             unlist(CPs$Predicted_Proximal_APA),
             dCPs$Predicted_Proximal_APA))
  dCPs$Predicted_Proximal_APA_Type <- 
      ifelse(!is.na(dCPs$Predicted_Proximal_APA) & 
                 is.na(dCPs$Predicted_Proximal_APA_Type),
             "Novel proximal APA",
             dCPs$Predicted_Proximal_APA_Type)
  dCPs$annotatedProximalCP <- unlist(dCPs$annotatedProximalCP)
  
  ## remove entries with no detected Predicted_Proximal_APA and 
  ## Predicted_Proximal_APA very close to Predicted_Distal_APA
  if (!output.all) 
  {
    dCPs <- dCPs[!is.na(dCPs$Predicted_Proximal_APA) & 
                   abs(dCPs$Predicted_Proximal_APA - 
                         dCPs$Predicted_Distal_APA) >= DIST2END, ]
  }
  dCPs
}
