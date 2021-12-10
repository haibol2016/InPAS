#' Calculate the depth weight for each sample or each experimental condition
#'
#' Calculate the depth weight for each sample of non-hugeData or each 
#'   experimental condition for hugeData: depth/mean(depth)
#'
#' @param metadata A data frame containing the metadata for a RNA-seq experiment,
#'  which can be extract from the SQLite database set up by [setup_sqlitedb()]
#' @param hugeData A logical(1), indicating whether it is huge data
#' @return A named numeric vector containing depth weight for each sample for 
#'   non-hugeData, or depth weight for each condition if hugeData.
#' @import S4Vectors RSQLite
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu

get_depthWeight <- function(metadata, hugeData) {
  tags <- metadata$tag
  depth <- metadata$depth
  condition <- metadata$condition
  
  if (hugeData){
    if (length(unique(condition)) > 1) {
      depth <- split(depth, condition)
      ## group mean depth
      depth <- sapply(depth, mean)
      depth.weight <- depth / mean(depth)
      names(depth.weight) <- names(depth)
    } else {
      depth.weight <- 1
      names(depth.weight) <- unique(condition)
    }
  } else if (length(unique(condition)) == 1) {
    depth.weight <- depth / mean(depth)
    names(depth.weight) <- tags
  }
  depth.weight
}

