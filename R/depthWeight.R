#' calculate the depth weight for each example
#'
#' calculate the depth weight for each example: depth/mean(depth)
#'
#' @param coverage a list. output of [coverageFromBedGraph()]
#' @param hugeData is it a huge dataset?
#' @param groupList a list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#'
#' @return  a numeric vector with depth weight
#' @import S4Vectors
#' @keywords internal
#'

depthWeight <- function(coverage, hugeData, groupList = NULL) {
  n <- names(coverage)
  if (hugeData) {
    depth <- numeric(length(coverage))
    for (i in 1:length(coverage)) {
      cvg <- NULL
      load(coverage[[i]]) 
      
      d <- sapply(cvg, function(.cvg) {
        sum(as.double(runValue(.cvg)) * runLength(.cvg))
      })
      depth[i] <- sum(d)
      rm(cvg)
    }
    if (!is.null(names(groupList))[1]) {
      names(depth) <- names(coverage)
      groups <- rep(names(groupList), sapply(groupList, length))
      names(groups) <- unlist(groupList)
      groups <- groups[match(names(depth), names(groups))]
      if (any(is.na(groups))) {
        stop("all the tags in groupList must have correct names as it in coverage")
      }
      depth <- split(depth, groups)
      
      ## group mean depth
      depth <- sapply(depth, mean)
      n <- names(depth)
    }
  } else {
    depth <- sapply(coverage, function(cvg) {
      d <- sapply(cvg, function(.cvg) {
        sum(as.double(runValue(.cvg)) * runLength(.cvg))
      })
      sum(d)
    })
    depth.weight <- depth / mean(depth)
    names(depth.weight) <- n
    depth.weight
  }
}

