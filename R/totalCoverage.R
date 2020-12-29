#' calculate the total coverage
#'
#' For huge dataset, it will read in the coverage from temp files and merge them
#' by groups
#'
#' @param coverage coverage for each sample, outputs of [coverageFromBedGraph()]
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be removed
#'   from the genome
#' @param hugeData hugeData or not
#' @param groupList a list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#'
#' @return a coverage list
#' @keywords internal
#' @import S4Vectors

totalCoverage <- function(coverage, genome,
                          hugeData = FALSE,
                          removeScaffolds = FALSE,
                          groupList = NULL) {
  if (!hugeData) {
    return(coverage)
  }
  seqnames <- trimSeqnames(genome, removeScaffolds)
  seqLen <- seqLen(genome, removeScaffolds)

  if (!is.null(names(groupList))[1]) {
    cov <- vector("list", length = length(groupList))
    names(cov) <- names(groupList)
    for (i in 1:length(cov)) {
      cov[[i]] <- list()
      for (s in seqnames) {
        cov[[i]][[s]] <- Rle(0, seqLen[s])
      }
    }

    for (i in 1:length(coverage)) {
      cvg <- NULL
      load(coverage[[i]])
      gp <- names(coverage)[i]
      gp <- which(sapply(groupList, function(.ele) any(gp %in% .ele)))[1]
      if (!is.null(gp)) {
        for (s in seqnames) {
          if (s %in% names(cvg)) {
            cov[[gp]][[s]] <- cov[[gp]][[s]] + cvg[[s]]
          }
        }
      }
      rm(cvg)
    }
    idx <- rep(FALSE, length(seqnames))
    names(idx) <- seqnames
    for (i in 1:length(cov)) {
      for (s in seqnames) {
        if (nrun(cov[[i]][[s]]) != 1 || runValue(cov[[i]][[s]][1]) != 0) {
          idx[s] <- TRUE
        }
      }
    }
    for (i in 1:length(cov)) {
      cov[[i]] <- cov[[i]][idx]
    }
  } else {
    cov <- vector("list", length = length(coverage))
    names(cov) <- names(coverage)
    for (i in 1:length(coverage)) {
      cvg <- NULL
      load(coverage[[i]])
      gp <- names(coverage)[i]
      cov[[gp]] <- list()
      for (s in seqnames) {
        if (s %in% names(cvg)) {
          cov[[gp]][[s]] <- cvg[[s]]
        }
      }
      rm(cvg)
    }
  }
  cov
}
