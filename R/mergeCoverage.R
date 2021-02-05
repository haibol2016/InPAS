#' merge read coverage from different samples
#'
#' @param coverage coverage for each sample, output of [coverageFromBedGraph()]
#' @param groupList a list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be removed
#'   from the genome
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply
#'
#' @return A list, coverage for a group of samples
#' @import BiocParallel
#' @keywords internal

mergeCoverage <- function(coverage, 
                          groupList,
                          genome, 
                          removeScaffolds = FALSE,
                          BPPARAM = NULL) {
  seqnames <- trimSeqnames(genome, removeScaffolds = removeScaffolds)
  seqLen <- seqLen(genome, removeScaffolds = removeScaffolds)
  seqRle <- list()
  for (i in seq_along(seqnames)) {
    seqRle[[seqnames[i]]] <- Rle(0, seqLen[i])
  }

  gpl <- rep(names(groupList), sapply(groupList, length))
  names(gpl) <- unlist(groupList)
  coverageL <- split(coverage[names(gpl)], gpl)

  cov <- list()

  for (j in 1:length(coverageL)) {
    covFiles <- coverageL[[j]]
    x <- 1:length(covFiles)
    y <- split(x, ceiling(x / 10))

    cov[[j]] <- seqRle

    for (i in 1:length(y)) {
      if (!is.null(BPPARAM)) {
        cv <- bptry(bplapply(covFiles[y[[i]]], function(.ele) {
          cvg <- NULL
          load(.ele)
          cvg <- cvg[seqnames]
          names(cvg) <- seqnames
          idx <- sapply(cvg, is.null)
          cvg[idx] <- seqRle[idx]
          cvg
        }, BPPARAM = BPPARAM))
        while (!all(bpok(cv))) {
          cv <- bptry(bplapply(covFiles[y[[i]]], function(.ele) {
            cvg <- NULL
            load(.ele)
            cvg <- cvg[seqnames]
            idx <- sapply(cvg, is.null)
            cvg[idx] <- seqRle[idx]
            cvg
          }, BPREDO = cv, BPPARAM = BPPARAM))
        }
        if (length(cv) < 10) {
          for (m in (length(cv) + 1):10) {
            cv[[m]] <- seqRle
          }
        }
        cov[[j]] <- bptry(bplapply(seqnames, function(s) {
          cv[[1]][[s]] + cv[[2]][[s]] + cv[[3]][[s]] +
            cv[[4]][[s]] + cv[[5]][[s]] + cv[[6]][[s]] +
            cv[[7]][[s]] + cv[[8]][[s]] + cv[[9]][[s]] +
            cv[[10]][[s]] + cov[[j]][[s]]
        }, BPPARAM = BPPARAM))
        while (!all(bpok(cov[[j]]))) {
          cov[[j]] <- bptry(bplapply(seqnames, function(s) {
            cv[[1]][[s]] + cv[[2]][[s]] + cv[[3]][[s]] +
              cv[[4]][[s]] + cv[[5]][[s]] + cv[[6]][[s]] +
              cv[[7]][[s]] + cv[[8]][[s]] + cv[[9]][[s]] +
              cv[[10]][[s]] + cov[[j]][[s]]
          }, BPREDO = cov[[j]], BPPARAM = BPPARAM))
        }
      } else {
        cv <- lapply(covFiles[y[[i]]], function(.ele) {
          cvg <- NULL
          load(.ele)
          cvg <- cvg[seqnames]
          names(cvg) <- seqnames
          idx <- sapply(cvg, is.null)
          cvg[idx] <- seqRle[idx]
          cvg
        })
        if (length(cv) < 10) {
          for (m in (length(cv) + 1):10) {
            cv[[m]] <- seqRle
          }
        }
        cov[[j]] <- lapply(seqnames, function(s) {
          cv[[1]][[s]] + cv[[2]][[s]] + cv[[3]][[s]] +
            cv[[4]][[s]] + cv[[5]][[s]] + cv[[6]][[s]] +
            cv[[7]][[s]] + cv[[8]][[s]] + cv[[9]][[s]] +
            cv[[10]][[s]] + cov[[j]][[s]]
        })
      }
      names(cov[[j]]) <- seqnames
    }
  }
  names(cov) <- names(coverageL)
  cov
}
