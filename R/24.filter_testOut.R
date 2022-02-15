#' filter 3' UTR usage test results
#'
#' filter results of [test_dPDUI()]
#'
#' @param res a [UTR3eSet-class] object, output of [test_dPDUI()]
#' @param gp1 tag names involved in group 1. gp1 and gp2 are used for filtering
#'   purpose if both are specified; otherwise only other specified thresholds
#'   are used for filtering.
#' @param gp2 tag names involved in group 2
#' @param outdir A character(1) vector, a path with write permission for storing
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param background_coverage_threshold  background coverage cut off value. for
#'   each group, more than half of the long form should greater than
#'   background_coverage_threshold. for both group, at least in one group, more
#'   than half of the short form should greater than
#'   background_coverage_threshold.
#' @param P.Value_cutoff cutoff of P value
#' @param adj.P.Val_cutoff  cutoff of adjust P value
#' @param dPDUI_cutoff  cutoff of dPDUI
#' @param PDUI_logFC_cutoff cutoff of PDUI log2 transformed fold change
#'
#' @return A data frame converted from an object of
#'    [GenomicRanges::GRanges-class].
#' @seealso [test_dPDUI()]
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' library(limma)
#' path <- system.file("extdata", package = "InPAS")
#' load(file.path(path, "eset.MAQC.rda"))
#' tags <- colnames(eset@PDUI)
#' g <- factor(gsub("\\..*$", "", tags))
#' design <- model.matrix(~ -1 + g)
#' colnames(design) <- c("Brain", "UHR")
#' contrast.matrix <- makeContrasts(
#'   contrasts = "Brain-UHR",
#'   levels = design
#' )
#' res <- test_dPDUI(
#'   eset = eset,
#'   method = "limma",
#'   normalize = "none",
#'   design = design,
#'   contrast.matrix = contrast.matrix
#' )
#' filter_testOut(res,
#'   gp1 = c("Brain.auto", "Brain.phiX"),
#'   gp2 = c("UHR.auto", "UHR.phiX"),
#'   background_coverage_threshold = 2,
#'   P.Value_cutoff = 0.05,
#'   adj.P.Val_cutoff = 0.05,
#'   dPDUI_cutoff = 0.3,
#'   PDUI_logFC_cutoff = .59
#' )
filter_testOut <- function(res,
                           gp1,
                           gp2,
                           outdir = getInPASOutputDirectory(),
                           background_coverage_threshold = 2,
                           P.Value_cutoff = 0.05,
                           adj.P.Val_cutoff = 0.05,
                           dPDUI_cutoff = 0.2,
                           PDUI_logFC_cutoff = log2(1.5)) {
  if (missing(res) || !is(res, "UTR3eSet")) {
    stop("res is required and must be an object of UTR3eSet")
  }
  if (!is.character(outdir) || length(outdir) != 1) {
    stop("An explicit output directory is required")
  } else {
    outdir <- file.path(outdir, "010.filtered.dPDUI")
    if (!dir.exists(outdir)) {
      dir.create(outdir,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }
    outdir <- normalizePath(outdir)
  }
  if (missing(gp1) || missing(gp2)) {
    testRes <- res$testRes
    testRes <- testRes[!(is.na(testRes[, "P.Value"]) |
      is.na(testRes[, "adj.P.Val"])), ]
    if (!"dPDUI" %in% colnames(testRes)) {
      message("dPDUI_cutoff will be ignored")
    } else {
      testRes <- testRes[abs(testRes[, "dPDUI"]) >=
        dPDUI_cutoff, , drop = FALSE]
    }
    if (!missing(PDUI_logFC_cutoff)) {
      if (!"logFC" %in% colnames(testRes)) {
        message("PDUI_logFC_cutoff will be ignored")
      } else {
        testRes <- testRes[abs(testRes[, "logFC"]) >=
          PDUI_logFC_cutoff, , drop = FALSE]
      }
    }
    testRes <- testRes[testRes[, "adj.P.Val"] <= adj.P.Val_cutoff &
      testRes[, "P.Value"] <= P.Value_cutoff, , drop = FALSE]
    testRes <- testRes[!is.na(rownames(testRes)), ]
    if (nrow(testRes) == 0) {
      return(NULL)
    }
    res <- res$usage[match(rownames(testRes), res$usage$transcript)]
    mcols(res) <- cbind(as.data.frame(mcols(res)), testRes)
    return(res)
  }

  ## In At least one group, more than half of the long/short form
  ## should be greater than background_coverage_threshold
  coln <- colnames(res$short)
  if (!all(c(gp1, gp2) %in% coln)) {
    stop("gp1 and gp2 must be the tags of samples (colnames of res$short)")
  }
  short <- res$short
  long <- res$long
  short.abg <- short > background_coverage_threshold
  long.abg <- long > background_coverage_threshold
  short.abg.gp1 <- short.abg[, gp1, drop = FALSE]
  short.abg.gp2 <- short.abg[, gp2, drop = FALSE]
  long.abg.gp1 <- long.abg[, gp1, drop = FALSE]
  long.abg.gp2 <- long.abg[, gp2, drop = FALSE]
  halfLgp1 <- ceiling(length(gp1) / 2)
  halfLgp2 <- ceiling(length(gp2) / 2)
  short.abg.gp1 <- rowSums(short.abg.gp1)
  short.abg.gp2 <- rowSums(short.abg.gp2)
  long.abg.gp1 <- rowSums(long.abg.gp1)
  long.abg.gp2 <- rowSums(long.abg.gp2)
  long.abg <- long.abg.gp1 >= halfLgp1 | long.abg.gp2 >= halfLgp2
  short.abg <- short.abg.gp1 >= halfLgp1 | short.abg.gp2 >= halfLgp2
  abg.gp1 <- long.abg.gp1 >= halfLgp1 | short.abg.gp1 >= halfLgp1
  abg.gp2 <- long.abg.gp2 >= halfLgp2 | short.abg.gp2 >= halfLgp2
  abg <- abg.gp1 & abg.gp2 & long.abg & short.abg

  PDUI.gp1 <- res$PDUI[, gp1, drop = FALSE]
  PDUI.gp2 <- res$PDUI[, gp2, drop = FALSE]
  PDUI.gp1 <- rowMeans(PDUI.gp1)
  PDUI.gp2 <- rowMeans(PDUI.gp2)
  if (!"dPDUI" %in% colnames(res$testRes)) {
    message("dPDUI is calculated by gp2 - gp1.")
    dPDUI <- PDUI.gp2 - PDUI.gp1
  } else {
    dPDUI <- res$testRes[, "dPDUI"]
  }
  adPDUI <- abs(dPDUI) >= dPDUI_cutoff

  if (!missing(PDUI_logFC_cutoff)) {
    if (!"logFC" %in% colnames(res$testRes)) {
      logFC <- log2(PDUI.gp2 + .Machine$double.xmin) -
        log2(PDUI.gp1 + .Machine$double.xmin)
    } else {
      logFC <- res$testRes[, "logFC"]
    }
    alFC <- abs(logFC) >= PDUI_logFC_cutoff
  } else {
    alFC <- rep(TRUE, length(abg))
  }

  res <- as(res, "GRanges")
  res$dPDUI <- dPDUI
  res$PASS <- res$adj.P.Val <= adj.P.Val_cutoff &
    res$P.Value <= P.Value_cutoff &
    abg & adPDUI & alFC
  res <- data.frame(sort(res))
  res <- res[, !colnames(res) %in% c(
    "width", "annotatedProximalCP",
    "truncated", "fit_value",
    "feature", "exon"
  )]
  saveRDS(res, file = file.path(outdir, "filtered.dPDUI.RDS"))
  res
}
